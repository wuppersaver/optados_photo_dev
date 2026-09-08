# Photoemission: making the answer independent of the DFT setup

Implementation brief. Written 7 Sep 2026 from a read-through of `src/photo.f90` and
the modules it leans on (`optics`, `jdos_utils`, `dos_utils`, `cell`, `electronic`,
`parameters`), on branch `photoemission-fixes` at commit `12455c6`. Line numbers
below are from that commit; they will drift, so re-grep the quoted code before
editing rather than trusting the number.

The question asked was: where does changing something about the DFT setup that
should not matter (vacuum thickness, number of layers, k-point grid, whether the
k-list is symmetry reduced, number of bands) change the photoemission numbers?
Nine places were found, and a tenth later. They are listed worst first, each with the evidence, the
change, and how to verify it. Nothing here is implemented yet.

A tenth was added on 8 Sep 2026 (item 10), found while auditing item 1: the
symmetry star that every placing routine walks is the operation list, and CASTEP
does not write time reversal into it. It changes no `.odo` number, and it does
not arise on any structure in the current campaign -- but it changes item 1's
clean-up step, so read the two together.

Companion document: `Photoemission_outstanding_work.md` records the previous round
of fixes and the test-fixture plan. Read its section 4.2 and 4.3 before touching
the test suite; the traps recorded there (symlinked fixtures, benchmark file
names, `jobconfig` as the authority on what each test compares) all still apply.

---

## 0. Ground rules for whoever does this

**Binaries.** Four untracked binaries sit in `optados/`. The test suite's
`test-suite/tests/userconfig` runs `optados.x.x86_64` (serial),
`optados.mpi.x86_64` and `optados.debug.x86_64`. `optados.x` is a stale build that
rejects `task : photoemission`; ignore it. Rebuild all three configurations after
every change. `make.system` sets `BUILD := fast` and `COMMS_ARCH := mpi` with
`:=`, so command-line variables override them:

```
make optados BUILD=fast  COMMS_ARCH=serial   # optados.x.x86_64
make optados BUILD=fast  COMMS_ARCH=mpi      # optados.mpi.x86_64
make optados BUILD=debug COMMS_ARCH=serial   # optados.debug.x86_64
```

Confirm each binary's mtime is newer than `src/photo.o` before running anything.
`make.system` is currently modified in the working tree; do not commit it.

**Tests.** From `test-suite/`, `./run_tests` interactively, or non-interactively
with the `photoemission_only` category (`testopt_photo_*`, 34 tests). The
photoemission `.odo` comparison has a relative tolerance of 1e-6 on `qe_total`,
`qe_bulk` and `mte`, so every fix below that changes a number will fail its
benchmark until re-blessed. Re-bless by copying the file named on the test's
`output =` line in `tests/jobconfig` to `benchmark.out.default.inp=Cu100_5L.odi`
in that test directory. Seven tests compare a `.dat`, fifteen compare
`Cu100_5L.opt_err`; do not bless those from the `.odo`.

**Fixture.** `test-suite/checkpoints/Cu100_5L_photo/` is a 5-layer Cu(100) slab,
9x9x1 grid reduced to five irreducible k-points with weights 4, 4, 8, 4, 1 out of
21. The QE is dominated (77 %) by one k-point. It is deterministic, not physical;
use it to check that a change does what the arithmetic says, not to judge
whether a number is reasonable.

**Scratch runs.** To run the fixture by hand outside the harness:

```
mkdir work && cd work
for f in bands dome_fmt fem_fmt ome_fmt pdos_fmt tmprob_fmt gkgrid_fmt; do
  bunzip2 -kc <repo>/test-suite/checkpoints/Cu100_5L_photo/Cu100_5L.$f.bz2 > Cu100_5L.$f
done
for f in dome fem ome pdos tmprob gkgrid; do <repo>/od2od -i ${f}_fmt -o ${f}_bin Cu100_5L; done
cp <repo>/test-suite/checkpoints/Cu100_5L_photo/Cu100_5L-out.cell .
cp <repo>/test-suite/tests/testopt_photo_3s_general/Cu100_5L.odi .
<repo>/optados.x.x86_64 Cu100_5L
```

The `.odi` files have no trailing newline. Appending keywords with `>>` glues the
first one onto the `iprint` line and the run dies with "Problem reading keyword
iprint". Add a newline first.

**Order of work.** Fixes 1 and 2 are local and mechanical; do them first, one
commit each. Fix 3 is a one-line default change. Fixes 4 and 5 change the
meaning of inputs and need the user's sign-off on the design before the
benchmarks are re-blessed. Fix 6 changes non-photo tests too. 7 to 9 are minor.

---

## 1. Momentum tensor and constant-binding-energy maps count the k-point weight twice

**Severity: high. Confirmed numerically.**

### What is wrong

Every routine that walks the symmetry star of a k-point weights each image by

```fortran
k_prefactor = kpoint_weight(N_k)*total_ks/photo_n_symm()
```

where `total_ks = kpoint_grid_dim(1)*kpoint_grid_dim(2)`. The contribution it
multiplies, `qe_tsm(...)` or `qe_osm(...)` or a freshly built `temp_contribution`,
already contains `kpoint_weight(N_k)`. So each k-point enters as

    w_k * (w_k * N_grid) = w_k * m_k

with `m_k` the star multiplicity. The intent was clearly to spread the
irreducible point's weight over its images, which needs a factor `1/n_symm`
per image and nothing else. The `devel_flag : no_symmetry` branch a few lines
above each site already uses exactly that.

Consequences:

- The tensor's total is inflated by the QE-weighted mean star multiplicity, and
  no longer equals the total QE the way the E_kin vs p map does.
- The tensor's *shape* changes between a symmetry-reduced and a full k-list of
  the same physics, and with the grid, because general points (large `m_k`) are
  over-weighted relative to points on symmetry lines and at Gamma.
- A global rescale at the end (`qe_norm = total_be_kmat_contribs/total_weighted`)
  hides none of this, because `total_be_kmat_contribs` is accumulated with the
  same wrong prefactor.

### Evidence

Fixture, `testopt_photo_3s_general` input with `photo_phi_halfwidth : 180`,
`photo_output : p_tensor ekin_ptrans_map`, `photo_pmat_bin_width : 0.02`:

| run | total QE (.odo) | sum of `*_ptensor.dat` | ratio | sum of `*_Ekin_ptrans_map.dat` |
|---|---|---|---|---|
| ops present (16) | 2.757e-4 | 4.797e-3 | 17.4 | 2.757e-4 |
| ops block emptied | 4.988e-4 | 8.219e-3 | 16.5 | 4.988e-4 |
| `devel_flag : no_symmetry` | 2.757e-4 | 2.645e-4 | 0.96 | 2.757e-4 |

17.4 is exactly the QE-weighted mean of `w_k * 81` over the five points, with
the fixture's 27/7 weight renormalisation included. The 0.96 in the last row is
the `E9.1E3` output format rounding, not a residual error. (The QE itself
differs between rows 1 and 2 because `make_weights` symmetrises the polarisation
over the ops; that is expected for a reduced k-list and is item 9.)

### Sites

All in `src/photo.f90`. Each `else` branch pairs with a `no_symmetry` branch
three lines above it that is already correct.

| routine | line |
|---|---|
| `full_momentum_tensor`, 3step explicit | 5540 |
| `full_momentum_tensor`, 3step bulk | 5594 |
| `full_momentum_tensor`, 1step | 5648 |
| `accumulate_gk_tensor` | 5815 |
| `const_binding_energy_map`, 3step explicit | 6319 |
| `const_binding_energy_map`, 3step bulk | 6368 |
| `const_binding_energy_map`, 1step | 6416 |
| `const_binding_energy_map_gkgrid`, three sites | 6636, 6716, 6790 |

`grep -n 'kpoint_weight(N_k)\*total_ks' src/photo.f90` finds all ten.

### Change

Replace every `kpoint_weight(N_k)*total_ks/photo_n_symm()` with
`1.0_dp/photo_n_symm()`. In `const_binding_energy_map_gkgrid` the three sites
have no `no_symmetry` branch and set `k_prefactor` once per k-point outside the
spin loop; the replacement is the same.

Then remove what becomes dead: the `total_ks` argument of `accumulate_gk_tensor`
(5770) at its declaration and every call.

**Do not delete `total_ks` itself, nor its `kpoint_grid_dim` import**, in
`full_momentum_tensor` (5401, 5441), `const_binding_energy_map` (6235, 6251) or
`const_binding_energy_map_gkgrid` (6558, 6575). It is dead the moment the
prefactor changes, and item 10 needs it back immediately: `kpoint_weight(N_k) *
total_ks` is the only expression in the code that recovers the star multiplicity
`m_k`, and item 10's detector is built on it. Leave it in place with a comment
saying so, or item 10 begins by reverting part of this one. Build with
`-Wall -Wunused-parameter -Wunused-variable`, which the branch keeps silent.

Once every branch computes the same prefactor, the `if (index(devel_flag,
'no_symmetry'))` split in each site differs only in whether `current_k` is
rotated. Keep the flag (it is the documented way to see a single image) but
collapse the branches to one `k_prefactor` line.

### Why this is right

For a k-point with weight `w_k = m_k/N_grid` the star has `m_k` images and the
total emission from the star is `w_k * Q_k`, which is what `qe_tsm` holds. Each
image must therefore receive `w_k * Q_k / m_k`. The loop over `n_symm`
operations hits each distinct image `n_symm/m'_k` times (`m'_k` distinct
in-plane images; the z-mirror of a symmetric slab doubles up every in-plane
block), so a per-hit factor of `1/n_symm` gives each distinct image
`w_k Q_k / m'_k` and the star sums to `w_k Q_k`. With `m'_k = m_k` that is
exact; with `m'_k < m_k` (time reversal adding images the point group does not)
the star is still correctly normalised, only fewer images are placed, which is
the existing behaviour of `phi_star_fraction` and is outside this fix.

### Verify

- Re-run the three scratch cases above. All three tensor sums must now equal the
  masked QE to within the output format, and rows 1 and 3 must agree with each
  other.
- The test whose benchmark moves: `testopt_photo_1s_constbindmap_serial` and
  `_mpi` (they compare `*_ref_0.00_const_map.dat`). `testopt_photo_1s_ekinptrans_*`
  and `testopt_photo_3s_bindenergy` do not use the prefactor and must be
  bit-identical; use that as the regression check.
- No `.odo` number changes anywhere.

---

## 2. The intraband (Drude) term is not projected onto the box

**Severity: high when `optics_intraband` is on; no effect otherwise.**

### What is wrong

`calc_photo_optics` (`src/photo.f90:1340` onwards) forms, per box, the
projected interband weights

```fortran
projected_matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2) = &
  matrix_weights(n_eigen, n_eigen_final, N_k, N_spin, N2)* &
  (pdos_weights_boxes(n_eigen, N_spin, N_k, box)/pdos_weights_k_band(n_eigen, N_spin, N_k))
```

and hands them to the JDOS. The intraband term, inside the same box loop at
line 1403, does not:

```fortran
dos_matrix_weights(N_geom, n_eigen, :, :) = matrix_weights(n_eigen, n_eigen, :, :, N_geom)
```

`calc_epsilon_2` then divides both by `box_volumes(box)`. So the interband
eps_2 of a box is that box's share of the slab, per box volume, and is
layer-count independent; the intraband eps_2 of a box is the *whole slab's*
Drude weight per box volume, and doubles when the slab is doubled. Every box
also gets the identical intraband term, and it is recomputed identically
`num_boxes` times. The comment at line 1410 asserts a projection that the code
does not do; it was written about a different, also wrong, per-atom division.

### Change

Inside the `optics_intraband` block, replace the assignment with the projected
diagonal, using the same guard as the interband loop:

```fortran
do N_geom = 1, size(matrix_weights, 5)
  do N_k = 1, num_kpoints_on_node(my_node_id)
    do N_spin = 1, nspins
      do n_eigen = 1, nbands
        if (pdos_weights_k_band(n_eigen, N_spin, N_k) .le. 0.0_dp) then
          dos_matrix_weights(N_geom, n_eigen, N_k, N_spin) = 0.0_dp
        else
          dos_matrix_weights(N_geom, n_eigen, N_k, N_spin) = &
            matrix_weights(n_eigen, n_eigen, N_k, N_spin, N_geom)* &
            pdos_weights_boxes(n_eigen, N_spin, N_k, box)/pdos_weights_k_band(n_eigen, N_spin, N_k)
        end if
      end do
    end do
  end do
end do
```

Note the index order of `dos_matrix_weights` is `(geom, band, kpt, spin)` while
`matrix_weights` is `(band, band, kpt, spin, geom)`; the existing line already
handles that with array sections, the loop above spells it out. Equivalent and
shorter: `projected_matrix_weights(n_eigen, n_eigen, ...)` is already exactly
this product for the diagonal (the interband loop does not `cycle` on the
diagonal), so the block can read

```fortran
dos_matrix_weights(N_geom, n_eigen, :, :) = projected_matrix_weights(n_eigen, n_eigen, :, :, N_geom)
```

which is the one-line fix. Prefer it. Rewrite the comment at 1410 to say the
intraband weight is the box's share of the diagonal OME, projected exactly like
the interband weights, and that this is what makes per-box eps_2 independent of
how many layers the slab has.

### Verify

No photo test sets `optics_intraband`, so the suite cannot see this. Check it
by hand on the fixture: run `testopt_photo_3s_general` with
`optics_intraband : true` and `iprint : 3`, which writes per-box
`*_epsilon_*.dat`. Before the fix every box has the same intraband eps_2; after
it they differ and their sum over boxes, weighted by box volume, equals the
unprojected value times the fraction of the slab that is explicit. Then add a
test: copy `testopt_photo_3s_general`, add `optics_intraband : true` and
`optics_drude_broadening`, register it in `jobconfig`, bless it.

---

## 3. The default Fermi level moves with the k-point sampling

**Severity: medium. Changes no benchmark (all 34 pin `efermi`).**

`src/parameters.f90:316`:

```fortran
if (.not. pdis) then
  efermi_choice = "optados"
else
  efermi_choice = "file"
end if
```

`optados` integrates OptaDOS's own DOS over the spectral k-points to find E_F.
A spectral run samples the k-points that resolve the emission, not the ones
that converge the charge density, and the two can differ by a lot: on the
fixture the DOS value is -1.36 eV against the SCF's -1.06 eV (see
`Photoemission_outstanding_work.md` section 6.1(d)). That 0.3 eV goes straight
into `evacuum_eff = efermi + photo_work_function` and hence into which states
clear the vacuum level, the Fermi-Dirac factors, the binding-energy axis, and
the `cu_curve` IMFP.

### Change

Make `file` the default for photoemission as it already is for `pdis`:

```fortran
if (pdis .or. photo) then
  efermi_choice = "file"
else
  efermi_choice = "optados"
end if
```

`photo` is set at line 248 from the task string, before `efermi` is read at
line 322, so the logical is available here. Add a line to the parameter report
that prints which choice is in force for a photo run. Update the user guide
entry for `efermi` to say that for photoemission the CASTEP SCF value is used
unless overridden, and why.

### Verify

All 34 photo tests pin `efermi : -1.0615`, so they must be bit-identical.
Remove the pin from a scratch copy of `testopt_photo_3s_general` and confirm
the `.odo` header reports "Set fermi energy from file" and E_F = -1.0617 eV.

---

## 4. The surface box height is a user number, not a property of the structure

**Severity: medium. Design change; needs the user's decision before re-blessing.**

### What is wrong

In `analyse_geometry`, `z_top = photo_slab_max` (line 521) is the top of box 1.
With inferred layers every other box is bounded by centroid midpoints, so its
height is the interlayer spacing; box 1's height is
`photo_slab_max - midpoint(1,2)`, i.e. `(photo_slab_max - z_1) + gap/2`, and
`photo_slab_max` is where the user judged the electron density to have died
away. That number sets:

- `box_volumes(1)`, hence eps_2 of box 1 (goes as 1/volume), hence its
  absorption coefficient `absorp_photo(1)`, hence `I_layer(2)` and every deeper
  layer through the recursion in `calc_absorp_layer`. The reflectivity was made
  insensitive to this in `slab_reflectivity`; the attenuation was not.
- `slab_half_height = z_top - z_middle`, which scales the k_z sub-cell length
  and so the adaptive smearing width in both `calculate_delta` and the
  photoemission JDOS.
- With `photo_imfp_model : layers`, the height-weighted mean IMFP of every atom.

`photo_slab_min` is dead except for a single-layer slab: line 532 sets
`z_middle` from it and line 602 overwrites `z_middle` unconditionally on the
inferred-layer path. Its only surviving job is to be required by `param_read`.

A related robustness point: all four slab keywords are absolute Cartesian z in
Angstrom. Adding vacuum usually means re-centring the slab, and then every one
of them has to be edited by hand or the run stops (or worse, silently keeps a
different layer explicit). `calc_electron_esc` already measures depth from the
top atom, which is the frame these keywords should be in.

### Change (recommended)

On the inferred-layer path (`SLAB_MODE_BOUNDS`, `.not. single_layer`), set the
top of box 1 by the same midpoint rule as every other boundary:

```fortran
z_top = layer_centroid(1) + 0.5_dp*typical_gap
```

so all boxes have the interlayer spacing as their height and `photo_slab_max`
becomes a sanity bound only: keep the existing "atom lies above z_top" error,
and add a warning if `photo_slab_max` differs from the derived `z_top` by more
than, say, half a `typical_gap`, saying the value was not used. Print the
derived value in the geometry block. For the single-layer case and for
`SLAB_MODE_LAYERS` nothing changes: there the user has said what the box is.

Second part, optional and larger: accept the slab keywords relative to the
topmost atom (a new keyword `photo_slab_frame : absolute | top_atom`, default
`absolute` for compatibility). Not required for invariance; record as follow-up
if not done.

### Consequence

Every 3step and 1step benchmark moves, because box 1's volume changes from
`(18.282 - z_1) + gap/2` to `gap`. Measure before blessing: run
`testopt_photo_3s_general` before and after and record the QE and the per-box
volumes in the `.odo` in the commit message. If the user prefers to keep
`photo_slab_max` as the box top, the alternative is to leave the geometry alone
and only fix the dead `photo_slab_min` (drop the requirement, or use it). Ask.

### Verify

With the fix, moving `photo_slab_max` by 0.5 Angstrom in the fixture input must
leave the QE unchanged to the last digit (it currently changes it). The
`testopt_photo_3s_rm_slab_max_fail` and `_rm_slab_min_fail` tests check that the
keywords are still required; decide whether `photo_slab_min` should stay
required and update or retire that test accordingly.

---

## 5. `photo_momentum : crystal` uses the listed k-point unfolded

**Severity: medium for non-rectangular surface cells; nil for Cu(100).**

### What is wrong

`calc_angle` (`src/photo.f90:2200`, 2216) takes `kpoint_r_cart(1:2, N_k)` as the
transverse momentum. `cell_calc_kpoint_r_cart` converts the fractional
coordinates written in `.bands` with no folding. CASTEP's irreducible list
lives in the fractional box `(-1/2, 1/2]`, which coincides with the 2D
Wigner-Seitz cell only for rectangular surface lattices. For a hexagonal
surface (Cu(111), the case the branch history mentions) the corners of that
box lie outside the first BZ, and for those points `E_transverse` is
overestimated, the escape-cone test `E_kinetic < E_transverse` rejects states
that physically emit at `k_par - G`, and `theta_arpes` is wrong for the rest.
Which representative CASTEP picks depends on the grid, so the answer moves
with the k-point setup. `gkgrid` mode carries `k+G` explicitly and is unaffected.

### Change

Add a module function in `photo.f90`:

```fortran
pure function fold_to_wigner_seitz_2d(k) result(kf)
  ! Return k + G with the smallest |k + G| over the in-plane reciprocal
  ! lattice, searching n1, n2 in -2..2. recip_lattice rows 1 and 2 are b1, b2;
  ! their z components are zero because analyse_geometry has already required
  ! a and b to lie in the xy plane.
```

and apply it to every read of `kpoint_r_cart(1:2, ...)` in `photo.f90`. Sites:
`calc_angle` 2216-2217 and 2240, `kinetic_energy_momentum_map` 4778, 4823,
4859, 4893, `full_momentum_tensor` 5477, 5536, 5539 (and the two sibling
blocks), `const_binding_energy_map` 6318, 6367, 6415 and their `no_symmetry`
partners. Cleanest is to fold once into a module array `kpar_cart(2, nk)` right
after each `cell_calc_kpoint_r_cart` call and read that instead; `grep -n
'kpoint_r_cart' src/photo.f90` lists every site.

### Verify

Cu(100) is rectangular, so all 34 benchmarks must be bit-identical. For a real
check, take any hexagonal slab from the user's Cu(111) runs, list `|k_par|`
before and after folding for the irreducible set, and confirm the emission
count in the `.odo` does not fall when the grid is changed from odd to even.

---

## 6. The JDOS grid length follows the number of bands

**Severity: low to medium. Changes jdos and optics benchmarks that leave the default.**

`src/jdos_utils.f90:216`, in `setup_energy_scale`:

```fortran
jdos_max_energy = efermi - max_band_energy
```

This is negative (the sign is absorbed by an `abs()` two lines down), and its
magnitude is the largest transition *from the Fermi level*. Transitions from
deeper occupied states to the top band are larger and are silently cut off the
grid, and the cut moves whenever the DFT band count changes. In a photo run
that grid is what `calc_epsilon_1`'s Kramers-Kronig sum runs over, so
`n`, `kappa`, `R` and `alpha` at the photon energy all depend on `nbands`.

### Change

Derive the default from the full band range:

```fortran
max_band_energy = maxval(band_energy); call comms_reduce(max_band_energy, 1, 'MAX')
min_band_energy = minval(band_energy); call comms_reduce(min_band_energy, 1, 'MIN')
jdos_max_energy = max_band_energy - min_band_energy
```

with both broadcast, positive, and reported. For photo runs, add a one-line
warning when `jdos_max_energy` was not given, saying the optical constants
depend on it and it should be set explicitly and converged.

### Consequence

None of the three `testopt_jdos_*` or seven `testopt_optics_*` inputs sets
`jdos_max_energy` (checked: `grep jdos_max_energy` over their `.odi` files
matches nothing), so all ten get a longer grid and a different eps_1 and must
be re-blessed. Record the old and new grid length and the largest relative
change in eps_1 at a representative energy in the commit message, so the
change is traceable the way the item 4 regrid was in
`Photoemission_outstanding_work.md` section 2.

---

## 7. `kpoint_grid_dim` is guessed from the k-point list

**Severity: low. Note and warn; no numeric change.**

`cell_find_MP_grid` infers the MP grid from the list in `.bands`. It feeds
`step(:)` in every adaptive and linear broadening width (`calculate_delta`
3197, `calculate_jdos`, `calculate_dos_at_e`) and the k broadening of the maps
(4762, 5443). For a `spectral_kpoints_list` it can misfire; `kpoint_mp_grid` in
the `.odi` overrides it. After fix 1 it no longer enters any prefactor.

Change: when `photo` is set and `kpoint_mp_grid` was not given, print the
inferred grid next to the "Grid size" line with a note that it was inferred.
Nothing else.

---

## 8. `total_field_emission` grows with the number of layers

**Severity: low. Printed only.**

`src/photo.f90:1804` sums the tunnelling probability times occupation over
*every* band and k-point and divides by the cell area. More layers means more
occupied states per area, so the printed number scales with slab thickness and
says nothing intensive. It is not used anywhere else (`grep total_field_emission`).

Change: either drop the print, or restrict the sum to states within a few kT
of E_F and label it as an estimate per unit area of the states that can
tunnel. Drop is simpler; `testopt_photo_3s_field` compares the `.odo` through
`test-suite/tools/parsers/parse_photo_odo.py`, which extracts `fermi_fb`,
`fermi_ab`, `fermi_lb`, `qe_bulk`, `qe_total`, `mte`, `layer_imfp` and the three
`ds_*` values and nothing else, so the benchmark is unaffected either way.

---

## 9. Smaller things noticed on the way

- **DS model bin truncation.** `src/photo.f90:2780`,
  `delta_index_photon = int(temp_photon_energy/delta_e)`, is the same
  truncation the previous round fixed for the photon sweep (items 4 and 5 in
  `Photoemission_outstanding_work.md`). Use `nint`. Only `testopt_photo_ds_general`
  can move.
- **The `_nosymm` fixtures are inconsistent input.** They pair the
  symmetry-reduced five-point k-list with an emptied `symmetry_ops` block. That
  is a legitimate code-path test for the zero-operation guard, but the numbers
  they produce are not comparable to the symmetric ones (QE 4.99e-4 against
  2.76e-4, because `make_weights` no longer symmetrises the polarisation). Say
  so in the checkpoint `README.md` so nobody reads the difference as a bug.
- **Comment at `src/photo.f90:1410`** is wrong as described in fix 2; rewrite it
  when doing that fix.

---

## 10. The symmetry star is taken to be what the operation list generates, and time reversal is not in it

**Severity: high for any azimuth-resolved output on a slab whose two faces
differ. No effect on `qe_total`, and no effect at all on a centrosymmetric slab.**

### What is wrong

CASTEP reduces the k-mesh with the point group **and time reversal**. The
`symmetry_ops` block it writes contains only the point group. So for a structure
whose in-plane point group does not contain `-I`, the k-point weights encode a
star that the operation list cannot reproduce:

    m_k   = kpoint_weight(N_k) * total_ks     the star CASTEP used
    m'_k  = |orbit of k under the n_symm ops| the star OptaDOS can place

with `m_k = 2 m'_k` for every k that is not its own `-k` modulo **G**.

Every routine that places a star -- `full_momentum_tensor`,
`const_binding_energy_map`, `const_binding_energy_map_gkgrid`,
`accumulate_gk_tensor` -- and also `phi_star_fraction`, treats the operation list
as the whole star. After fix 1 the **total** is still exactly `w_k Q_k`
regardless (the per-hit `1/n_symm` cancels the `n_symm/m'_k` orbit
over-counting), so no `.odo` number moves. What is wrong is the **placement**:
half the images are missing from the map, and the half that are placed are
twice as bright as they should be.

`photo.f90` says as much already, in `phi_star_fraction`'s own description of a
star as what the operations generate; this item is about the case where that is
not the whole star.

### Evidence

Built for this purpose, 8 Sep 2026, and run through CASTEP 23.1 `--dryrun`:

| slab | ops in block | point group | irreducible k from 9x9x1 |
|---|---|---|---|
| 16-layer Cu(111) 1x1 | 12 | D3d, `-3m` | 12 |
| same + 1 H at the fcc hollow | **6** | C3v, `3m` (P3m1) | **12** |

Same grid, same number of irreducible points, an **identical multiset of
weights** -- but half the operations. CASTEP reduced the mesh as if it had 12,
and the missing six are time reversal, recorded nowhere in the files OptaDOS
reads.

Orbits computed under the operations actually written (rotations converted to
fractional; note the block is Cartesian -- see `photo_symm_2d`):

| slab | k-points with `m_k = m'_k` | with `m_k = 2 m'_k` |
|---|---|---|
| 16-layer Cu(111) | 12 of 12 | 0 |
| + H adsorbate | 5 of 12 | **7 of 12** |

The seven carry **74 % of the BZ weight** (4x0.0741 + 3x0.1481). The five
unaffected are Gamma, K, and three points on mirror lines, each its own `-k`
modulo **G**. The ratio was 1 or 2 and never anything else, which is not luck:
time reversal is an order-2 extension of the point group, so orbit-stabiliser
allows only those two values.

**It does not arise anywhere in the current campaign.** All 411 `-out.cell`
files under `~/structures` were checked: 387 contain `-I` in-plane, and the
other 24 have a single operation (no `symmetry_ops` block), where each k stands
for itself and `m_k = m'_k = 1`. The condition to watch for is a slab whose two
faces differ -- an adsorbate on one side, a polar termination, a supported film.
An alkali-antimonide photocathode slab would qualify.

### Change

Detection is nearly free, because fix 1 removes an expression that already
contains both numbers. Once per k-point:

```fortran
m_k     = kpoint_weight(N_k)*total_ks          ! already computed today
m_prime = photo_star_size(N_k)                 ! |orbit| under the n_symm ops
```

`photo_star_size` folds each image with `k - floor(k + 0.5)` and counts
distinct results -- the idiom `cell_find_MP_grid` already uses.

**Implemented 8 Sep 2026, and more simply than first written here.** The per-k
`n_tr` above is unnecessary. Doubling is harmless for a k-point that does not
need it: when `-k` is already in the orbit the negated images coincide with those
already placed, so each distinct image is hit twice and the halved prefactor
returns it exactly the weight it had. One global flag therefore suffices, and the
ten placing loops need no per-k branch and no restructuring at all -- only a
change of which list they walk:

```fortran
pure function photo_star_n_ops() result(n_ops)   ! n_symm, doubled when needed
pure function photo_star_op(i)   result(mat)     ! photo_symm_2d(i), negated beyond n_symm
```

`photo_n_symm()` -> `photo_star_n_ops()` and `photo_symm_2d(nsymm_op)` ->
`photo_star_op(nsymm_op)` throughout the placing routines and in
`phi_star_fraction`; 37 call sites, all mechanical. Negating the operation is the
same as negating k, since `{-R k} = {R (-k)}` and the group is closed. When the
star is already complete `photo_star_n_ops()` returns `photo_n_symm()` unchanged,
so nothing whatsoever moves on a centrosymmetric structure.

Two `devel_flag` values, matching the `no_symmetry` convention:

- `no_time_reversal` -- detect and report, place nothing extra. The A/B control.
- `force_time_reversal` -- double even when the star is complete. This exists to
  test the claim the design rests on: forcing it on a centrosymmetric fixture
  must reproduce the run bit for bit. Verified on `Cu100_5L`: **zero differing
  lines**.

`phi_star_fraction` needs the same gate, or the code carries two different
definitions of "star".

**Self-check, and why it is worth more than the correction.** `m_k` must be a
positive integer and `m_k/m'_k` must be exactly 1 or 2. Both are theorems, not
tolerances. Anything else means `kpoint_grid_dim` is wrong, and the right
response is to warn and decline to correct rather than guess. This is the first
automatic test `cell_find_MP_grid` has ever had.

### On trusting `kpoint_grid_dim`

`cell_find_MP_grid` (`cell.f90:94`) infers the grid from the **spacings** of the
time-reversed, folded, per-dimension value set. Two defects bear on this item.

**1. The documented aliasing bug** (its own comment, AJM 15/9/2014). Time
reversal maps the set `S` to `S` union `-S`. An unshifted grid is already
TR-symmetric so nothing changes; a shifted grid generically doubles, and when
the shift is an odd multiple of `1/(4n)` the union is *exactly uniform* with
spacing `1/(2n)` -- algebraically the same set of numbers as an unshifted `2n`
grid. No spacing-based test can separate them, because the distinguishing
information is destroyed by the union, not overlooked. Reproduced by
transcribing the routine: true 2x1x1 shifted by 1/8 or 3/8 -> reports 4; true
4x1x1 shifted by 1/16 -> reports 8.

**2. `int` where `nint` is meant, undocumented.**

```fortran
kpoint_grid_dim(idim) = int(1.0_dp/min_img)        ! cell.f90, both branches
```

For a true 3x1x1 shifted by 1/12 the uniform TR set has spacing `1/6`, and
`1.0/min_img` evaluates to `5.999999999988`, so `int` returns **5** -- neither
the true 3 nor the aliased 6. Rounding in the file happens to protect the
unshifted cases (a 6 dp k-point makes `min_img` slightly small, so `1/min_img`
lands just *above* the integer and truncation is harmless), which is why this
has never shown up. Use `nint`. Same class as item 9.

Neither defect blocks this item: a wrong grid makes `m_k` wrong by the same
factor, and the integer-and-ratio self-check above fires rather than silently
mis-correcting. `kpoint_mp_grid` in the `.odi` overrides the inference and is
the reliable route; item 7 already proposes reporting which of the two was used.

Worth stating plainly: on a symmetry-reduced k-list the inference is not
reliable in principle. The per-dimension projection of an irreducible set need
not be TR-symmetric even for an unshifted grid, so the pre- versus post-TR count
does not rescue it either. Prefer the user's value; validate whatever is used.

### Verify

- `Cu111_H` (16 Cu + 1 H at the fcc hollow, `--dryrun`, 9x9x1): 6 operations,
  C3v, 12 k-points, and 7 of them with `m_k = 2 m'_k`. The control, the same
  slab without the H, must give 12 operations and `m_k = m'_k` throughout.
- On any current campaign structure the detector must never fire.
- With the correction on, the const-binding-energy map of `Cu111_H` must gain
  six-fold structure where it had three-fold, at unchanged total.
- `testopt_photo_*` benchmarks must be **bit-identical**: no `.odo` number
  depends on this.

**The existing fixture cannot test this.** `Cu100_5L` carries weights
4/21, 4/21, 8/21, 4/21, 1/21 against `kpoint_mp_grid : 9 9 1`, so
`m_k = w_k*81` is 15.43, 30.86, 3.86 -- not integers. Its five-point list is a
truncated subset with renormalised weights and `m_k` is not recoverable from it;
the self-check would correctly refuse. A case-B test needs its own fixture with
genuine CASTEP weights.

---

## What was checked and found sound

So the implementer does not re-audit it:

- `cell_area = cell_volume/real_lattice(3,3)` and `box_volumes = box_heights*cell_area`
  are vacuum independent; `analyse_geometry` refuses tilted `a`, `b`.
- Per-atom emission is `states per cell area` times `pdos fraction on the atom`,
  and both scale oppositely with layer count, so per-layer QE converges with
  slab thickness by construction. The bulk term repeats the deepest explicit
  box with `bulk_repeat`, not the box height.
- The k_z sub-cell length is replaced by `pi/slab_half_height` and
  `pi/bulk_repeat`, so no `recip_lattice(3,3)` reaches any width.
- Escape depths are measured from the topmost atom (`calc_electron_esc`).
- OMEs from `ome_bin` are per-cell normalised and vacuum independent;
  `make_weights` has no volume factor.
- ~~The FEM tensor lookup is referenced to `evacuum_eff`, which the docblock at
  `make_foptical_weights` 3376 records as measured to remove the vacuum
  dependence of the cell-averaged potential.~~

  **RETRACTED 2026-09-08. This was wrong, and it was the only item here that
  was cleared by citing a code comment rather than by checking against the
  file it describes.** The `.fem_bin` energy axis is *absolute*, on the same
  eigenvalue scale as `band_energy`: CASTEP writes
  `Ef_min = fem_ef_origin + SPECTRAL_FEM_EF_MIN` (`spectral.f90:552`) with
  `fem_ef_origin` the SCF Fermi energy. Subtracting `evacuum_eff` therefore
  moved an already-consistent lookup off by roughly the work function, and the
  invariance the docblock claimed to create is automatic without it.

  Measured: with the term removed, 1-step QE is exactly invariant under a +3 eV
  shift of the eigenvalue zero (ratio 1.0000 at 51 photon energies); with it,
  QE moves by 3-14x. All 1-step campaign results are invalid; 3-step is
  untouched, as it never opens `fem_bin`. Fixed at `photo.f90:3397`. Full
  evidence in `evidence_dependence/One_step_FEM_energy_reference_bug.md`.
- `kpoint_weight` from `.bands` sums to one; `electrons_per_state` handles spin.

## Suggested commit sequence

1. Fix 1 (prefactor) + re-bless the two `constbindmap` benchmarks.
2. Fix 2 (intraband projection) + new intraband test.
3. Fix 3 (efermi default) + user guide line.
4. Fix 9, DS `nint`, together with the `int` -> `nint` in `cell_find_MP_grid`
   (item 10) -- same class, one commit. Re-bless `ds_general` if it moves.
5. Fix 6 (JDOS grid) + re-bless affected jdos/optics tests.
6. Fixes 7 and 8.
6b. Fix 10 (time-reversed star). After fix 1, since it builds on the prefactor
   it replaces and on the `total_ks` that fix 1 must therefore keep. Land the
   detector and its self-check first, on its own, and confirm it stays silent on
   every campaign structure; the correction and the `Cu111_H` fixture second.
   No benchmark may move in either commit.
7. Fix 5 (Wigner-Seitz folding), after confirming Cu(100) is bit-identical.
8. Fix 4 (surface box), last, after the user has chosen the design, with the
   full photo re-bless and the measured before/after QE in the message.

Each commit message should state which benchmarks moved and by how much, in the
style of the existing history on this branch.
