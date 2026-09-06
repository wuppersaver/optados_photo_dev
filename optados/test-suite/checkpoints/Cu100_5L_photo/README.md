# Cu100_5L photoemission checkpoint

`Cu100_5L.cell` + `Cu100_5L.param` regenerate every fixture in this directory.
Run CASTEP on them, then convert each `*_bin` to `*_fmt` with `od2od` and `bzip2`
the result. The `.bands` and `-out.cell` come straight from the same run.

Generated 6 Sep 2026 with CASTEP 27.1, `aec4cfa88 reuse-pad-zero`.

## Things that are deliberate, and must not be "corrected"

**`SPECTRAL_TASK: PHOTOOPTICS`, not `All`.** `All` sets `calc_core`/`calc_dome`/
`calc_ome` only. `fem_bin`, `tmprob_bin` and `gkgrid_bin` — which every
photoemission test needs — are written only by `PHOTOOPTICS`. The `.param` that
used to live here said `All` and therefore could not produce the fixtures beside
it.

**The five spectral k-points are the ones that emit.** Emission needs
`hbar^2 k_par^2 / 2m < hv - work_function`, which at hv = 6 eV and phi = 4.2202 eV
is `|k_par| < 0.683 1/Ang`. Of the fifteen k-points in the 9x9x1 irreducible
wedge, ten lie outside that cone and contribute *exactly* zero; the five inside
carry the whole quantum efficiency. Their measured shares on the full wedge:

| k-point | \|k_par\| | share of QE |
|---|---|---|
| (-1/9,  0  ) | 0.270 1/Ang |  2.3 % |
| (-1/9, 1/9 ) | 0.382 1/Ang | 77.0 % |
| (-1/9, 2/9 ) | 0.605 1/Ang | 13.1 % |
| ( 0  ,-2/9 ) | 0.541 1/Ang |  7.5 % |
| ( 0  ,  0  ) | 0.000       |  0.2 % |

The previous fixture used a 3x3x1 list whose only in-cone point was Gamma. Since
Gamma has zero transverse momentum and maps to itself under all sixteen symmetry
operations, that fixture exercised neither the transverse-momentum machinery nor
the symmetry star, and one narrow gaussian carried the entire answer — which is
how a 15 % error in the JDOS energy grid survived in this suite undetected.

**The k-point weights are renormalised, and the QE is not physical.** In the
9x9x1 wedge these five carry 4/81, 4/81, 8/81, 4/81 and 1/81, summing to 7/27.
CASTEP requires spectral weights to sum to one (`cell_read_spectral_data` aborts
otherwise), so they are stored as 4/21, 4/21, 8/21, 4/21, 1/21. That preserves
the five points' relative weight exactly and rescales the sum by 27/7 = 3.857.

The resulting QE is therefore **deterministic but not physical, on purpose**. The
optics is sampled on these five points too, so it does not reproduce the
full-wedge value either. A regression fixture does not need the converged answer;
it needs the k-point loop, the weights, the transverse momentum, the escape cone,
the angular acceptance and the symmetry star all to do real work. Do not compare
these numbers with a converged calculation, and do not renormalise the weights
again.

**`SPECTRAL_FEM_EF_MAX : 15.0 eV`.** The one-step model looks the coherency
tensor up at the *free electron* final-state energy. For the topmost occupied band
at hv = 6 eV that is 2.04 eV, and the window is shifted down by the Fermi energy
because `fem_energy_info(2) = SPECTRAL_FEM_EF_MIN + E_F`. A nominal `EF_MAX` of 3
would only reach 1.94 eV and `make_foptical_weights` aborts. 15 gives
`ef_top = 13.94 eV`, headroom to hv ~ 18 eV, and 80 energy bins instead of the
default 120 — which is most of what pays for five k-points rather than three.

**65 bands is close to the floor, not padding.** Excluding the top bands and
re-running the three-step model on the full wedge:

| bands kept | QE | change |
|---|---|---|
| 65 | 6.932e-05 | — |
| 59 | 6.932e-05 |  0.00 % |
| 49 | 6.932e-05 |  0.00 % |
| 39 | 6.925e-05 | -0.10 % |
| 35 | 6.299e-05 | -9.1 % |
| 33 | 5.824e-05 | -16 % |

The three-step final states at hv = 6 eV sit at bands 32-33, but most of the loss
below 40 bands is *optics*: `epsilon_1(6 eV)` comes from a Kramers-Kronig integral
over the whole spectrum, so truncating the conduction manifold biases the
absorption, hence `I_layer`, hence the QE. Cutting bands to shrink the fixture is
not available.

## The `efermi` pin in the test inputs

Every photo `.odi` pins `efermi` rather than letting OptaDOS compute it. That is
required, not habit: the SCF runs on a dense 9x9x1 grid, but the spectral outputs
are written for five k-points, and OptaDOS's own Fermi energy from that sample is
wrong. Measured on the previous fixture set:

    CASTEP dense SCF (what the .bands carries)   -1.0617 eV
    OptaDOS from the 3 sampled k-points          -1.3614 eV
    difference                                    0.30 eV

300 meV goes straight into the emission threshold and the effective work function.
**Whenever these fixtures are regenerated, the pin must be updated to the Fermi
energy in the new `.bands`.** For this set that is `-0.039010 Ha = -1.0615 eV`.

## File versions

`tmprob_bin`, `gkgrid_bin` and `fem_bin` are **version 2**; `ome_bin`, `dome_bin`
and `pdos_bin` are version 1. Version 2 is not a relabelling — the transmission
probability is normalised, the G+k weights are normalised by the plane-wave norm,
and the `fem` file holds the coherency tensor rather than wavepacket matrix
elements. OptaDOS refuses version 1 for the first three. `od2od` preserves the
stamps; it did not before 6 Sep 2026, when it wrote 1.0 into everything and
silently downgraded any version-2 file passed through it.
