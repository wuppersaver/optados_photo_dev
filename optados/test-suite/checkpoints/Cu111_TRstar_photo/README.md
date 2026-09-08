# Cu(111) time-reversed-star checkpoint

A matched pair of Cu(111) 1x1 slabs that differ in one thing only: whether the
two faces are the same.

| seed | atoms | stacking | grid | point group | in-plane `-I`? | ops |
|---|---|---|---|---|---|---|
| `Cu111_star_sym` | 6 Cu | A C B A C B | 9x9x1 | D3d, `-3m` | yes | 12 |
| `Cu111_star` | 7 Cu | A C B A C B **C** | 9x9x1 | C3v, `3m` (P3m1) | **no** | **6** |
| `Cu111_star_even` | 6 Cu | A C B A C B | **4x4x1** | D3d, `-3m` | yes | 12 |

`Cu111_star_even` is the symmetric slab on an **even** grid. On a hexagonal cell
that mesh is not invariant under the point group, which is a defect in star
placement generally and not only in the time-reversal bookkeeping. It exists so
that message has a fixture too; see case E of
`documents/Photoemission_star_report_cases.md`.

`<seed>.cell` + `<seed>.param` regenerate both. The `.bands` and `-out.cell` come
straight from that run; nothing else is needed here, because the symmetry-star
report is written before OptaDOS opens its first matrix element.

Generated 8 Sep 2026 with CASTEP 23.1 (`425c8e1c2 photooptics`).

## What it is for

Item 10 of `documents/Photoemission_setup_invariance_fixes.md`. CASTEP reduces
the k-mesh using the point group **and time reversal**, but writes only the point
group into `symmetry_ops`. When the in-plane group has no `-I`, the k-point
weights describe a star twice as large as the operation list can reach, and
nothing in the files says so.

Both slabs give **12 irreducible k-points from the same 9x9x1 grid with an
identical multiset of weights**, despite one having half the operations. That is
the whole point of the fixture: the difference is invisible in the k-point list
and shows up only when the weights are compared against the orbits.

| seed | `m_k = m'_k` | `m_k = 2 m'_k` |
|---|---|---|
| `Cu111_star_sym` | 12 of 12 | 0 |
| `Cu111_star` | 5 of 12 | **7 of 12**, carrying 74 % of the weight |

The five unaffected points are Gamma, K, and three on mirror lines -- each its
own `-k` modulo **G**.

## Things that are deliberate, and must not be "corrected"

**The seventh atom sits on the hcp hollow, not the fcc one.** That is what breaks
the inversion centre while keeping the three-fold axis. Two near misses, both
made while building this:

- An adatom on the **fcc continuation** site adds nothing: `A C B A C B A` is
  simply a seven-layer fcc slab, still centrosymmetric about its middle layer.
  CASTEP returns D3d and 12 operations, unchanged.
- An **H** adatom breaks the symmetry correctly, but `analyse_geometry` then
  refuses the structure -- *"the inferred layers are not all the same material"* --
  before the star is ever reported. Keeping every atom Cu avoids needing
  `photo_layers_tops` just to reach the thing under test.

**It is not converged, and does not need to be.** The SCF settings are chosen for
speed. This fixture tests symmetry and k-point bookkeeping, which depend on the
k-point list, the weights and the operations -- not on the eigenvalues. Do not
raise the cutoff to make the numbers "better"; no number here means anything
physically.

**Six layers, not sixteen.** The in-plane group and the k-point set are identical
either way, so a thicker slab bought nothing and cost minutes per run.

**The runs stop for want of an `ome_bin`, and that is expected.** See the
`can_fail` note on `OPTADOS_PHOTO_STAR_OK` in `tests/userconfig`.
