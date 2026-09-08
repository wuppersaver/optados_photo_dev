# testopt_photo_1s_shift_invariance

1-step FEM lookup: invariance under a shift of the eigenvalue zero.

**Not yet runnable.** Written 2026-09-08 alongside the fix at `photo.f90:3397`.
This directory holds the design only — the shifted fixture has still to be
generated, and there is no `[testopt_photo_1s_shift_invariance]` section in
`tests/jobconfig`, so the harness does not pick it up. It takes about 0.4 s once
the fixture exists.

Named to match its siblings (`testopt_photo_1s_general`, `_layers`,
`_constbindmap_*`) and placed in `tests/` for the same reason, so that adding the
fixture and one `jobconfig` stanza is all that is left to do.

Unlike its siblings it is an **invariance** test, not a benchmark comparison: it
runs the code twice and compares the two runs against each other, so it needs no
blessed `benchmark.out.*` file of its own.

## What it asserts

The `.fem_bin` energy axis is absolute, on the same eigenvalue scale as
`band_energy`. Shifting the zero of that scale by a constant therefore cannot
change any physical result: every energy the 1-step model uses moves together.

So: take the 5L photoemission checkpoint, add one constant `DELTA` to
**everything on the eigenvalue scale**, and require **bit-identical** 1-step QE.

## Why this test and not a numerical benchmark

This is a class of bug no blessed-number benchmark can catch, because every
blessed number was blessed while the code was broken. The `testopt_photo_1s_*`
benchmarks all passed throughout, on results that were wrong by a factor of
3–14. An invariance test does not need a correct reference value — it compares
the code against itself under a transformation that must be a no-op.

It is also exact rather than tolerance-based: bit-identical, not "within 1e-6".

Had it existed, it would have failed on day one.

## What must be shifted, all by the same DELTA

Miss any one of these and the test fails for the wrong reason:

| file | field |
|---|---|
| `<seed>.bands` | every eigenvalue, **and** the `Fermi energy` header |
| `<seed>.fem_fmt` | `Ef_min` (slot 2) **and** `Ef_origin` (slot 5) |
| `<seed>.odi` | the pinned `efermi`, if the run sets one |

`Ef_step` and `Ef_broadening` are differences, not positions — they must **not**
move. Neither must the work function, which is `E_vac − E_F` and so is itself
invariant.

Note `Ef_origin` and `Ef_min` both shift: `Ef_min` already has `Ef_origin`
folded into it (`spectral.f90:552`), so they move together rather than one
compensating the other. Applying `Ef_origin` a second time as a correction is
the mistake §6 of the bug write-up warns against.

## Procedure

1. Copy the 5L fixture to `shifted/`.
2. Add `DELTA` (3.0 eV was used in the investigation) to the fields above.
3. Run OptaDOS 1-step on both, same `.odi` otherwise.
4. Compare the QE column of the two `.odo` files.

**Pass:** every value bit-identical.
**Fail:** any difference at all — the lookup is reading a shifted axis as if it
were absolute.

## Expected result either side of the fix

| binary | outcome |
|---|---|
| before `photo.f90:3397` | QE ratio 3–14x across the sweep — **fails** |
| after | ratio 1.0000 at all 51 photon energies — **passes** |

Measured in §5 of `evidence_dependence/One_step_FEM_energy_reference_bug.md`
with two locally built binaries differing by that one line.

## Companion assertion

`elec_read_foptical_mat` now refuses a `.fem_bin` whose `Ef_origin` disagrees
with the `.bands` Fermi energy by more than 1 meV. **The shifted fixture must
shift both**, or it will trip that check instead of testing the lookup — which
would still be a correct failure, just not the one this test is for.

## What it does not cover

Nothing about whether the *absolute* QE is right — only that it is independent
of an arbitrary energy origin. A systematic error affecting both runs equally
passes. Pair it with the `MTE <= hv - phi` bound check (§10.2 of the write-up,
not implemented) for the complementary direction.
