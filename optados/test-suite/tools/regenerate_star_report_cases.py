#!/usr/bin/env python3
"""Regenerate documents/Photoemission_star_report_cases.md.

Every block in that document is captured from a real OptaDOS run, not typed by
hand, so the document cannot drift away from the code. Re-run this after any
change to the report in ``photo_check_star`` and commit the result alongside.

    ./test-suite/tools/regenerate_star_report_cases.py

Needs the serial binary ``optados.x.x86_64`` built and, on this machine, the
Intel environment sourced first (``source /opt/intel/oneapi/2025.1/oneapi-vars.sh``)
or the link fails on MKL.

Each case is assembled from fixtures already in the repository -- the two photo
checkpoints and nothing else -- so a fresh clone reproduces it. The runs all stop
for want of an ``ome_bin``; that is deliberate, and harmless here, because the
report is written before the first matrix element is read.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HA_TO_EV = 27.2113845

HERE = Path(__file__).resolve().parent
TESTSUITE = HERE.parent
ROOT = TESTSUITE.parent
BINARY = ROOT / 'optados.x.x86_64'
OUTPUT = ROOT / 'documents' / 'Photoemission_star_report_cases.md'

TRSTAR = TESTSUITE / 'checkpoints' / 'Cu111_TRstar_photo'
CU5L = TESTSUITE / 'checkpoints' / 'Cu100_5L_photo'

ODI = """task                 : photoemission
photo_model          : dosds
photo_photon_energy  : 6.0
jdos_spacing         : 0.05
broadening           : adaptive
efermi               : {efermi}
optics_geom          : unpolar
optics_qdir          : 1 1 0
photo_slab_min       : {slab_min}
photo_slab_max       : {slab_max}
photo_work_function  : 4.60
photo_imfp_model     : const
photo_imfp_value     : 5.3
photo_elec_field     : 0.0
photo_temperature    : 300.0
photo_momentum       : crystal
iprint               : 1
{extra}"""

# label, seed, checkpoint, slab window, extra .odi lines, what it demonstrates
CASES = [
    ('A', 'Cu111_star_sym', TRSTAR, (4.0, 18.0), '',
     'A complete symmetry star: nothing to do',
     'The 6-layer Cu(111) slab, D3d, 12 operations, on an odd 9x9x1 grid. The '
     'operations reach every k-point the weights describe, so the section is '
     'three lines and no prose. **This is what every structure in the current '
     'campaign produces.**'),

    ('B', 'Cu111_star', TRSTAR, (4.0, 18.0), '',
     'Time reversal needed, and applied',
     'The same slab with a seventh Cu on the hcp hollow: C3v, 6 operations. '
     'CASTEP reduced the mesh using time reversal as well, so 7 of the 12 '
     'k-points stand for twice as many emission directions as the operation '
     'list can generate. The missing half is placed.'),

    ('C', 'Cu111_star', TRSTAR, (4.0, 18.0), 'devel_flag           : no_time_reversal\n',
     'Time reversal needed, deliberately not applied',
     'The same run with `devel_flag : no_time_reversal`. This is the A/B '
     'control for the correction, and the state the code was in before item 10.'),

    ('D', 'Cu100_5L', CU5L, (7.559, 18.282), '',
     'k-point weights are not those of a whole mesh',
     'The Cu(100) 5-layer photoemission fixture. Its five k-points carry '
     'weights renormalised over 21 against a 9x9x1 grid, so how many emission '
     'directions each stands for cannot be recovered. Nothing is corrected.'),

    ('E', 'Cu111_star_even', TRSTAR, (4.0, 18.0), '',
     'The mesh is not invariant under the point group',
     'The symmetric slab again, but on an **even** 4x4x1 grid. On a hexagonal '
     'cell the half-step offset an even Monkhorst-Pack grid carries suits a '
     'four-fold axis, not a three-fold, so rotating a k-point lands off the '
     'mesh. Measured: hexagonal grids 3, 5, 9, 21, 33, 51 are invariant and '
     '4, 6, 8, 12, 16, 32 are not; square lattices are invariant either way. '
     'This breaks every output built by placing a star, with or without '
     'item 10 -- the check is only what makes it visible.'),

    ('F', 'Cu111_star_sym', TRSTAR, (4.0, 18.0), 'kpoint_mp_grid       : 18 18 1\n',
     'The stated grid does not match the k-point list',
     'Case A with `kpoint_mp_grid : 18 18 1` in the `.odi` when the mesh is '
     'really 9x9x1 -- an ordinary typo. A run whose grid was *reconstructed* '
     'rather than stated prints only the first two lines of this block; there '
     'is no separate wording for it, because every route to a mis-reconstructed '
     'grid also breaks mesh invariance and is caught by case E first.'),
]


def fermi_ev(bands: Path) -> float:
    head = bands.read_text(errors='replace')[:400]
    match = re.search(r'Fermi energ\w*\s*\(in atomic units\)\s*([-\d.E+]+)', head)
    if match is None:
        raise SystemExit(f'no Fermi energy in {bands}')
    return float(match.group(1)) * HA_TO_EV


def run_case(seed, checkpoint, slab, extra, workdir: Path) -> str:
    shutil.copy(checkpoint / f'{seed}-out.cell', workdir)
    bands = workdir / f'{seed}.bands'
    bands.write_bytes(subprocess.run(
        ['bunzip2', '-kc', str(checkpoint / f'{seed}.bands.bz2')],
        capture_output=True, check=True).stdout)
    (workdir / f'{seed}.odi').write_text(ODI.format(
        efermi=f'{fermi_ev(bands):.4f}', slab_min=slab[0], slab_max=slab[1], extra=extra))

    subprocess.run([str(BINARY), seed], cwd=workdir,
                   capture_output=True, timeout=600)

    odo = workdir / f'{seed}.odo'
    if not odo.is_file():
        raise SystemExit(f'{seed}: no .odo produced')
    lines = [l.rstrip() for l in odo.read_text(errors='replace').splitlines()]
    try:
        first = next(i for i, l in enumerate(lines) if 'Symmetry Star Bookkeeping' in l)
    except StopIteration:
        raise SystemExit(f'{seed}: no star section in the .odo')

    # The box, up to and including its closing rule.
    last = next(i for i in range(first + 1, len(lines))
                if lines[i].strip().startswith('+---'))

    # Then any warning block that follows it. Warnings are a separate house-style
    # box, so they are part of what this section prints and must be captured too;
    # capturing only the informational box would document half the behaviour.
    i = last + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i < len(lines) and lines[i].strip().startswith('|') and not lines[i].strip('| '):
        i += 1                                     # the spacer line
    if i < len(lines) and lines[i].strip().startswith('!---'):
        end = next(j for j in range(i + 1, len(lines))
                   if lines[j].strip().startswith('!---'))
        return '\n'.join(lines[first:last + 1] + [''] + lines[i:end + 1])
    return '\n'.join(lines[first:last + 1])


def main() -> int:
    if not BINARY.is_file():
        raise SystemExit(f'{BINARY} not built. make optados BUILD=fast COMMS_ARCH=serial')

    parts = [
        '# Photoemission: every message the symmetry-star section can print',
        '',
        'Generated by `test-suite/tools/regenerate_star_report_cases.py`. Every',
        'block below is captured from a real run of `optados.x.x86_64`, so this',
        'document cannot drift away from the code. Re-run the script after any',
        'change to the report in `photo_check_star` and commit the result.',
        '',
        'Background: item 10 of `Photoemission_setup_invariance_fixes.md`.',
        '',
        'The section is written before the first matrix element is read, which is',
        'why these runs need only a `.bands` and an `-out.cell`, and why they all',
        'stop afterwards for want of an `ome_bin`.',
        '',
        '## The cases at a glance',
        '',
        '| | outcome | corrected? |',
        '|---|---|---|',
        '| A | complete star, nothing needed | n/a |',
        '| B | time reversal needed | yes |',
        '| C | time reversal needed | no, `devel_flag : no_time_reversal` |',
        '| D | weights not those of a whole mesh | no, cannot be established |',
        '| E | mesh not invariant under the point group | no, cannot be established |',
        '| F | stated grid does not match the k-point list | no, cannot be established |',
        '',
        'Only **A** should appear on a well-posed photoemission run. **B** is',
        'correct behaviour on a slab whose two faces differ. **D**, **E** and **F**',
        'all mean the bookkeeping could not be established, and nothing is changed;',
        'of those, **E** reports a defect that affects the run whether or not this',
        'check exists.',
        '',
    ]

    with tempfile.TemporaryDirectory() as tmp:
        for tag, seed, cp, slab, extra, title, blurb in CASES:
            work = Path(tmp) / tag
            work.mkdir()
            print(f'  case {tag}: {seed}', file=sys.stderr)
            block = run_case(seed, cp, slab, extra, work)
            parts += [f'## {tag} — {title}', '', blurb, '', '```', block, '```', '']

    OUTPUT.parent.mkdir(exist_ok=True)
    OUTPUT.write_text('\n'.join(parts))
    print(f'wrote {OUTPUT} ({len(parts)} blocks)', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
