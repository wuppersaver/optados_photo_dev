#!/usr/bin/env python3
"""Layer boundaries for an OptaDOS photoemission run, from the charge density.

Reads a CASTEP ``.den_fmt`` and the matching ``-out.cell`` and prints the
``photo_layers_tops`` and ``photo_slab_middle`` block to paste into an ``.odi``.

Why the density rather than the atom positions
----------------------------------------------
Left to itself OptaDOS puts the boundary between two layers at the midpoint of
their centroids, and takes the top of the slab from ``photo_slab_max``, which the
user picks by eye from where the density has fallen off.  Those are two unrelated
rules for the same kind of quantity, and the second one leaks: the box volumes
set the dielectric functions, so a half-angstrom change in ``photo_slab_max``
moved the reflectivity by 12 percent and the quantum efficiency by 8.

Boundaries at the minima of the planar-averaged rho(z), with the surface where
rho(z) falls below a fraction of its maximum, is one criterion for both.  It
tiles by construction, it handles rumpling, molecules and interfaces the same
way, and it gives volumes that mean something.

This is deliberately a standalone script with no dependency beyond numpy, so it
works for anyone with a single CASTEP run.  Wrap it in an AiiDA calcfunction if
you want provenance over a campaign; do not make AiiDA a requirement for using
it.

Usage
-----
    layer_boundaries.py <seed>                 looks for <seed>.den_fmt and
                                               <seed>-out.cell in the cwd
    layer_boundaries.py <seed> --surface 0.01  vacuum threshold, fraction of max
    layer_boundaries.py <seed> --plot out.png  overlay for checking by eye
"""
import argparse
import os
import sys

import numpy as np


def read_den_fmt(path):
    """Planar-averaged charge density.

    Returns (z, rho, c) with z in Angstrom and rho in arbitrary units -- only
    ratios and the positions of the minima are used, so the normalisation of the
    CASTEP column ("electrons/grid_point * number of grid_points") never has to
    be unpicked.
    """
    header = []
    with open(path) as handle:
        for _ in range(12):
            line = next(handle)
            header.append(line.strip().split())
            if 'END header:' in line:
                break
    spin = int(header[7][0])
    c = float(header[5][5])
    n_a, n_b, n_c = (int(x) for x in header[8][:3])
    per_plane = n_a * n_b * spin

    data = np.genfromtxt(path, skip_header=11)
    data = data[data[:, 2].argsort()]          # by the c grid index
    column = data[:, 3]
    planes = np.add.reduceat(column, range(0, len(column), per_plane)) / per_plane
    z = (np.arange(len(planes)) + 1) * (c / (n_c + 1))
    return z, planes, c


def atom_z_from_cell(path):
    """Cartesian z of every atom, in Angstrom, from a CASTEP -out.cell."""
    lattice, frac, cart = [], [], []
    block = None
    for line in open(path):
        low = line.lower().strip()
        if low.startswith('%block'):
            block = low.split()[-1]
            continue
        if low.startswith('%endblock'):
            block = None
            continue
        if block == 'lattice_cart' and len(line.split()) == 3:
            lattice.append([float(x) for x in line.split()])
        elif block == 'positions_frac' and len(line.split()) >= 4:
            frac.append(float(line.split()[3]))
        elif block == 'positions_abs' and len(line.split()) >= 4:
            cart.append(float(line.split()[3]))
    if cart:
        return np.sort(np.array(cart))
    if not frac or len(lattice) < 3:
        raise ValueError('%s: no positions found' % path)
    return np.sort(np.array(frac) * lattice[2][2])


def group_layers(atom_z):
    """Group atom z into layers, the way analyse_geometry does.

    The tolerance is half the median of the "large" gaps, where large means above
    10 percent of the biggest, so intra-layer rumpling does not drag the median
    towards zero.  Mirrored from the Fortran on purpose: the point of this script
    is to feed that code, not to disagree with it about what a layer is.
    """
    gaps = np.diff(atom_z)
    if len(gaps) == 0:
        return [np.array(atom_z)]
    large = gaps[gaps > 0.1 * gaps.max()]
    tol = 0.5 * float(np.median(large))
    layers, current = [], [atom_z[0]]
    for gap, z in zip(gaps, atom_z[1:]):
        if gap > tol:
            layers.append(np.array(current))
            current = []
        current.append(z)
    layers.append(np.array(current))
    return layers[::-1]                        # surface first, as OptaDOS orders


def minimum_between(z, rho, lower, upper):
    """z of the lowest planar-averaged density strictly between two planes."""
    mask = (z > lower) & (z < upper)
    if not mask.any():
        return 0.5 * (lower + upper)           # grid too coarse to resolve it
    window = np.where(mask)[0]
    return float(z[window[np.argmin(rho[window])]])


def surface_plane(z, rho, top_atom, fraction):
    """First plane above the outermost atoms where rho falls below the threshold.

    Same criterion people already apply by hand to choose photo_slab_max, which
    is the point: the surface and the interior boundaries stop being defined by
    two unrelated rules.
    """
    threshold = fraction * rho.max()
    above = np.where(z > top_atom)[0]
    for i in above:
        if rho[i] < threshold:
            return float(z[i])
    raise ValueError('density never falls below %g of its maximum above the '
                     'topmost atom -- too little vacuum, or raise --surface'
                     % fraction)


def analyse(seed, fraction=0.01, root='.'):
    den = os.path.join(root, seed + '.den_fmt')
    cell = os.path.join(root, seed + '-out.cell')
    z, rho, c = read_den_fmt(den)
    atom_z = atom_z_from_cell(cell)
    layers = group_layers(atom_z)
    centroids = [float(l.mean()) for l in layers]

    z_top = surface_plane(z, rho, atom_z.max(), fraction)

    interior = [minimum_between(z, rho, centroids[i + 1], centroids[i])
                for i in range(len(centroids) - 1)]

    n_boxes = (len(centroids) + 1) // 2
    return dict(seed=seed, z=z, rho=rho, c=c, atom_z=atom_z,
                centroids=centroids, z_top=z_top,
                interior=interior, n_boxes=n_boxes,
                z_middle=interior[n_boxes - 1])


def keyword_block(result):
    tops = [result['z_top']] + result['interior'][:result['n_boxes'] - 1]
    lines = ['! generated by layer_boundaries.py from %s.den_fmt' % result['seed'],
             '! boundaries at the minima of the planar-averaged charge density,',
             '! surface where it falls below the chosen fraction of its maximum.',
             '! photo_slab_max is the first entry and must not be set separately.',
             '%BLOCK photo_layers_tops']
    lines += ['  %12.5f' % t for t in tops]
    lines += ['%ENDBLOCK photo_layers_tops',
              '',
              'photo_slab_middle    : %12.5f' % result['z_middle']]
    return '\n'.join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('seed')
    ap.add_argument('--surface', type=float, default=0.01,
                    help='vacuum threshold as a fraction of max rho (default 0.01)')
    ap.add_argument('--root', default='.', help='directory holding the files')
    ap.add_argument('--plot', help='write an overlay of rho(z) and the boundaries')
    ap.add_argument('--compare', action='store_true',
                    help='also show the centroid midpoints OptaDOS would use')
    args = ap.parse_args(argv)

    r = analyse(args.seed, args.surface, args.root)
    print(keyword_block(r))

    if args.compare:
        mid = [0.5 * (r['centroids'][i] + r['centroids'][i + 1])
               for i in range(len(r['centroids']) - 1)]
        print('\n! interior boundary   density min   centroid midpoint   shift')
        for i, (d, m) in enumerate(zip(r['interior'], mid), 1):
            print('!   %2d %18.5f %15.5f %13.5f' % (i, d, m, d - m))
        print('! surface %25.5f' % r['z_top'])

    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.plot(r['z'], r['rho'], lw=1.2, label=r'planar-averaged $\rho(z)$')
        for b in r['interior']:
            ax.axvline(b, color='0.6', lw=0.8)
        ax.axvline(r['z_top'], color='C3', lw=1.4, label='surface')
        ax.axvline(r['z_middle'], color='C2', lw=1.4, label='explicit/bulk')
        ax.set_yscale('log')
        floor = max(r['rho'].min(), 1e-12)
        ax.plot(r['atom_z'], np.full_like(r['atom_z'], floor),
                marker='|', ms=12, ls='none', color='C1', label='atoms')
        ax.set_xlabel(r'$z$ ($\mathrm{\AA}$)')
        ax.set_ylabel(r'$\rho(z)$ (arb.)')
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print('\n! wrote %s' % args.plot)
    return 0


if __name__ == '__main__':
    sys.exit(main())
