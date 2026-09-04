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

Boundaries at the minima of the planar-averaged rho(z), with the surface placed
where **0.5 percent of one layer's charge remains above it**, is one criterion for
both.  It tiles by construction, it handles rumpling, molecules and interfaces
the same way, and it gives volumes that mean something.

The 0.5 percent is a convention, not a derivation -- the vacuum tail is an
exponential, so any level-set criterion moves the plane by about 0.9 Ang per
decade of the fraction.  It is the right convention because it reproduces the
rule of thumb it replaces, to within 0.07 Ang on every Cu slab tested, while
being stated in a form that transfers to a material whose density maximum sits
somewhere else.

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

    Kept as an alternative to the charge criterion, not as the default: a
    fraction of ``rho.max()`` is not transferable, because the maximum sits
    wherever the densest species in the cell happens to be. On a heterostructure
    the same fraction therefore means different things in different parts of the
    same slab.
    """
    threshold = fraction * rho.max()
    above = np.where(z > top_atom)[0]
    for i in above:
        if rho[i] < threshold:
            return float(z[i])
    raise ValueError('density never falls below %g of its maximum above the '
                     'topmost atom -- too little vacuum, or raise --surface'
                     % fraction)


def charge_per_layer(z, rho, centroids):
    """Charge in one interior layer, as a partition rather than an attribution.

    Taken between the two density minima bracketing a layer near the middle of
    the slab. Two reasons for the choice. It has to be an *interior* layer,
    because a surface layer genuinely holds different charge -- the spill-out is
    the very thing being measured at the other end. And it has to be a partition:
    charge cannot be attributed uniquely to a layer, since neighbouring layers
    overlap, but every electron does belong to exactly one box, and this is the
    same partition the boxes themselves use.
    """
    if len(centroids) < 3:
        raise ValueError('need at least three layers to take an interior one')
    i = len(centroids) // 2
    upper = minimum_between(z, rho, centroids[i], centroids[i - 1])
    lower = minimum_between(z, rho, centroids[i + 1], centroids[i])
    window = (z > lower) & (z < upper)
    return float(rho[window].sum() * (z[1] - z[0]))


def surface_plane_charge(z, rho, top_atom, fraction, q_layer):
    """Plane above which a given fraction of one layer's charge remains.

    The default criterion. It does **not** pin the surface down any harder than
    a density threshold does -- the vacuum tail is an exponential, so
    ``Q(z) = integral of rho from z ~ lambda*rho(z)`` carries the same
    z-dependence, and either criterion moves the plane by lambda*ln(10), about
    0.9 Ang, per decade of the fraction. What it buys is *meaning*: "half a
    percent of a layer's electrons lie outside the box" is a statement that can
    be defended and compared between materials, and it is local, so it survives
    a heterostructure where a fraction of the global rho.max() does not.

    0.5 % is the default because it reproduces the convention it replaces. The
    charge lying above "outermost atom plus one bulk repeat" measures 0.32 to
    0.41 % of a layer on the Cu slabs, and conversely 0.5 % lands 1.73 to 1.78 Ang
    above the outermost atoms against a 1.782 Ang repeat. The old rule of thumb
    was a good one; this states it in a form that transfers.

    Contamination from deeper layers is negligible here: layer 2 contributes
    exp(-repeat/lambda), about 1 %, of layer 1's tail, and layer 3 a hundredth of
    that, so a 1 % correction to a quantity set at 0.5 % moves the plane by a few
    thousandths of an Angstrom.
    """
    target = fraction * q_layer
    dz = z[1] - z[0]
    above = np.cumsum(rho[::-1])[::-1] * dz
    candidates = np.where((z > top_atom) & (above < target))[0]
    if not len(candidates):
        raise ValueError('the charge above the topmost atom never falls below '
                         '%g of a layer -- too little vacuum, or raise the '
                         'fraction' % fraction)
    return float(z[candidates[0]])


def analyse(seed, fraction=0.005, root='.', criterion='charge'):
    den = os.path.join(root, seed + '.den_fmt')
    cell = os.path.join(root, seed + '-out.cell')
    z, rho, c = read_den_fmt(den)
    atom_z = atom_z_from_cell(cell)
    layers = group_layers(atom_z)
    centroids = [float(l.mean()) for l in layers]

    if criterion == 'charge':
        q_layer = charge_per_layer(z, rho, centroids)
        z_top = surface_plane_charge(z, rho, atom_z.max(), fraction, q_layer)
    else:
        q_layer = None
        z_top = surface_plane(z, rho, atom_z.max(), fraction)

    interior = [minimum_between(z, rho, centroids[i + 1], centroids[i])
                for i in range(len(centroids) - 1)]

    n_boxes = (len(centroids) + 1) // 2
    return dict(seed=seed, z=z, rho=rho, c=c, atom_z=atom_z,
                centroids=centroids, z_top=z_top, q_layer=q_layer,
                criterion=criterion, fraction=fraction,
                interior=interior, n_boxes=n_boxes,
                z_middle=interior[n_boxes - 1])


def keyword_block(result):
    tops = [result['z_top']] + result['interior'][:result['n_boxes'] - 1]
    if result['criterion'] == 'charge':
        how = ('! surface where %.3f %% of one layer\'s charge remains above it.'
               % (100 * result['fraction']))
    else:
        how = ('! surface where rho falls below %.3f %% of its maximum.'
               % (100 * result['fraction']))
    lines = ['! generated by layer_boundaries.py from %s.den_fmt' % result['seed'],
             '! boundaries at the minima of the planar-averaged charge density,',
             how,
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
    ap.add_argument('--surface', type=float, default=0.005,
                    help='fraction of one layer\'s charge left above the surface '
                         'plane (default 0.005); with --criterion density it is '
                         'instead a fraction of max rho')
    ap.add_argument('--criterion', choices=('charge', 'density'), default='charge',
                    help='how the surface plane is placed (default charge)')
    ap.add_argument('--root', default='.', help='directory holding the files')
    ap.add_argument('--plot', help='write an overlay of rho(z) and the boundaries')
    ap.add_argument('--compare', action='store_true',
                    help='also show the centroid midpoints OptaDOS would use')
    args = ap.parse_args(argv)

    r = analyse(args.seed, args.surface, args.root, args.criterion)
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
