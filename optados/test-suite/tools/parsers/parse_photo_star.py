"""
Parser function parse() to pull the symmetry-star bookkeeping out of a .odo.

Item 10 of documents/Photoemission_setup_invariance_fixes.md: CASTEP reduces the
k-mesh with the point group AND time reversal, but writes only the point group
into symmetry_ops. photo_check_star establishes whether the operation list can
reach the whole star each k-point weight describes, and reports it near the top
of the .odo.

Deliberately a parser of its own rather than an addition to parse_photo_odo.
Extending that one would give every existing photoemission benchmark new fields
and force a re-bless of all thirty-four, for a quantity none of them is about.

The star report is written before any matrix element is read, so a run that goes
on to fail for want of an ome_bin still carries a complete and correct block.
That is what lets the tests here ship a fixture of nothing but a .bands and an
-out.cell.

The box is framed with ``|`` when nothing is wrong and with ``!`` when it carries
a warning, so every pattern here accepts either. Matching only ``|`` would make
exactly the interesting cases parse as empty.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

n_ops = re.compile(r"[|!]\s*Symmetry operations in the cell\s*:\s*(\d+)\s*[|!]")
n_grid = re.compile(r"[|!]\s*In-plane MP grid points \([^)]*\)\s*:\s*(\d+)\s*[|!]")
n_tr = re.compile(r"[|!]\s*k-points needing time-reversed partners\s*:\s*(\d+)\s*[|!]")
tr_weight = re.compile(r"[|!]\s*\.\.\. carrying a k-point weight of\s*:\s*([\d.]+)\s*[|!]")
disabled = re.compile(r"[|!]\s*Nothing has been changed")
bad_int = re.compile(r"[|!]\s*k-points with an unexpected weight\s*:\s*(\d+)\s*[|!]")
mesh_open = re.compile(r"[|!]\s*k-points the operations move off the mesh\s*:\s*(\d+)\s*[|!]")
bad_ratio = re.compile(r"[|!]\s*k-points with an unexpected image count\s*:\s*(\d+)\s*[|!]")


def parse(fname):
    """
    Open the file, parses it and return the values
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    for lno, l in enumerate(lines):

        match = n_ops.search(l)
        if match:
            retdict['star_n_ops'].append(int(match.groups()[0]))
            continue

        match = n_grid.search(l)
        if match:
            retdict['star_n_grid'].append(int(match.groups()[0]))
            continue

        match = n_tr.search(l)
        if match:
            retdict['star_n_tr'].append(int(match.groups()[0]))
            continue

        match = tr_weight.search(l)
        if match:
            retdict['star_tr_weight'].append(float(match.groups()[0]))
            continue

        match = bad_int.search(l)
        if match:
            retdict['star_bad_integer'].append(int(match.groups()[0]))
            continue

        match = mesh_open.search(l)
        if match:
            retdict['star_mesh_not_invariant'].append(int(match.groups()[0]))
            continue

        match = bad_ratio.search(l)
        if match:
            retdict['star_bad_ratio'].append(int(match.groups()[0]))
            continue

        # "Nothing needed time reversal" and "we could not tell" must not look
        # the same to testcode. The first reports star_n_tr = 0; the second
        # reports no count at all, and is caught here instead.
        if disabled.search(l):
            retdict['star_disabled'].append(1)
            continue
        ###############################################################

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)

    return retdict
