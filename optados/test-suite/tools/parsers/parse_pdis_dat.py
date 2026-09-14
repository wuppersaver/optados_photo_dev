"""
Parser function parse() for the .pdis.dat file of a pdispersion run.

The file is a header of '#' lines naming the projectors, then one block per
k-point: a 'K-point   N   kx ky kz' line and one row per band, the band energy
followed by the weight of each projector. Returns the k-point coordinates, every
band energy and every weight, in file order, so a k-point that moved or a band
that changed shows up as a mismatch.
"""
from __future__ import print_function

import inspect
from collections import defaultdict

from . import show_output


def parse(fname):
    """
    Open the file, parse it and return the values
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        for line in f:
            fields = line.split()
            if not fields or line.lstrip().startswith('#'):
                continue
            if fields[0] == 'K-point':
                retdict['kpoint_index'].append(int(fields[1]))
                retdict['kpoint'].extend(float(x) for x in fields[2:5])
                continue
            retdict['energy'].append(float(fields[0]))
            retdict['weight'].extend(float(x) for x in fields[1:])

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-" * 72)
    return retdict
