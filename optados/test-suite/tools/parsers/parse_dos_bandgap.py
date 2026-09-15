"""
Parser function parse() for the Fermi energy and the Bandgap Analysis block of
an OptaDOS .odo file.

The band gap is computed from a gather of every node's k-points onto the root
node, so it is the part of a DOS run that depends on the MPI decomposition. The
k-point coordinates of the VBM and CBM are extracted as well as the gaps: a
gather that put the right values at the wrong k-points would leave the gaps
alone and move these.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

e_fermi_ab = re.compile(r"Fermi\ energy\ \(Adaptive\ broadening\)\ \:\s*([0-9\.-]+)\s*")
e_fermi_insulator = re.compile(r"Fermi energy assuming insulator :\s*([0-9\.-]+)")
thermal = re.compile(r"Thermal Bandgap :\s*([0-9\.-]+)")
optical = re.compile(r"Spin :\s*(\d+)\s*:\s*([0-9\.-]+)\s*eV\s*<- OBg")
average = re.compile(r"Spin :\s*(\d+)\s*:\s*([0-9\.-]+)\s*eV\s*<- ABg")
weighted = re.compile(r"Weighted Average :\s*([0-9\.-]+)")
kpoint = re.compile(r"(At kpoint|Between VBM kpoint|and CBM kpoint)\s*:\s*"
                    r"([0-9\.-]+)\s+([0-9\.-]+)\s+([0-9\.-]+)")
count = re.compile(r"Spin :\s*(\d+)\s*:\s*(\d+)\s+(\d+)\s*\|")


def parse(fname):
    """
    Open the file, parse it and return the values
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    for l in lines:
        m = e_fermi_ab.search(l)
        if m:
            retdict["fermi_ab"].append(float(m.group(1)))
            continue
        m = e_fermi_insulator.search(l)
        if m:
            retdict["fermi_insulator"].append(float(m.group(1)))
            continue
        m = thermal.search(l)
        if m:
            retdict["thermal_gap"].append(float(m.group(1)))
            continue
        m = optical.search(l)
        if m:
            retdict["optical_gap"].append(float(m.group(2)))
            continue
        m = average.search(l)
        if m:
            retdict["average_gap"].append(float(m.group(2)))
            continue
        m = weighted.search(l)
        if m:
            retdict["weighted_average_gap"].append(float(m.group(1)))
            continue
        m = kpoint.search(l)
        if m:
            retdict["gap_kpoint"].extend(float(x) for x in m.groups()[1:])
            continue
        m = count.search(l)
        if m:
            retdict["vbm_cbm_multiplicity"].extend([int(m.group(2)), int(m.group(3))])
            continue

    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
