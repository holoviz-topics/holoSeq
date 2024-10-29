# Proof of cocept data are Arima HiC reads from the Arctic Ground Squirrel mUroPar1 VGP genomeArk repository
# processed with Dephine's Pretext workflow using Bellerophon to remove chimeric reads
# The paired bam output is converted to PAF with an awk script (!) and that's what is read
# in this code.
# The pairs are parsed to extract the haplotype designator from the contig name
# typically a suffix like H1 extracted in getHap - rewrite that to suit your names.
# Contig ordering really matters for the plots to make any sense.
# Ideally, curators name them so they sort alphanumerically without effort.
# There's a sorthapqname function that is used here. Designed for the VGP data
# Will need to be replaced for other naming conventions.
# It's a mess.
# Sorting by contig name is based on VGP conventions - SUPER_ first, then scaffolds
# One problem to watch out for is that any differences in ordering of the X and Y contigs can make all sorts of
# artifacts appear such as the kaleidoscopic patterns seen in Pretextviewer.
# Ross Lazarus October 2024

from collections import OrderedDict
from functools import cmp_to_key
import gzip
import math
import numpy as np
import os
import sys

# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/huge.paf"

holoSeqHeaders = ["@v1HoloSeq1D", "@v1HoloSeq2D"]


def rotatecoords(x, y, radians=0.7853981633974483, origin=(0, 0)):
    # https://gist.github.com/LyleScott/d17e9d314fbe6fc29767d8c5c029c362
    offset_x, offset_y = origin
    adjusted_x = x - offset_x
    adjusted_y = y - offset_y
    cos_rad = math.cos(radians)
    sin_rad = math.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y
    return qx, qy


def getHap(contig):
    """
    function to return suffix H1 from chrH1 - adjust to suit.
    """
    return contig[-2:]


def sorthapqname(s1, s2):
    """
    fugly hack to sort super contigs before anything else
    then by contig number or if they are the same offset
    ('SUPER_2H1', 226668729) , ('SUPER_1H1', 284260672), ('SUPER13_unloc_5H1',..), (Scaffold_1116H2, ...)
    """
    if s1[0] == s2[0]:  # simplest case - same contig, sort on offset
        return s1[1] - s2[1]  # neg if left sorts before
    s11, s12 = s1[0].split("_", 1)
    s1n = s12.split("_")[-1][:-2]
    s1_super = (s11.upper() == "SUPER") or (s11.upper().startswith("CHR"))
    s21, s22 = s2[0].split("_", 1)
    s2n = s22.split("_")[-1][:-2]
    s2_super = (s21.upper() == "SUPER") or (s21.upper().startswith("CHR"))
    if s1n.isdigit():
        s1nn = int(s1n)
    else:
        s1nn = ord(s1n[0]) * 1000
    if s2n.isdigit():
        s2nn = int(s2n)
    else:
        s2nn = ord(s2n[0]) * 1000
    if s1_super == s2_super:
        nunder1 = len(s1[0].split("_"))
        nunder2 = len(s2[0].split("_"))  # _unloc or whatever
        if nunder1 == nunder2:
            return s1nn - s2nn
        else:
            return nunder1 - nunder2
    elif s1_super:
        return -1
    elif s2_super:
        return 1
    else:
        return s1nn - s2nn


def export_mapping(hsId, outFileName, haps, hnames, hlens, x, y, anno):
    """
    @v1HoloSeq2D for example
    """

    def prepHeader(haps, hnames, hlens, hsId):
        """
        holoSeq output format
        """
        h = ["%s%s %s %d" % ('@', haps[i], hnames[i], hlens[i]) for i in range(len(hlens))]
        h.insert(0, hsId)
        return h

    hdr = prepHeader(haps, hnames, hlens, hsId)

    print('haps =', haps)
    print('header=', hdr)
    with gzip.open(outFileName, mode='wb') as ofn:
        ofn.write(str.encode('\n'.join(hdr) + '\n'))
        if len(anno) == len(x):
            for i in range(len(x)):
                row = str.encode("%d %d %s\n" % (x[i], y[i], anno[i]))
                ofn.write(row)
        else:
            for i in range(len(x)):
                row = str.encode("%d %d\n" % (x[i], y[i]))
                ofn.write(row)

inFile = sys.argv[1]
print('inFile=', inFile)
hlstarts = OrderedDict()
hqstarts = OrderedDict()
hlens = {}
haps = []
hlsorts = {}
hqsorts = {}
with open(inFile, "r") as f:
    for i, rows in enumerate(f):
        row = rows.split()
        c1 = row[0]
        c2 = row[5]
        for hap in [getHap(c1), getHap(c2)]:
            if hap not in haps:
                haps.append(hap)
                hlens[hap] = {}
                hlsorts[hap] = []
                hqsorts[hap] = []
        hp = getHap(c1)   
        if not hlens[hp].get(c1, None):
            hlens[hp][c1] = int(row[1])
            hlsorts[hp].append((int(row[1]), c1))
            hqsorts[hp].append((c1, int(row[1])))
        hp = getHap(c2)   
        if not hlens[hp].get(c2, None):
            hlens[hp][c2] = int(row[6])
            hlsorts[hap].append((int(row[6]), c2))
            hqsorts[hap].append((c2, int(row[6])))
print(haps)
for hap in haps:
    cum = 1
    hlsorts[hap].sort(reverse=True)
    hqsorts[hap].sort(key=cmp_to_key(sorthapqname))
    hlstarts[hap] = OrderedDict()
    hqstarts[hap] = OrderedDict()
    for clen, contig in hlsorts[hap]:
        hlstarts[hap][contig] = cum
        cum += clen
    cum = 1
    for contig, clen in hqsorts[hap]:
        hqstarts[hap][contig] = cum
        cum += clen
print('hqstarts=',hqstarts)
h1starts = [hqstarts[haps[0]][x] for x in hqstarts[haps[0]].keys()]
h1names = list(hqstarts[haps[0]].keys())
if len(haps) > 1:
    h2starts = [hqstarts[haps[1]][x] for x in hqstarts[haps[1]].keys()]
    h2names = list(hqstarts[haps[1]].keys())
# have the axes set up so prepare the three plot x/y vectors
# for a second pass to calculate all the coordinates.
# adding tooltips just does not scale so abando - see the tooltip old version
cis1 = {"x": [], "y": []}
cis2 = {"x": [], "y": []}
trans1 = {"x": [], "y": []}
with open(inFile, "r") as f:
    for rows in f:
        row = rows.split()
        c1 = row[0]
        c2 = row[5]
        H1 = getHap(c1)
        H2 = getHap(c2)
        if H1 != H2:  # trans
            if H1 == haps[0]:  # x is h1 for trans - otherwise ignore
                trans1["x"].append(hqstarts[H1][c1] + int(row[2]))
                trans1["y"].append(hqstarts[H2][c2] + int(row[7]))
            else:
                trans1["y"].append(hqstarts[H1][c1] + int(row[2]))
                trans1["x"].append(hqstarts[H2][c2] + int(row[7]))
        else:  # cis
            if H1 == haps[0]:
                cis1["x"].append(hqstarts[H1][c1] + int(row[2]))
                cis1["y"].append(hqstarts[H1][c2] + int(row[7]))
            else:
                cis2["x"].append(hqstarts[H1][c1] + int(row[2]))
                cis2["y"].append(hqstarts[H1][c2] + int(row[7]))
outPrefix = os.path.basename(inFile)
hap = haps[0]
export_mapping(
    holoSeqHeaders[1],
    "%s_cis1.hseq.gz" % outPrefix,
    [hap for i in range(len(h1names))],
    h1names,
    h1starts,
    cis1["x"],
    cis1["y"],
    [],
)
