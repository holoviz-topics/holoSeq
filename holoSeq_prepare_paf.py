# for Mashmap paf, python holoSeq_prepare_paf.py --inFile  hg002_2k99.paf --title "hg002 Mashmap" --hap_indicator None --contig_sort length
# works ok!
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

import argparse

from collections import OrderedDict
from functools import cmp_to_key
import gzip
import math
import numpy as np
import os

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
     help="None, Suffix (H[1,2]) Dashsuffix (_H...)"
    """
    hi = args.hap_indicator
    if hi == "None":
        return('H1')
    elif hi == "Suffix":
        return contig[-2:]
    elif hi == "Dashsuffix":
        return contig.split('_')[-1]


def VGPsortfunc(s1, s2):
    """
    fugly hack to sort super contigs before anything else
    then by contig number or if they are the same offset
    ('SUPER_2H1', 226668729) , ('SUPER_1H1', 284260672), ('SUPER13_unloc_5H1',..), (Scaffold_1116H2, ...)
    or chr10_H1 etc
    """
    if s1[0] == s2[0]:  # simplest case - same contig, sort on offset
        return s1[1] - s2[1]  # neg if left sorts before
    sorts = [{},{}]
    for i, (contig, length) in enumerate([s1,s2]):
        if '_' in contig:
            conta, contb = [x.upper() for x in contig.split("_", 1)]
        else:
            conta = contig.upper()
            contb = ''
        if conta.startswith('CHR'):
            contign = conta.replace("CHR","").replace("_","")
        elif conta.startswith('SUPER') or conta.startswith('SCAFFOLD') and contb > '':
            if "UNLOC" in contig.upper():
                unloc = True
                contign = conta.split("_")[0]
            else:
                unloc = False
                contign = contb
        else:
            contign = contb.split("_")[-1][:-2]
        if unloc:
            contign = contb.split("_")[0]
        if contign.isdigit():
            nval = int(contign)
        else:            
            nval = ord(contign[0])
        is_super = not unloc and (conta in ["SUPER", "CHR", 'SCAFFOLD'])
        sorts[i]['is_super'] = is_super
        sorts[i]['n'] = nval
        sorts[i]['conta'] = conta
    if sorts[0]['conta'] == sorts[1]['conta']:
        return sorts[0]['n']  - sorts[1]['n']
    elif sorts[0]['conta'] in ['SUPER', 'CHR']:
        return -1
    elif sorts[1]['conta'] in ['SUPER', 'CHR']:
        return 1
    else:
        nunder1 = len(sorts[0]['conta'].split("_"))
        nunder2 = len(sorts[1]['conta'].split("_"))  # _unloc or whatever
        if nunder1 == nunder2:
            return sorts[0]['n']  - sorts[1]['n'] 
        else:
            return nunder1 - nunder2
    


def Lengthsortfunc(s1, s2):
    """
   
    """
    return s1[1] - s2[1]  # neg if left sorts before

def export_mapping(hsId, outFileName, haps, hnames, hlens, x, y, anno, title):
    """
    @v1HoloSeq2D for example
    """

    def prepHeader(haps, hnames, hlens, hsId, title):
        """
        holoSeq output format
        """
        h = ["@%s %s %d" % (haps[i], hnames[i], hlens[i]) for i in range(len(hlens))]
        h.insert(0, '@title %s' % title)
        h.insert(0, hsId)
        return h

    hdr = prepHeader(haps, hnames, hlens, hsId, title)

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

parser = argparse.ArgumentParser(description="", epilog="")
parser.add_argument("--inFile", help="PAF with paired alignments", default="mUroPar1.paf")
parser.add_argument("--title", help="Title for the plot", default="Plot title goes here")
parser.add_argument("--contig_sort", help="VGPname, name, length", default="VGPname")
parser.add_argument("--hap_indicator", help="None, Suffix (H[1,2]) Dashsuffix (_H...)", default="None")
parser.add_argument("--version", "-V", action="version", version='0.1')
args = parser.parse_args()
inFile = args.inFile
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
for hap in haps:
    cum = 1
    if args.contig_sort == "VGPname":
        hqsorts[hap].sort(key=cmp_to_key(VGPsortfunc))
    elif args.contig_sort == "name":
        hqsorts[hap].sort()
    hlsorts[hap].sort(reverse=True)
    hlstarts[hap] = OrderedDict()
    hqstarts[hap] = OrderedDict()
    for clen, contig in hlsorts[hap]:
        hlstarts[hap][contig] = cum
        cum += clen
    cum = 1
    for contig, clen in hqsorts[hap]:
        hqstarts[hap][contig] = cum
        cum += clen
if args.contig_sort == "length":
    starts = hlstarts
else:
    starts = hqstarts
h1starts = [starts[haps[0]][x] for x in starts[haps[0]].keys()]
h1names = list(starts[haps[0]].keys())
print('h1names=',h1names[:30])
if len(haps) > 1:
    if args.contig_sort == "length":
        starts = hlstarts
    else:
        starts = hqstarts
    h2starts = [starts[haps[1]][x] for x in starts[haps[1]].keys()]
    h2names = list(starts[haps[1]].keys())
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
    args.title
)
