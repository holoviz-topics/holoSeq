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
from pathlib import Path

import gzip
import math
import numpy as np
import os

import pybigtools

from holoseq_gff import parseGFF

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
        return "H1"
    elif hi == "Suffix":
        return contig[-2:]
    elif hi == "Dashsuffix":
        return contig.split("_")[-1]


def VGPsortfunc(s1, s2):
    """
    fugly hack to sort super contigs before anything else
    then by contig number or if they are the same offset
    ('SUPER_2H1', 226668729) , ('SUPER_1H1', 284260672), ('SUPER13_unloc_5H1',..), (Scaffold_1116H2, ...)
    or chr10_H1 etc
    """
    if s1[0] == s2[0]:  # simplest case - same contig, sort on offset
        return s1[1] - s2[1]  # neg if left sorts before
    sorts = [{}, {}]
    for i, (contig, length) in enumerate([s1, s2]):
        if "_" in contig:
            conta, contb = [x.upper() for x in contig.split("_", 1)]
        else:
            conta = contig.upper()
            contb = ""
        if conta.startswith("CHR"):
            contign = conta.replace("CHR", "").replace("_", "")
        elif conta.startswith("SUPER") or conta.startswith("SCAFFOLD") and contb > "":
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
        is_super = not unloc and (conta in ["SUPER", "CHR", "SCAFFOLD"])
        sorts[i]["is_super"] = is_super
        sorts[i]["n"] = nval
        sorts[i]["conta"] = conta
    if sorts[0]["conta"] == sorts[1]["conta"]:
        return sorts[0]["n"] - sorts[1]["n"]
    elif sorts[0]["conta"] in ["SUPER", "CHR"]:
        return -1
    elif sorts[1]["conta"] in ["SUPER", "CHR"]:
        return 1
    else:
        nunder1 = len(sorts[0]["conta"].split("_"))
        nunder2 = len(sorts[1]["conta"].split("_"))  # _unloc or whatever
        if nunder1 == nunder2:
            return sorts[0]["n"] - sorts[1]["n"]
        else:
            return nunder1 - nunder2


def Lengthsortfunc(s1, s2):
    """ """
    return s1[1] - s2[1]  # neg if left sorts before


class gffConvert:
    """
        gff = {"start": [], "end": [], "strand": [], "type": [], "qualifiers": []}
    """
    def __init__(self, gff):
        dgff = parseGFF(gff)
        ftypes = np.unique(gff['type'])
        df = pd.DataFrame.from_dict(dgff)
    

class bwConvert:
    """
    points = [(50*i, 100+random.random()) for i in range(10000)]
    hv.Curve(points).opts(interpolation='steps-post').opts(width=1000)
    bigwig to barchart - will autoscale to line?
    """

    def __init__(self, inFname, outFname, args):
        self.inFname = inFname
        fakepath = "in.bw"
        if os.path.isfile(fakepath):
            os.remove(fakepath)
        p = Path(fakepath)
        p.symlink_to(inFname)  # required by pybigtools (!)
        bwf = pybigtools.open(fakepath)
        chrlist = bwf.chroms()
        hnames = list(chrlist.keys())
        hlens = [chrlist[x] for x in hnames]
        hstarts = [0]
        cum = 0
        for clen in hlens:
            cum += clen
            hstarts.append(cum)
        data = {}
        for i,chr in enumerate(hnames):
            cstart = hstarts[i]
            data[chr] = {}
            bw = bwf.records(chr)
            # Return the records of a given range on a chromosome. The result is an iterator of tuples. For BigWigs, these tuples are in the format (start: int, end: int, value: float).
            data[chr]["xstart"] = [x[0] + cstart for x in bw]
            bw = bwf.records(chr) 
            data[chr]["xend"] = [x[1] + cstart for x in bw]
            bw = bwf.records(chr) 
            data[chr]["xval"] = [x[2] for x in bw]
        self.export_mapping("@v1HoloSeq1D line", outFname, hnames, hstarts, data, args)

    def export_mapping(self, hsId, outFname, hnames, hstarts, data, args):
        """
        for bigwig
        @v1HoloSeq2D for example
        """

        def prepHeader(hnames, hstarts, hsId, args):
            """
            holoSeq output format
            """
            h = ["@%s %d" % (hnames[i], hstarts[i]) for i in range(len(hnames))]
            metah = [hsId, "@@title %s" % args.title, "@@datasource %s" % "bigwig", "@@datafile %s" % self.inFname, "@@refURI %s" % args.refURI]
            
            return metah + h

        hdr = prepHeader(hnames, hstarts, hsId, args)

        with gzip.open(outFname, mode="wb") as ofn:
            ofn.write(str.encode("\n".join(hdr) + "\n"))
            for chr in data.keys():
                for i in range(len(data[chr]["xstart"])):
                    row = str.encode(
                        "%d %d\n" % (data[chr]["xstart"][i], data[chr]["xval"][i])
                    )
                    ofn.write(row)


class pafConvert:
    """
    paf to xy and axis metadata
    """

    def __init__(self, inFname, args):
        self.inFname = inFname
        hlstarts = {}
        hqstarts = {}
        hlens = {}
        haps = []
        hlsorts = {}
        hqsorts = {}
        with open(inFname) as f:
            for i, rows in enumerate(f):
                row = rows.strip().split()
                if len(row) > 7:
                    c1 = row[0]
                    c2 = row[5]
                    hap = getHap(c1)
                    if not hap in haps:
                        haps.append(hap)
                        hlens[hap] = {}
                        hlsorts[hap] = []
                        hqsorts[hap] = []
                        print('added hap', hap)
                    if not hlens[hap].get(c1, None):
                        hlens[hap][c1] = int(row[1])
                        hlsorts[hap].append((int(row[1]), c1))
                        hqsorts[hap].append((c1, int(row[1])))
                    hap = getHap(c2)
                    if not hap in haps:
                        haps.append(hap)
                        hlens[hap] = {}
                        hlsorts[hap] = []
                        hqsorts[hap] = []
                        print('added hap', hap)
                    if not hlens[hap].get(c2, None):
                        print('row=',row, "added", c2)
                        hlens[hap][c2] = int(row[6])
                        hlsorts[hap].append((int(row[6]), c2))
                        hqsorts[hap].append((c2, int(row[6])))
        for hap in haps:
            cum = 0
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
            cum = 0
            for contig, clen in hqsorts[hap]:
                hqstarts[hap][contig] = cum
                cum += clen
        if args.contig_sort == "length":
            starts = hlstarts
        else:
            starts = hqstarts
        h1starts = [starts[haps[0]][x] for x in starts[haps[0]].keys()]
        h1names = list(starts[haps[0]].keys())
        print("h1names=", h1names[:30])
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
        print('hqstarts=',hqstarts)
        if args.contig_sort == "length":
            starts = hlstarts
        else:
            starts = hqstarts
        rowi = 0
        with open(inFname) as f:
            for rows in f:
                rowi += 1
                row = rows.strip().split()
                if len(row) > 7:
                    c1 = row[0]
                    c2 = row[5]
                    H1 = getHap(c1)
                    H2 = getHap(c2)
                    if H1 != H2:  # trans
                        #print('row',rowi,'c2=',c2,'H2=', H2)
                        if H1 == haps[0]:  # x is h1 for trans - otherwise ignore
                            trans1["x"].append(starts[H1][c1] + int(row[2]))
                            trans1["y"].append(starts[H2][c2] + int(row[7]))
                        else:
                            trans1["y"].append(starts[H1][c1] + int(row[2]))
                            trans1["x"].append(starts[H2][c2] + int(row[7]))
                    else:  # cis
                        if H1 == haps[0]:
                            cis1["x"].append(starts[H1][c1] + int(row[2]))
                            cis1["y"].append(starts[H1][c2] + int(row[7]))
                        else:
                            cis2["x"].append(starts[H2][c1] + int(row[2]))
                            cis2["y"].append(starts[H2][c2] + int(row[7]))
        print('ncis1=',len(cis1["x"]))
        if (len(cis1["x"])) > 0:
            hap = haps[0]
            ofn = "%s_cis%s_hseq.gz" % (hap, inFname)
            self.export_mapping(
                holoSeqHeaders[1],
                ofn,
                [hap for i in range(len(h1names))],
                h1names,
                h1starts,
                cis1["x"],
                cis1["y"],
                [],
                args,
            )
        if (len(cis2["x"])) > 0:
            ofn = "%s_cis%s_hseq.gz" % (hap, inFname)
            hap = haps[1]
            self.export_mapping(
                holoSeqHeaders[1],
                ofn,
                [hap for i in range(len(h2names))],
                h2names,
                h2starts,
                cis2["x"],
                cis2["y"],
                [],
                args,
            )
        if (len(trans1["x"])) > 0:
            ofn = "%s_trans_hseq.gz" % (inFname)
            self.export_mapping(
                holoSeqHeaders[1],
                ofn,
                [haps[0] for i in range(len(h1names))]
                + [haps[1] for i in range(len(h2names))],
                h1names + h2names,
                h1starts + h2starts,
                trans1["x"],
                trans1["y"],
                [],
                args,
            )

    def export_mapping(self, hsId, outFname, haps, hnames, hstarts, x, y, anno, args):
        """
        @v1HoloSeq2D for example
        """

        def prepHeader(haps, hnames, hstarts, hsId, args):
            """
            holoSeq output format
            """
            h = ["@%s %s %d" % (haps[i], hnames[i], hstarts[i]) for i in range(len(hnames))]
            metah = [hsId, "@@title %s" % args.title, "@@datasource %s" % "bigwig", "@@datafile %s" % self.inFname, "@@refURI %s" % args.refURI]
            
            return metah + h

        hdr = prepHeader(haps, hnames, hstarts, hsId, args)

        with gzip.open(outFname, mode="wb") as ofn:
            ofn.write(str.encode("\n".join(hdr) + "\n"))
            if len(anno) == len(x):
                for i in range(len(x)):
                    row = str.encode("%d %d %s\n" % (x[i], y[i], anno[i]))
                    ofn.write(row)
            else:
                for i in range(len(x)):
                    row = str.encode("%d %d\n" % (x[i], y[i]))
                    ofn.write(row)


parser = argparse.ArgumentParser(description="", epilog="")
parser.add_argument(
    "--inFile", help="PAF with paired alignments", nargs="+", default=[]
)
parser.add_argument(
    "--title", help="Title for the plot", default="Plot title goes here"
)
parser.add_argument("--contig_sort", help="VGPname, name, length", default="VGPname")
parser.add_argument(
    "--refURI",
    help="URI for the genome reference sequence used for the coordinates for metadata",
    default="Unknown",
)
parser.add_argument(
    "--hap_indicator", help="None, Suffix (H[1,2]) Dashsuffix (_H...)", default="None"
)
parser.add_argument("--version", "-V", action="version", version="0.1")
args = parser.parse_args()
for f in args.inFile:
    ps = Path(f).suffix.lower()
    print("inFile=", f, ps)

    if ps == ".paf":
        p = pafConvert(f, args)
    elif ps in [".bw", ".bigwig"]:
        outf = "%s_bw.hseq.gz" % f
        p = bwConvert(f, outf, args)
    else:
        print(f, "unknown type - cannot process")
