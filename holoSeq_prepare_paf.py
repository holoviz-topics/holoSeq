# for Mashmap paf, python holoSeq_prepare_paf.py --inFile  hg002_2k99.paf --title "hg002 Mashmap" --hap_indicator None --contig_sort length
# for HiC pairs
# python holoSeq_prepare_paf.py --inFile mUroPar1H1H2.paf --xclenfile mUroPar1H1suffix.len --yclenfile mUroPar1H2suffix.len --contig_sort VGPname --hap_indicator Suffix --title "VGP mUroPar1 HiC data"
#
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
import random
import re
import os

import pybigtools


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
    function to return suffix H1 from "chrH1" - adjust to suit the specific conventions used by
    whoever assembles your genomes.
     help="None, Suffix (H[1,2]) Dashsuffix (_H...)"
    """
    hi = args.hap_indicator
    if hi == "None":
        return "H1"
    elif hi == "Suffix":
        return contig[-2:]
    elif hi == "Dashsuffix":
        return contig.split("_")[-1]


def getContigs(lenFile):
    # samtools faidx will make one of these from a genome fasta
    # whitespace delimited contig names and lengths.
    contigs = []
    seen = {}
    with open(lenFile) as lf:
        for i, row in enumerate(lf):
            c, clen = row.split()[:2]
            if seen.get(c, None):
                print("Contig %s seen again at row %d of %s" % (c, i, lenFile))
            else:
                seen[c] = c
                contigs.append((c, int(clen)))
    return contigs

def VGPsortfunc(s1, s2):
    """
    # fugly hack to sort super contigs before anything else
    # then by contig number or if they are the same offset
    # ('SUPER_2H1', 226668729) , ('SUPER_1H1', 284260672), ('SUPER13_unloc_5H1',..), (Scaffold_1116H2, ...)
    # or chr10_H1 etc
    # always work with uppercase
    # ^([^_]+)_(\S+)H[12]{1}$ gives group1 SUPER or Scaffold and group2 X or 22
    #     ^([^_]+)_([^_]+)_unloc_(\S+)H[12]{1}$ gives super 11 x for super_11_unloc_XH2
    """

    ss1 = re.compile(r'^([^_]+)_(\S+)H[12]{1}$')
    ss2 = re.compile(r'^([^_]+)_([^_]+)_UNLOC_(\S+)H[12]{1}$')

    def matchme(s):
        c1 = n1 = n2 = None
        found = ss2.search(s)
        if found:
            c1, n1, n2 = found.groups()[:3]
        else:
            found = ss1.search(s)
            if found:
                c1, n1 = found.groups()[:2]
        if n1:
            if n1.isdigit():
                n1 = int(n1)
            else:
                n1 = ord(n1[0])
        if n2:
            if n2.isdigit():
                n2 = int(n2)
            else:
                n2 = ord(n2[0])
        #print('c1, n1, n2 =', c1, n1, n2)
        return c1, n1, n2

    if s1[0] == s2[0]:  # simplest case - same contig, sort on offset
        return s1[1] - s2[1]  # neg if left sorts before
    u1 = "unloc" in s1[0].lower()
    u2 = "unloc" in s2[0].lower()
    if u1 and not u2:
        return 1 # unloc goes after
    elif u2 and not u1:
        return -1 # unloc goes after
    isSuper1 = (not u1) and (("super" in s1[0].lower()) or ("chr" in s1[0].lower()))
    isSuper2 = (not u2) and (("super" in s2[0].lower()) or ("chr" in s2[0].lower()))
    isScaff1 = (not u1) and (("scaffold" in s1[0].lower()))
    isScaff2 = (not u2) and (("scaffold" in s2[0].lower()))
    if isSuper1 and not isSuper2:
        return -1
    elif isSuper2 and not isSuper1:
        return 1
    # Must parse
    
    c1, n11, n12 = matchme(s1[0].upper())
    c2, n21, n22 = matchme(s2[0].upper())
    if not c1 or not c2:
        print('no match for ss1 and/or ss2 =', ss1, ss2)
        return 0
    else:
        if isSuper1 or isScaff1: # must both be supers or scaffolds
            return n11 - n21
        else: # must both be unlocs
            #print('both unlocs s1, s2, n11, n12, n21, n22', s1, s2, n11, n12, n21, n22)
            if n11 == n21:
                return n12 - n22
            else:
                return n11 - n21


def VGPsortfunc0(s1, s2):
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
        if "UNLOC" in contig.upper():
            unloc = True
            contign = conta.split("_")[0]
        else:
            unloc = False
            contign = contb
        if conta.startswith("CHR"):
            contign = conta.replace("CHR", "").replace("_", "")
        if contign.isdigit():
            nval = int(contign)
        elif contign[:-2].isdigit():  # suffix H2 is VGP standard
            nval = int(contign[:-2])
        else:
            nval = ord(contign[0])  # ensure X and Y sort to the end
        is_super = not unloc
        sorts[i]["is_super"] = is_super
        sorts[i]["n"] = nval
        sorts[i]["conta"] = conta
        sorts[i]["unloc"] = unloc
    if sorts[0]["is_super"] and not sorts[1]["is_super"]:
        return -1
    elif sorts[1]["is_super"] and not sorts[0]["is_super"]:
        return 1
    elif sorts[0]["is_super"] and sorts[1]["is_super"]:
        return sorts[0]["n"] - sorts[1]["n"]
    elif sorts[0]["conta"] == sorts[1]["conta"]:
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


def contsort(contigs, args):
    # sort and return offsets to starts of each contig
    if args.contig_sort.lower() == "vgpname":
        contigs.sort(key=cmp_to_key(VGPsortfunc))
    elif args.contig_sort.lower() == "name":
        contigs.sort()
    elif args.contig_sort.lower() == "length":
        contigs.sort(key=cmp_to_key(Lengthsortfunc), reverse=True)
    scont = OrderedDict()
    cum = 0
    for i, (cont, clen) in enumerate(contigs):
        scont[cont] = cum
        cum += int(clen)
    #print('scont=', scont)
    return scont


class gffConvert:
    """
            Only care about mRNA cds and stop codons initally. Turn into segments. Filter so only data in contigs is retained from input.
    SUPER_1 miniprot        mRNA    139006290       139072696       22660   -       .       ID=MP000006;Rank=1;Identity=0.9979;Positive=0.9984;Target=XP_026244093.1 1 4350
        fix positions as we go with a lookup contig -> cumulated offset
        can pop open https://www.ncbi.nlm.nih.gov/protein/XP_026244093.1
    """

    def __init__(self, gff, outFname, contigs, args):
        mrnaseen = {}
        stopcodons = {}
        mrna = []
        cds = {}
        comment = "#"
        self.hsId = "@v1HoloSeq2D"
        self.inFname = gff
        with open(gff) as g:
            for i, row in enumerate(g):
                if not row.startswith(comment):
                    rs = row.split()
                    (id, name, kind, startp, endp, score, strand, phase, text) = rs[:9]
                    startp = int(startp)
                    endp = int(endp)
                    offset = contigs.get(id, None)
                    if not offset and args.addH1:
                        offset = contigs.get(id + "H1", None)
                        id = id + "H1"
                    if not offset:
                        print(
                            "Ignored gff3 id %s missing from supplied xcontigs, in row %d of %s with addH1=%s"
                            % (id, i, gff, args.addH1)
                        )
                    else:
                        if kind.lower() == "mrna":
                            anno = text.split(";")
                            tanno = [
                                x.strip()[7:]
                                for x in anno
                                if x.lower().startswith("target=")
                            ]
                            target = tanno[0]
                            if target:
                                if mrnaseen.get(target, None):
                                    print(
                                        "Seeing mrna target %s again at row %d"
                                        % (target, i)
                                    )
                                else:
                                    mrnaseen[target] = target
                                mrna.append(
                                    (
                                        id,
                                        target,
                                        startp + offset,
                                        endp + offset,
                                        strand,
                                        score,
                                    )
                                )
                            else:
                                print("no target found in %s at row %d" % (text, i))
                        elif kind.lower() == "stop_codon" and target:
                            stopcodons[target] = startp + offset
                        elif kind.lower() == "cds":
                            if cds.get(target, None):
                                cds[target].append(
                                    (id, startp + offset, endp + offset, strand, score)
                                )
                            else:
                                cds[target] = [
                                    (id, startp + offset, endp + offset, strand, score)
                                ]

        self.export_mapping(outFname, contigs, mrna, stopcodons, cds, args)

    def export_mapping(self, outFname, contigs, mrna, stopcodons, cds, args):
        """
        for GFF
        @v1HoloSeq2D for example
        """

        def prepHeader(contigs, args):
            """
            holoSeq output format
            """
            h = ["@%s %d" % (k, contigs[k]) for k in contigs.keys()]
            metah = [
                self.hsId,
                "@@GFF",
                "@@title %s" % args.title,
                "@@datasource %s" % "gff",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % args.xclenfile,
            ]

            return metah + h

        hdr = prepHeader(contigs, args)

        with gzip.open(outFname, mode="wb") as ofn:
            ofn.write(str.encode("\n".join(hdr) + "\n"))
            lastseg = ("", 0, 0)
            y = 100
            for i, m in enumerate(mrna):
                (contig, targ, startp, endp, strand, score) = m
                if startp >= lastseg[1] and startp <= lastseg[2]:
                    y -= 2
                else:
                    y = 100
                    lastseg = (targ, startp, endp)
                stopc = stopcodons.get(targ, -1)
                row = str.encode(
                    f"mrna {targ} {contig} {startp} {endp} {y} {y} {strand} {score} {stopc}\n"
                )
                ofn.write(row)
                for cd in cds[targ]:
                    row = str.encode(
                        f"cds {targ} {cd[0]} {cd[1]} {cd[2]} {y} {y} {cd[3]} {cd[4]}\n"
                    )
                    ofn.write(row)


class bwConvert:
    """
    points = [(50*i, 100+random.random()) for i in range(10000)]
    hv.Curve(points).opts(interpolation='steps-post').opts(width=1000)
    bigwig to barchart - will autoscale to line?
    """

    def __init__(self, inFname, outFname, args, contigs):
        self.inFname = inFname
        self.hsId = "@v1HoloSeq2D"
        fakepath = "in.bw"
        if os.path.isfile(fakepath):
            os.remove(fakepath)
        p = Path(fakepath)
        p.symlink_to(inFname)  # required by pybigtools (!)
        bwf = pybigtools.open(fakepath)
        bchrlist = bwf.chroms()
        bwchrs = list(bchrlist.keys())
        data = {}
        for i, bchr in enumerate(bwchrs):
            cchr = bchr
            if (not contigs.get(bchr, None)) and args.addH1:
                cchr = cchr + "H1"
            if contigs.get(cchr, None):
                cstart = contigs[cchr]
                data[cchr] = {}
                bw = bwf.records(bchr)
                # Return the records of a given range on a chromosome. The result is an iterator of tuples. For BigWigs, these tuples are in the format (start: int, end: int, value: float).
                data[cchr]["xstart"] = [x[0] + cstart for x in bw]
                bw = bwf.records(bchr)
                data[cchr]["xend"] = [x[1] + cstart for x in bw]
                bw = bwf.records(bchr)
                data[cchr]["xval"] = [x[2] for x in bw]
            else:
                print(
                    "Bigwig contig %s not found in supplied X axis lengths file" % cchr
                )
        self.export_mapping(outFname, contigs, data, args)

    def export_mapping(self, outFname, contigs, data, args):
        """
        for bigwig
        @v1HoloSeq2D for example
        """

        def prepHeader(contigs, args):
            """
            holoSeq output format
            """
            h = ["@%s %d" % (k, contigs[k]) for k in contigs.keys()]
            metah = [
                self.hsId,
                "@@bigwig" "@@title %s" % args.title,
                "@@datasource %s" % "bigwig",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % args.xclenfile,
            ]

            return metah + h

        hdr = prepHeader(contigs, args)

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
    Assumes pairs of points representing HiC contact pairs, or Mashmap sequence similarity hits. HiC data typically comes from pairs of haplotypes and is used to help assemble all
    the contigs into chromosomes

    These coordinate pairs are of 3 types - both ends on one or the other reference sequence, or one end on each for HiC pairs.
    Let's call these cis if the same reference and trans if between the two references. The ones in cis are likely to be consistent-ish between the two haplotypes, but the trans ones turn
    out to be very interesting...

    For HiC data, there should really be little difference - whether a haplotype contacts another haplotype or itself depends on the 3D folding and since the haplotypes are wound together into
    a helix, then the whole total contig length - about 2 meters for mammals are folded up into a 10μ^3 ball.

    python holoSeq_prepare_paf.py --inFile mUroPar1H1H2.paf --xclenfile mUroPar1H1suffix.len --yclenfile mUroPar1H2suffix.len --contig_sort VGPname --hap_indicator Suffix
    """

    def __init__(self, inFname, args, xcontigs, ycontigs):
        self.inFname = inFname
        # have the axes set up so prepare the three plot x/y vectors
        # for a second pass to calculate all the coordinates.
        # adding tooltips just does not scale so abando - see the tooltip old version
        cis1 = {"x": [], "y": []}
        cis2 = {"x": [], "y": []}
        trans1 = {"x": [], "y": []}
        haps = []
        rowi = 0
        # print('ycon %s' % (ycontigs))
        with open(inFname) as f:
            for rowi, rows in enumerate(f):
                row = rows.strip().split()
                if len(row) > 7:
                    c1 = row[0]
                    c2 = row[5]
                    n1 = int(row[2])
                    n2 = int(row[7])
                    H1 = getHap(c1)
                    if rowi == 0:
                        haps.append(H1)
                    H2 = getHap(c2)
                    if H2 not in haps:
                        haps.append(H2)
                    if H1 != H2:  # trans
                        if H1 == haps[0]:  # x is h1 for trans - otherwise ignore
                            xstart = xcontigs[c1]
                            ystart = ycontigs[c2]
                            trans1["x"].append(xstart + n1)
                            trans1["y"].append(ystart + n2)
                        else:
                            xstart = xcontigs[c2]
                            ystart = ycontigs[c1]
                            trans1["y"].append(ystart + n1)
                            trans1["x"].append(xstart + n2)
                    else:  # cis
                        if H1 == haps[0]:
                            xstart = xcontigs[c1]
                            ystart = xcontigs[c2]
                            cis1["x"].append(xstart + n1)
                            cis1["y"].append(ystart + n2)
                        else:
                            xstart = ycontigs[c1]
                            ystart = ycontigs[c2]
                            cis2["x"].append(xstart + n1)
                            cis2["y"].append(ystart + n2)
        print("ncis1=", len(cis1["x"]))
        print("ncis2=", len(cis2["x"]))
        if (len(cis1["x"])) > 0:
            hap = haps[0]
            ofn = "%s_cis%s_hseq.gz" % (inFname, hap)
            self.export_mapping(
                holoSeqHeaders[1],
                ofn,
                haps,
                xcontigs,
                ycontigs,
                cis1["x"],
                cis1["y"],
                [],
                args,
            )
        if (len(cis2["x"])) > 0:
            hap = haps[1]
            ofn = "%s_cis%s_hseq.gz" % (inFname, hap)
            hap = haps[1]
            self.export_mapping(
                holoSeqHeaders[1],
                ofn,
                haps,
                ycontigs,
                ycontigs,
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
                haps,
                xcontigs,
                ycontigs,
                trans1["x"],
                trans1["y"],
                [],
                args,
            )

    def export_mapping(
        self, hsId, outFname, haps, xcontigs, ycontigs, x, y, anno, args
    ):
        """
        @v1HoloSeq2D for example
        """

        def prepHeader(haps, hsId, xcontigs, ycontigs, args):
            """
            holoSeq output format
            """
            h = ["@%s %s %d" % (getHap(k), k, xcontigs[k]) for k in xcontigs.keys()]
            if len(haps) > 1:
                h += [
                    "@%s %s %d" % (getHap(k), k, ycontigs[k]) for k in ycontigs.keys()
                ]
            metah = [
                hsId,
                "@@heatmap",
                "@@title %s" % args.title,
                "@@datasource %s" % "bigwig",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % args.xclenfile,
                "@@yclenfile %s" % args.yclenfile,
            ]

            return metah + h

        hdr = prepHeader(haps, hsId, xcontigs, ycontigs, args)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument(
        "--inFile",
        help="PAF with paired alignments, bigwig or gff3. Separate multiple with space. Must have matching .paf, .bw/.bigwig, or .gff/.gff3 extension to parse",
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--xclenfile",
        help="X axis contig names and lengths, whitespace delimited",
        required=True,
    )
    parser.add_argument(
        "--addH1",
        help="Bigwig and gff contigs can have H1 added if that matches the supplied xclenfile contig names. Not recommended - best to map against the right fasta",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--yclenfile",
        help="Optional Y axis contig names and lengths, whitespace delimited for different reference sequences",
        required=False,
    )
    parser.add_argument(
        "--title", help="Title for the plot", default="Plot title goes here"
    )
    parser.add_argument(
        "--contig_sort", help="VGPname, name, length, none", default="length"
    )
    parser.add_argument(
        "--refURI",
        help="URI for the genome reference sequence used for the coordinates for metadata",
        default="Unknown",
    )
    parser.add_argument(
        "--hap_indicator",
        help="None, Suffix (H[1,2]) Dashsuffix (_H...)",
        default="None",
    )
    parser.add_argument("--version", "-V", action="version", version="0.1")
    args = parser.parse_args()
    xcontigs = getContigs(args.xclenfile)
    sxcontigs = contsort(xcontigs, args)
    if args.yclenfile:
        ycontigs = getContigs(args.yclenfile)
        sycontigs = contsort(ycontigs, args)
    else:
        sycontigs = sxcontigs
    for f in args.inFile:
        ps = Path(f).suffix.lower()
        print("inFile=", f, ps)

        if ps == ".paf":
            p = pafConvert(f, args, sxcontigs, sycontigs)
        elif ps in [".bw", ".bigwig"]:
            outf = "%s.hseq.gz" % f
            p = bwConvert(f, outf, args, sxcontigs)
        elif ps in [".gff3", ".gff"]:
            outf = "%s.hseq.gz" % f
            p = gffConvert(f, outf, sxcontigs, args)
        else:
            print(f, "unknown type - cannot process")
