# for Mashmap paf, python holoSeq_prepare_paf.py --inFile  hg002_2k99.paf --title "hg002 Mashmap" --hap_indicator None --contig_sort length
# for HiC pairs
# python holoSeq_prepare_paf.py --inFile mUroPar1H1H2.paf --xclenfile mUroPar1H1suffix.len --yclenfile mUroPar1H2suffix.len --contig_sort VGPname --hap_indicator Suffix --title "VGP mUroPar1 HiC data"
# panel serve holoseq_display.py --show --args --inFile mUroPar1H1H2.paf_cisH1_hseq.gz mUroPar1H1H2.paf_cisH2_hseq.gz mUroPar1H1H2.paf_trans_hseq.gz  --size 1000
#
# python holoSeq_prepare_paf.py --inFile mUroPar1_protH1.gff --xclenfile mUroPar1H1suffix.len --contig_sort VGPname --title "mUroPar1 NCBI protein GFF"

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
import itertools
import logging
import math
import random
import re
import os

import pybigtools

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("holoseq_display")

# inFile = "galaxy_inputs/paf/bothmap.paf.tab.tabular"
inFile = "/home/ross/rossgit/holoviews-examples/huge.paf"

holoSeqHeaders = ["@v1HoloSeq1D", "@v1HoloSeq2D"]


def rotatecoords(x, y, radians=0.7853981633974483, origin=(0, 0)):
    # https://gist.github.com/LyleScott/d17e9d314fbe6fc29767d8c5c029c362
    # to rotate so the diagonal becomes the x axis
    # this inflates by sqrt(2) so would need to rescale 1D tracks - not sure that's ideal but worth trying.
    # xcis1r, ycis1r = rotatecoords(xcis1, ycis1, radians=0.7853981633974483, origin=(max(xcis1),max(ycis1)))
    # pafxycis1 = pd.DataFrame(np.vstack([xcis1r,ycis1r]).T, columns = ['x', 'y'])
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
    haps = []
    with open(lenFile) as lf:
        for i, row in enumerate(lf):
            row = [x.strip() for x in row.strip().split()]
            if len(row) > 1:
                c, clen = row[:2]
                if seen.get(c, None):
                    log.debug("Contig %s seen again at row %d of %s" % (c, i, lenFile))
                else:
                    seen[c] = c
                    contigs.append((c, int(clen)))
                h = getHap(c)
                if not h in haps:
                    haps.append(h)
    return contigs, haps


def VGPsortfunc(s1, s2):
    """
    # big fugly hack to sort super contigs before anything else
    # then by contig number or if they are the same offset
    # ('SUPER_2H1', 226668729) , ('SUPER_1H1', 284260672), ('SUPER13_unloc_5H1',..), (Scaffold_aa16H2, ...)
    # or chr10H1 etc
    # always work with uppercase
    # ^([^_]+)_(\S+)H[12]{1}$ gives group1 SUPER or Scaffold and group2 X or 22
    #     ^([^_]+)_([^_]+)_unloc_(\S+)H[12]{1}$ gives super aa x for super_aa_unloc_XH2
    """

    ss1 = re.compile(
        r"^([^_]+)_(\S+)H[12]{1}$"
    )  # should match VGP non-unloc - super/chr/scaffold. Always uppercase at entry.
    ss2 = re.compile(r"^([^_]+)_([^_]+)_UNLOC_(\S+)H[12]{1}$")  # for the unloc

    def intOrd(n):
        "x, y, z, w typically sex chromosomes - there be some weird beasts"
        if n:
            n = n.replace("_", "").replace("chr", "")  # in case chr_aa or something
            if n.isdigit():
                return int(n)
            else:
                return ord(n[0])
        else:
            return None

    def matchme(s):
        found = ss2.search(s)
        if found:
            c1, n1, n2 = found.groups()[:3]
        else:
            found = ss1.search(s)
            if found:
                c1, n1 = found.groups()[:2]
                n2 = None
            else:
                c1 = n1 = n2 = None
        n1 = intOrd(n1)
        n2 = intOrd(n2)
        return c1, n1, n2

    if s1[0] == s2[0]:  # simplest case - same contig, sort on offset
        return s1[1] - s2[1]  # neg if left sorts before
    else:  # deal with chr - assume no dashes?
        if "chr" in s1[0].lower():
            if "chr" in s2[0].lower():  # punt for chr22H2 or so
                n1 = intOrd(s1[0].lower())
                n2 = intOrd(s2[0].lower())
                return n1 - n2
            else:
                return -1  # put chr first
        elif "chr" in s2[0].lower():
            return 1
    u1 = "unloc" in s1[0].lower()
    u2 = "unloc" in s2[0].lower()
    if u1 and not u2:
        return 1  # u1 unloc goes after
    elif u2 and not u1:
        return -1  # u1 goes before unloc
    isSuper1 = (not u1) and (("super" in s1[0].lower()))
    isSuper2 = (not u2) and (("super" in s2[0].lower()))
    isScaff1 = (not u1) and (("scaffold" in s1[0].lower()))
    isScaff2 = (not u2) and (("scaffold" in s2[0].lower()))
    if isSuper1 and not isSuper2:
        return -1
    elif isSuper2 and not isSuper1:
        return 1
    # Must parse

    c1, naa, nab = matchme(s1[0].upper())
    c2, nba, nbb = matchme(s2[0].upper())
    if not c1 or not c2:
        log.debug("no match for ss1 %s and/or ss2 %s" % (ss1, ss2))
        return 0
    else:
        if isSuper1 or (
            isScaff1 and isScaff2
        ):  # if a super must both be supers or if both are scaffolds
            return naa - nba
        else:  # must both be unlocs
            if naa == nba:
                return nab - nbb
            else:
                return naa - nba


def Lengthsortfunc(s1, s2):
    """ """
    return s1[1] - s2[1]  # neg if left sorts before


def contsort(contigs, args):
    # sort and return offsets to starts of each contig
    # hstarts = list(itertools.accumulate(hlens))
    if args.contig_sort.lower() == "vgpname":
        contigs.sort(key=cmp_to_key(VGPsortfunc))
    elif args.contig_sort.lower() == "name":
        contigs.sort()
    elif args.contig_sort.lower() == "length":
        contigs.sort(key=cmp_to_key(Lengthsortfunc), reverse=True)
    clens = [x[1] for x in contigs]
    cnames = [x[0] for x in contigs]
    cstarts = list(itertools.accumulate(clens))
    cstarts.insert(0, 0)  # first one starts at 0
    scont = OrderedDict(zip(cnames, cstarts))
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
        segs = {}
        comment = "#"
        self.hsId = "@v1HoloSeq1D"
        self.inFname = gff
        print("contigs=", str(contigs)[:1000])
        with open(gff) as g:
            for i, row in enumerate(g):
                if not row.startswith(comment):
                    (id, name, kind, startp, endp, score, strand, phase, text) = [
                        x.strip() for x in row.split()[:9]
                    ]
                    if not segs.get(id, None):
                        segs[id] = []
                    if kind.lower() in ["cds", "mrna"]:
                        anno = text.split(";")
                        tanno = [
                            x.strip()[7:]
                            for x in anno
                            if x.lower().startswith("target=")
                        ]
                        target = tanno[0]
                    startp = int(startp)
                    endp = int(endp)
                    offset = contigs.get(id, -1)
                    if args.addH1 and offset < 0:
                        offset = contigs.get(id + "H1", -1)
                        id = id + "H1"
                    if offset < 0:
                        print(
                            "Ignored gff3 id %s missing from supplied xcontigs, in row %d %s of %s with addH1=%s"
                            % (id, i, row, gff, args.addH1)
                        )
                    else:
                        if kind.lower() == "mrna":
                            if target:
                                if mrnaseen.get(target, None):
                                    print(
                                        "Seeing mrna target %s again at row %d"
                                        % (target, i)
                                    )
                                else:
                                    mrnaseen[target] = target
                                segs[id].append(
                                    (
                                        startp + offset,
                                        endp + offset,
                                        strand,
                                        score,
                                        target,
                                        "mrna",
                                    )
                                )
                            else:
                                print("no target found in %s at row %d" % (text, i))
                        elif kind.lower() == "stop_codon":
                            segs[id].append((startp + offset, target, "stopc"))
                        elif kind.lower() == "cds":
                            segs[id].append(
                                (
                                    startp + offset,
                                    endp + offset,
                                    strand,
                                    score,
                                    target,
                                    "cds",
                                )
                            )

        self.export_mapping(outFname, contigs, segs, args)

    def export_mapping(self, outFname, contigs, segs, args):
        """
        for GFF
        @v1HoloSeq2D for example
        A default  Y value of 100 is set for each mRNA's extent, but often there are dozens of different named sequences in the databases that will overlap.
        Assuming the GFF is sorted, overlapping mRNA is identified as ending or starting in the previous extent and the y value is decremented to avoid overlap
        Y is reset if no overlap
        """

        def prepHeader(contigs, args):
            """
            holoSeq output format
            """
            h = ["@%s %s %d" % (getHap(k), k, contigs[k]) for k in contigs.keys()]
            metah = [
                self.hsId,
                "@@GFF 1",
                "@@title %s" % args.title,
                "@@datasource GFF",
                "@@datafile %s" % self.inFname,
                "@@refURI %s" % args.refURI,
                "@@xclenfile %s" % args.xclenfile,
            ]

            return metah + h

        def ranges_overlap(x1, x2, y1, y2):
            # buried in https://stackoverflow.com/questions/6821156/how-to-find-range-overlap-in-python
            if x1 == x2 or y1 == y2:
                return False
            return x1 <= y2 and y1 <= x2

        hdr = prepHeader(contigs, args)

        with gzip.open(outFname, mode="wb") as ofn:
            ofn.write(str.encode("\n".join(hdr) + "\n"))
            y = 100
            for con in contigs.keys():
                lastseg = ("", 0, 0)
                subs = segs.get(con, [])
                if len(subs) > 0:
                    subs.sort(key=lambda x: x[0])
                    for i, m in enumerate(subs):
                        kind = m[-1]
                        if kind == "mrna":
                            (startp, endp, strand, score, targ, _) = m
                            if ranges_overlap(lastseg[1], lastseg[2], startp, endp):
                                y -= 1
                            else:
                                y = 100
                                lastseg = (targ, startp, endp)
                            row = str.encode(
                                f"mrna {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                            )
                            ofn.write(row)
                        elif kind == "cds":
                            (startp, endp, strand, score, targ, _) = m
                            row = str.encode(
                                f"cds {targ} {con} {startp} {endp} {y} {y} {strand} {score}\n"
                            )
                            ofn.write(row)
                        elif kind == "stopc":
                            (startp, targ, _) = m
                            row = str.encode(f"stopc {targ} {con} {startp}\n")
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
            h = ["@%s %s %d" % (getHap(k), k, contigs[k]) for k in contigs.keys()]
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

    def __init__(self, inFname, args, xcontigs, ycontigs, haps):
        self.inFname = inFname
        # have the axes set up so prepare the three plot x/y vectors
        # for a second pass to calculate all the coordinates.
        # adding tooltips just does not scale so abandoned - see the tooltip old version
        rowi = 0
        ncis1 = ncis2 = ntrans = 0
        # print('ycon %s' % (ycontigs))
        hsId = holoSeqHeaders[1]
        self.inFname = inFname
        self.prepPafGZ(hsId, haps, xcontigs, ycontigs, args)
        with open(inFname) as f:
            for rowi, rows in enumerate(f):
                row = rows.strip().split()
                if len(row) > 7:
                    c1 = row[0]
                    c2 = row[5]
                    n1 = int(row[2])
                    n2 = int(row[7])
                    H1 = getHap(c1)
                    H2 = getHap(c2)
                    if H1 != H2:  # trans
                        if H1 == haps[0]:  # x is h1 for trans - otherwise ignore
                            xstart = xcontigs[c1]
                            ystart = ycontigs[c2]
                            row = str.encode("%d %d\n" % (xstart + n1, ystart + n2))
                            self.transf.write(row)
                            ntrans += 1
                        else:
                            xstart = xcontigs[c2]
                            ystart = ycontigs[c1]
                            row = str.encode("%d %d\n" % (xstart + n2, ystart + n1))
                            self.transf.write(row)
                            ntrans += 1
                    else:  # cis
                        if H1 == haps[0]:
                            xstart = xcontigs[c1]
                            ystart = xcontigs[c2]
                            row = str.encode("%d %d\n" % (xstart + n1, ystart + n2))
                            self.cis1f.write(row)
                            ncis1 += 1
                        else:
                            xstart = ycontigs[c1]
                            ystart = ycontigs[c2]
                            row = str.encode("%d %d\n" % (xstart + n1, ystart + n2))
                            self.cis2f.write(row)
                            ncis2 += 1
        self.cis1f.close()
        self.cis2f.close()
        self.transf.close()
        print("ncis1=%d, ncis2=%d, ntrans=%d" % (ncis1, ncis2, ntrans))

    def prepPafGZ(self, hsId, haps, xcontigs, ycontigs, args):
        """
        @v1HoloSeq2D for example
        """

        def prepHeader(haps, hsId, xcontigs, ycontigs, args, outf):
            """
            holoSeq output format - prepare gzip output channels
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

            outs = "\n".join(metah + h) + "\n"
            outf.write(str.encode(outs))

        self.cis1f = gzip.open("%s_cis%s_hseq.gz" % (self.inFname, haps[0]), mode="wb")
        prepHeader(haps[0], hsId, xcontigs, ycontigs, args, self.cis1f)
        self.cis2f = gzip.open("%s_cis%s_hseq.gz" % (self.inFname, haps[1]), mode="wb")
        prepHeader(haps[1], hsId, xcontigs, ycontigs, args, self.cis2f)
        self.transf = gzip.open("%s_trans_hseq.gz" % (self.inFname), mode="wb")
        prepHeader(haps, hsId, xcontigs, ycontigs, args, self.transf)


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
    haps = []
    yhaps = []
    xcontigs, xhaps = getContigs(args.xclenfile)
    sxcontigs = contsort(xcontigs, args)
    if args.yclenfile:
        ycontigs, yhaps = getContigs(args.yclenfile)
        sycontigs = contsort(ycontigs, args)
    else:
        sycontigs = sxcontigs
    for h in xhaps + yhaps:
        if not h in haps:
            haps.append(h)
    for f in args.inFile:
        ps = Path(f).suffix.lower()
        print("inFile=", f, ps)

        if ps == ".paf":
            p = pafConvert(f, args, sxcontigs, sycontigs, haps)
        elif ps in [".bw", ".bigwig"]:
            outf = "%s.hseq.gz" % f
            p = bwConvert(f, outf, args, sxcontigs)
        elif ps in [".gff3", ".gff"]:
            outf = "%s.hseq.gz" % f
            p = gffConvert(f, outf, sxcontigs, args)
        else:
            print(f, "unknown type - cannot process")
