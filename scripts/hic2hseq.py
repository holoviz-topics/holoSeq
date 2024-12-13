"""Convert a .hic file into HoloSeq hseq format.

HoloSeq hseq format is described in holoSeq/docs/HoloSeqOverview.md

# Conversion from .hic format to hseq

This script converts a .hic file into the HoloSeq hseq format using the
[straw](https://github.com/aidenlab/straw/) library. The hic format
stores a symmetric matrix, and the converter only provides the lower
triangle of counts to distinguish it from asymmetric matrices.

## Usage:

1. Install the package with the required dependencis:
   ```
   pip install holoSeq[converters]
   ```

2. Download a .hic file. For example, in the series
   [GSE207951](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207951),
   click on "(custom)" in the table at the bottom, check a single .hic
   file, then click the "Download" button. Untar it.

3. Run the conversion script:
   ```
   python3 scripts/hic2hseq.py input.hic output.hseq.gz
   ```

4. Use the output file in the HoloSeq environment to display the matrix:
   ```
   panel serve scripts/holoseq_display.py --show --args --inFile output.hseq.gz --size 1000
   ```

## Example Data:

A pre-converted HTAN file for [GSM6326543 is
available (530MB)](https://pub-867b121072f54b4a9eecdf01cd27246b.r2.dev/GSM6326543_A001C007.hg38.nodups.pairs.hseq.gz).

See hseq conversion output of GSM6326543 at `holoSeq/docs/assets/README_GSM6326543.png`
"""

import gzip
import logging
import os
import argparse
from typing import List

import hicstraw

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)


class GzipOut:
    """Context manager for writing gzipped files."""

    def __init__(self, fn: str):
        self.fn = fn
        self.f = gzip.open(fn, "wb")

    def write(self, x: str) -> None:
        self.f.write(str.encode(x))

    def __enter__(self):
        return self

    def __exit__(self, *_) -> None:
        self.f.close()


def cumulative_sum(x: List[int]) -> List[int]:
    if not x:
        return x
    y = [x[0]]
    for n in x[1:]:
        y.append(n + y[-1])
    return y


def convert_hic_to_hseq(hicfn: str, ostream: GzipOut, max_chrom: int, title: str) -> None:
    """Convert a .hic file into HoloSeq hseq format.

    Parameters:
        hicfn (str): Path to the input .hic file.
        ostream (GzipOut): GzipOut instance for writing the output.
        max_chrom (int): Maximum number of chromosomes to include.
        title (str): Title for the output matrix.
    """
    if not os.path.isfile(hicfn):
        raise FileNotFoundError(f"Input file '{hicfn}' does not exist.")

    hic = hicstraw.HiCFile(hicfn)
    chroms = hic.getChromosomes()

    if max_chrom < 1 or max_chrom > len(chroms):
        max_chrom = len(chroms)

    resolution = min(hic.getResolutions())
    offsets = cumulative_sum([c.length for c in chroms])
    cnames = [c.name for c in chroms]

    ostream.write(
        f"@v1HoloSeq2D\n@title {title}\n"
        + "".join(f"@H1 {chrom} {offset}\n" for chrom, offset in zip(cnames, offsets))
    )

    xbase = 0
    for i, cx in enumerate(chroms[:max_chrom]):
        ybase = 0
        for cy in chroms[: i + 1]:
            logger.info(f"Processing {cx.name} vs. {cy.name}")

            try:
                for result in hicstraw.straw(
                    "observed", "NONE", hicfn, cx.name, cy.name, "BP", resolution
                ):
                    x = xbase + result.binY
                    y = ybase + result.binX
                    ostream.write(f"{x} {y}\n" * int(result.counts))
            except Exception as e:
                logger.error(f"Error processing {cx.name} vs. {cy.name}: {e}")

            ybase += cy.length
        xbase += cx.length


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("hicfile", help="Input .hic format filename or URL")
    ap.add_argument("hseqgz", help="Output filename (will be gzipped)")
    ap.add_argument("--title", default="", help="Matrix title (default: input hic path)")
    ap.add_argument(
        "--max-chrom",
        type=int,
        default=0,
        help="Maximum number of chromosomes to convert (default: all)",
    )

    argv = ap.parse_args()

    argv.title = argv.title or argv.hicfile

    try:
        with GzipOut(argv.hseqgz) as ostream:
            convert_hic_to_hseq(argv.hicfile, ostream, argv.max_chrom, argv.title)
        logger.info("Conversion completed successfully.")
    except Exception as e:
        logger.error(f"Failed to convert .hic to HoloSeq format: {e}")
