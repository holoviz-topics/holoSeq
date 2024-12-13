import argparse
from pathlib import Path

import panel as pn

from holoseq.view import maker

pn.extension()

parser = argparse.ArgumentParser(description="", epilog="")
parser.add_argument("--inFile", help="gzipped hseq coordinates and contigs", nargs="+")
parser.add_argument("--size", help="Display size in pixels. Default is 1000", default=1000)

args = parser.parse_args()
path = Path(args.inFile[0]).resolve()
width = int(args.size)

holoseq_panel, title = maker(path=path, width=width)
pn.panel(holoseq_panel).servable(title=title)
