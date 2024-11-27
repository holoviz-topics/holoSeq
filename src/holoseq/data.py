import array
import gzip
from pathlib import Path

from holo_seq.config import VALID_HSEQ_FORMATS


def load_2d(path: Path):
    """Reduces memory footprint by ~50%."""
    haploids = {}
    x_coords = array.array("l")
    y_coords = array.array("l")
    annotations = []
    num_dimensions = -1
    title = "Plot"

    sentinel = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            # Check if we have a valid holoSeq file using the first line as a header, and raise an
            # error if it is not.
            if sentinel == 0:
                hseq_format = line[:12]
                if hseq_format not in VALID_HSEQ_FORMATS:
                    raise IOError("Not a valid holoSeq file.")

                num_dimensions = VALID_HSEQ_FORMATS.index(hseq_format) + 1
                sentinel += 1

            else:
                if line[0] == "@":
                    row = line[1:]

                    # If the line is a title, extract it.
                    if row[:6] == "title ":
                        title = row[6:]

                    else:
                        tokens = [token.strip() for token in row.split()]
                        if len(tokens) == 3:
                            haploid_name, contig_name, position = tokens

                            if not haploids.get(haploid_name, None):
                                haploids[haploid_name] = {
                                    "contig_names": [],
                                    "positions": array.array("l"),
                                }
                            haploids[haploid_name]["contig_names"].append(contig_name)
                            haploids[haploid_name]["positions"].append(int(position))
                        else:
                            raise IOError("Not a valid holoSeq file.")

                else:
                    tokens = [token.strip() for token in line.split()]
                    num_tokens = len(tokens)
                    if num_dimensions == 2:
                        if num_tokens < 2:
                            raise IOError("Not a valid holoSeq file.")
                        else:
                            if tokens[0].isdigit() and tokens[1].isdigit():
                                x_coords.append(int(tokens[0]))
                                y_coords.append(int(tokens[1]))
                                if num_tokens > 2:
                                    annotations.append(tokens[2:])
                            else:
                                raise IOError("Not a valid holoSeq file.")
                    else:
                        if tokens[0].isdigit():
                            x_coords.append(int(tokens[0]))
                            if num_tokens > 1:
                                annotations.append(tokens[1:])
                        else:
                            raise IOError("Not a valid holoSeq file.")

    return haploids, x_coords, y_coords, annotations, title
