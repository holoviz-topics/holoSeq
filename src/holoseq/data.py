import array
import gzip
from pathlib import Path

from holoseq.config import VALID_HSEQ_FORMATS


def is_valid_header(header: str) -> tuple[bool, int, str]:
    tokens = [token.strip() for token in header.split()]

    check = False
    hseq_format = tokens[0]
    if hseq_format not in VALID_HSEQ_FORMATS:
        msg = f"Not a valid holoSeq file. Header must start with one of: {VALID_HSEQ_FORMATS}."
        raise IOError(msg)
    else:
        check = True

    num_dimensions = VALID_HSEQ_FORMATS.index(hseq_format) + 1
    plot_type = "bar"
    if num_dimensions == 1 and len(tokens) > 1:
        plot_type = tokens[1]

    return check, num_dimensions, plot_type


def load(path: Path):
    haploids = {}
    x_coords = array.array("l")
    y_coords = array.array("l")
    annotations = []
    metadata = {}
    gff_data = []
    title = "Plot"

    is_gff = False
    num_dimensions = -1
    plot_type = "bar"

    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            # Check if we have a valid holoSeq file using the first line as a header, and raise an
            # error if it is not.
            if i == 0:
                _, num_dimensions, plot_type = is_valid_header(line)
                continue

            if line[0] == "@":
                row = line[1:]
                if row[0] == "@":
                    tokens = [token.strip() for token in row[1:].split()]
                    metadata[tokens[0]] = tokens[1:]
                    if tokens[0] == "GFF":
                        is_gff = True

                else:
                    tokens = [token.strip() for token in row.split()]

                    if len(tokens) >= 3:
                        haploid_name, contig_name, position = tokens[:3]
                        if not haploids.get(haploid_name, None):
                            haploids[haploid_name] = {
                                "contig_names": [],
                                "positions": array.array("l"),
                            }

                        if num_dimensions == 2:
                            haploids[haploid_name]["contig_names"].append(contig_name)
                            haploids[haploid_name]["positions"].append(int(position))
                    else:
                        msg = (
                            "NOT A VALID holoSeq FILE.\n"
                            f"Line {i} of {str(path)} does not include a `reference name` "
                            "`contig name`, and `contig length`."
                        )
                        raise IOError(msg)
            else:
                tokens = [token.strip() for token in line.split()]
                num_tokens = len(tokens)

                if num_dimensions == 2:
                    if num_tokens < 2:
                        msg = (
                            "NOT A VALID holoSeq FILE.\n"
                            f"Line {i} of {str(path)} needs at least two valid integer "
                            "coordinates to be a valid 2D holoSeq file."
                        )
                        raise IOError(msg)

                    if num_tokens >= 2:
                        coords_check = tokens[0].isdigit() and tokens[1].isdigit()
                        if coords_check:
                            x_coords.append(int(tokens[0]))
                            y_coords.append(int(tokens[1]))
                            if num_tokens > 2:
                                annotations.append(tokens[2:])
                        else:
                            msg = (
                                "NOT A VALID holoSeq FILE.\n"
                                f"Line {i} of {str(path)} needs at least two valid integer "
                                "coordinates to be a valid 2D holoSeq file."
                            )
                            raise IOError(msg)
                else:
                    if is_gff:
                        gff_data.append(tokens)

                    else:
                        if tokens[0].isdigit():
                            x_coords.append(int(tokens[0]))
                            if num_tokens > 1:
                                y_coords.append(int(tokens[1]))
                            if num_tokens > 2:
                                annotations.append(tokens[2:])
                        else:
                            msg = (
                                "NOT A VALID holoSeq FILE.\n"
                                f"Line {i} of {str(path)} needs at least one valid integer to be "
                                "a valid 1D holoSeq file."
                            )
                            raise IOError("Not a valid holoSeq file.")

    return haploids, x_coords, y_coords, annotations, title, plot_type
