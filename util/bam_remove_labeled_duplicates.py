#!/usr/bin/env python3
# encoding: utf-8

import os, sys, re
import logging
import argparse
import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="mark duplicates in bam",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_bam",
        "-i",
        dest="input_bam",
        required=True,
        type=str,
        help="input bam file, coordinate sorted",
    )

    parser.add_argument(
        "--output_bam",
        "-o",
        dest="output_bam",
        required=True,
        type=str,
        help="output bam file",
    )

    parser.add_argument(
        "--output_dups_bam",
        dest="output_dups_bam",
        required=False,
        type=str,
        help="output bam file",
    )

    parser.add_argument(
        "--debug",
        "-d",
        dest="debug",
        action="store_true",
        default=False,
        help="debug mode",
    )

    args = parser.parse_args()

    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam

    output_dups_bam_filename = args.output_dups_bam

    if args.debug:
        logger.setLevel(logging.DEBUG)

    bamreader = pysam.AlignmentFile(input_bam_filename, "rb")

    if ((not "SO" in bamreader.header.as_dict()["HD"])) or bamreader.header.as_dict()[
        "HD"
    ]["SO"] != "coordinate":
        raise RuntimeError(
            "Error, file: {} must be coordinate sorted".format(input_bam_filename)
        )

    bamwriter = pysam.AlignmentFile(output_bam_filename, "wb", template=bamreader)

    if output_dups_bam_filename:
        dups_bamwriter = pysam.AlignmentFile(
            output_dups_bam_filename, "wb", template=bamreader
        )

    duplicate_counter = 0
    for read in bamreader.fetch():
        if read.is_duplicate:
            duplicate_counter += 1
            if dups_bamwriter:
                dups_bamwriter.write(read)
        else:
            bamwriter.write(read)

    logger.info("Done. Removed {} duplicates".format(duplicate_counter))

    sys.exit(0)


if __name__ == "__main__":
    main()
