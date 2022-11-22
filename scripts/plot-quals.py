#!/usr/bin/env python

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input_file",
        required=True,
    )
    parser.add_argument(
        "--output_file",
        required=True,
    )
    args = parser.parse_args()
    quals = [record.qual for record in VariantFile(args.input_file)]
    plt.hist(quals)
    try:
        plt.savefig(args.output_file)
    except Exception as e:
        print(e)

if __name__ == "__main__":
    main()