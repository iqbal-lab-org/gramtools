import os
import time
import argparse
import logging
import subprocess

import build
import infer


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--infer", help="",
                        action="store_true")
    parser.add_argument("--build", help="",
                        action="store_true")

    parser.add_argument("--prg", help="",
                        type=str)
    parser.add_argument("--fast", help="",
                        type=str)
    parser.add_argument("--vcf", help="",
                        type=str)

    parser.add_argument("--ksize", help="",
                        type=int)
    parser.add_argument("--output", help="",
                        type=str)

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if args.build:
        build.run(args)

    elif args.infer:
        infer.run(args)
