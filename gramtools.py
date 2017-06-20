import os
import time
import argparse
import logging
import subprocess

from py_interface import build, quasimap


def setup_logging(level):
    log = logging.getLogger('gramtools')
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(level)
    return log


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--debug", help="",
                        action="store_true")

    parser.add_argument("--quasimap", help="",
                        action="store_true")
    parser.add_argument("--profile", help="",
                        action="store_true")
    parser.add_argument("--build", help="",
                        action="store_true")

    parser.add_argument("--gram-files", help="",
                        type=str)
    parser.add_argument("--reference", help="",
                        type=str)
    parser.add_argument("--vcf", help="",
                        type=str)

    parser.add_argument("--ksize", help="",
                        type=int)

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    setup_logging(level)

    if args.build:
        build.run(args)

    elif args.quasimap:
        quasimap.run(args)
