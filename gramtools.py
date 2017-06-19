import os
import time
import argparse
import logging
import subprocess


def setup_logging():
    log = logging.getLogger('gramtools')
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)
    return log

log = setup_logging()


from py_interface import build, infer





def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--infer", help="",
                        action="store_true")
    parser.add_argument("--build", help="",
                        action="store_true")

    parser.add_argument("--gram-files", help="",
                        type=str)
    parser.add_argument("--fast", help="",
                        type=str)
    parser.add_argument("--vcf", help="",
                        type=str)

    parser.add_argument("--ksize", help="",
                        type=int)

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if args.build:
        build.run(args)

    elif args.infer:
        infer.run(args)
