"""Align tracks."""


import sys
import logging

from argparse import ArgumentParser


def script(chmm_path, h3k27ac_path, out_path, buffer):
    chmm_file = open(chmm_path)
    # Read in all lines at once
    contents = chmm_file.read()
    # Dump the first 2 lines, don't need them
    track = contents.split('\n')[2:]
    logging.info("Found %d enhancer bins total" %
                 sum(1 for x in track if x == '7'))

    # Open the bedGraph for acetylation
    with open(h3k27ac_path) as bed_file:
        enhancers = []
        for line in bed_file:
            tokens = line.strip().split('\t')
            start = int(tokens[1])
            end = int(tokens[2])
            if track[(start + 50)/200] == '7' or track[(end - 50)/200] == '7':
                for i in range(start, end):
                    if track[(i + 50)/200] == '7' or track[(i-50)/200] == '7':
                        enhancers.append((i, float(tokens[3])))

    # We label which enhancer each index belongs to
    indexed_enh = []
    ind = 0
    for i, enh in enumerate(enhancers):
        if enh[0] != enhancers[i-1][0] + 1:
            ind += 1
        indexed_enh.append((enh[0], enh[1], ind))

    # Write our results
    with open(out_path, 'w') as out_file:
        for enh in indexed_enh:
            out_file.write('\t'.join((str(x) for x in enh)) + '\n')

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))


def parse_args(args):
    parser = ArgumentParser(description=__doc__.strip())
    parser.add_argument('chmm_path', metavar='CHMM')
    parser.add_argument('h3k27ac_path', metavar='ENH')
    parser.add_argument('out_path', metavar='OUT')
    parser.add_argument('--buffer', metavar='BUFF', type=int, default=50)
    return parser.parse_args(args)

if __name__ == '__main__':
    sys.exit(main())
