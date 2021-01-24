# Sequence bioinformatics, WS 20/21, Daniel Huson

from optparse import OptionParser
from typing import Tuple, List
from parse_header import parse_header

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.setrecursionlimit(10**6)

__author__ = "Emil Paulitz and Raphael Olipitz"


def main():
    """
     Draw heat map for virus alignment header lines

   Usage: YourNameHere_heatmap.py [options] input-file

    Options:
        -h, --help              show this help message and exit
    """

    parser = OptionParser("%prog infile", description="Draw heatmap for virus alignment headers",
                          epilog="Author(s): " + __author__)

    options, args = parser.parse_args()

    if len(args) != 1:
        raise IOError("Must specify input file, got:", args)

    file_loc = args[0]

    taxa, proteins, table = parse_header_lines(file_loc)

    if (len(taxa) != 0 and len(proteins) != 0):
        plt.figure(dpi = (max(len(taxa), len(proteins)) * 6))
        maxLenLabels = max([len(label) for label in taxa + proteins])
        FontDeterminingVar = max(len(taxa), len(proteins), maxLenLabels/10)
        plt.rcParams.update({'font.size': min(18, 72/FontDeterminingVar)})
        plt.pcolor(table)
        plt.yticks(np.arange(0.5, len(taxa), 1), taxa)
        plt.xticks(np.arange(0.5, len(proteins), 1), proteins, rotation=10, ha='right')
        for i in range(len(taxa)):
            for j in range(len(proteins)):
                if (table[i, j] == 0):
                    num = 0
                elif (round(table[i, j],2) == 0):
                    num = 0.01
                else:
                    num = round(table[i, j],2)
                plt.text(j + 0.5, i + 0.5, num, ha="center", va="center", color="w")
        plt.tight_layout()
        savename = './heatmap{}.svg'.format(''.join([c for c in file_loc if c.isnumeric()]))
        plt.savefig(savename)#, dpi = 'figure')
        #plt.show()
    else:
        print('No taxa or no proteins found. Check input.')


def parse_header_lines(file_loc: str) -> Tuple[List[str], List[str], np.array]:
    """Parses the input lines

                Parameters
                ----------
                file_loc : str
                    The file location

                Returns
                -------
                Tuple[List[str], List[str], np.array]
                    the taxon names, protein names and percentage matrix
                """

    # dict of taxon-dict of protein-names to store the number of occurences
    occ = {}

    proteins_set = set()

    with open(file_loc, 'r') as file_obj:
        for line in file_obj.readlines():

            _, taxon, protein = parse_header(line)
            proteins_set.add(protein)

            if (taxon in occ.keys()):
                if (protein in occ[taxon]):
                    occ[taxon][protein] += 1
                else:
                    occ[taxon][protein] = 1
            else:
                occ[taxon] = {}
                occ[taxon][protein] = 1

    total_sum = sum([occ[t][p] for t in occ.keys() for p in occ[t].keys()])

    # fix order of taxa and proteins
    taxa = list(occ.keys())
    proteins = list(proteins_set)

    perc = []
    for taxon in taxa:
        taxon_list = []
        for protein in proteins:
            if (protein in occ[taxon].keys()):
                taxon_list.append(occ[taxon][protein]/total_sum)
            else:
                taxon_list.append(0)
        perc.append(taxon_list)

    # please implement
    return (taxa, proteins, np.array(perc))

if __name__ == '__main__':
    main()

