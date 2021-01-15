import math
from optparse import OptionParser
from typing import Tuple, Set, Dict, List

# Sequence bioinformatics, WS 20/21, Daniel Huson

__author__ = "Emil Paulitz and Raphael Olipitz"


def main():
    """
     Proof-of-concept implementation of minimap

   Usage: minimap.py [options] target-file query-file

    Options:
        -h, --help              show this help message and exit
        -k K                    k-mer size
        -w W                    window size
    """

    parser = OptionParser("%prog [options] infile", description="Proof-of-concept implementation of minimap",
                          epilog="Author(s): " + __author__)

    parser.add_option("-k", default=15, action="store", dest="k", help="k-mer size")
    parser.add_option("-w", default=10, action="store", dest="w", help="window size")
    parser.add_option("-c", default=False, action="store_true", dest="check",
                      help="check sub sequence extraction, reverse complement and hash computation")
    (options, args) = parser.parse_args()

    if len(args) != 2:
        raise IOError("Must specify target-file and query-file, got:", args)

    target_file = args[0]
    query_file = args[1]

    if options.check:  # check that hash function, sub sequence extraction and reverse complement appear to be correct:
        s = 'CACGGTAGA'
        if h(sk(s, 1, 5, 0)) != 107:  # hash value of ACGGT
            print(sk(s, 1, 5, 0), "=", h(sk(s, 1, 5, 0)), "should be: 107")  # result should be 107
        if h(sk(s, 1, 5, 1)) != 91:  # hash value of ACCGT (reverse complement of ACGGT)
            print(sk(s, 1, 5, 1), "=", h(sk(s, 1, 5, 1)), "should be 91")  # result should be 91

    k = options.k
    w = options.w

    targets = read_fastA(target_file)
    queries = read_fastA(query_file)

    for tar in targets:
        print("target:", tar[0], "length=", len(tar[1]))

    H = index(targets, w, k)

    print("Index size:", len(H))

    for query in queries:
        print("Query:", query[0])

        matches = map(H, query[1], w, k, 500)
        print("Matches:", len(matches))

        for t, r, i_min, i_max, ii_min, ii_max in matches:
            print("match:", " target=", t + 1, "range:", i_min + 1, "-", i_max + 1, "to", ii_min + 1, "-", ii_max + 1,
                  "reverse:", r)
            print(query[1][i_min:i_max])
            print(sk(targets[t][1], ii_min, ii_max - ii_min, r))



def sk(s: str, i: int, k: int, r: int) -> str:
    """Gets the sub-sequence of s that starts at position i and has length k

            Parameters
            ----------
            s : str
                The sequence
            i : int
                Start pos
            k : int
                Length
            r : int
                if 1, return reverse complement

            Returns
            -------
            str
                the sub sequence
            """

    # please implement


    String_to_return = ""

    #return supstring
    if r == 0:

        for chr in range(i , k-1):

            String_to_return += s[chr]

    #return reverse complement
    if r == 1:

        for chr in range(i, k - 1):

            if s[chr] == "A":

                String_to_return += "T"

            if s[chr] == "T":
                String_to_return += "A"

            if s[chr] == "C":
                String_to_return += "G"

            if s[chr] == "G":
                String_to_return += "C"




    return String_to_return


def h(s: str) -> int:
    """Gets the hash value of a DNA sequence

                Parameters
                ----------
                s : str
                    The sequence


                Returns
                -------
                int
                    the none-negative hash value
                """

    # please implement

    #hashvalue of the DNA bases.

    hash_value_of_A = 0
    hash_value_of_C = 1
    hash_value_of_G = 2
    hash_value_of_T = 3

    intString = ""
    hash_value = 0

    #computes int complement of the string.

    for chr in s:

        if chr == "A":

            intString += str(hash_value_of_A)

        if chr == "C":

            intString += str(hash_value_of_C)


        if chr == "G":

            intString += str(hash_value_of_G)

        if chr == "T":

            intString += str(hash_value_of_T)


    k = len(intString)

    while k > 0:

        hash_value += int(intString[k]) * (4**(k - 1))






    return hash_value


def minimizerSketch(s: str, w: int, k: int) -> Set[Tuple[int, int, int]]:
    """computes all minimizers for a sequence

                    Parameters
                    ----------
                    s : str
                        The sequence
                    w : int
                        Window size
                    k : int
                        k-mer size

                    Returns
                    -------
                    Set[Tuple[int,int,int]]
                        Set of minimizers (h,i,r), where h is hash value, i is position and r is strand
                    """
    M = set()

    # please implement

    return M


def index(fastA: [Tuple[str, str]], w: int, k: int) -> Dict[int, List[Tuple[int, int, int]]]:
    """computes the hash index of minimizers for a given list of sequences

                     Parameters
                     ----------
                     fastA : List[Tuple[str, str]]
                         List of fastA records, each consisting of a header and sequence
                     w : int
                         Window size
                     k : int
                         k-mer size

                     Returns
                     -------
                     Dict[int, List[Tuple[int, int, int]]]:
                         Dictionary mapping hash values to (t,i,r) where t is the number of the target sequence,
                         i is the position in the sequence and r is the strand
                     """

    H = dict()

    # please implement

    return H


def map(H: Dict[int, List[Tuple[int, int, int]]], q: str, w: int, k: int, eps: int) -> List[
    Tuple[int, int, int, int, int, int]]:
    """ maps a query to target sequences represented in a given hash table

                     Parameters
                     ----------
                     H : Dict[int, List[Tuple[int, int, int]]]
                        Dictionary mapping hash values to (t,i,r) where t is the number of the target sequence,
                        i is the position in the sequence and r is the strand
                     q : str
                        Query sequence
                     w : int
                         Window size
                     k : int
                         k-mer size
                     eps: int
                        Epsilon band width

                     Returns
                     -------
                     List[Tuple[int, int, int, int, int, int]]:
                         List of tuples (t,r,i_min,i_max,ii_min,ii_max) where t is the number of the target sequence,
                         r is the relative strand, i_min and i_max are the positions in the query and
                         ii_min and ii_max are the positions in the target sequence
                     """

    result = []

    # please implement

    return result


def read_fastA(file_loc: str) -> [Tuple[str, str]]:
    """Gets list of headers and sequences in fastA format

        Parameters
        ----------
        file_loc : str
            The file location of the spreadsheet

        Returns
        -------
        list
            a list of tuples, each containing a header and sequence
        """
    ins = open(file_loc)

    records = []

    header = ""
    seq = ""
    for line in ins:
        line = line.strip()
        if line.startswith(">"):
            if header != "":
                records.append((header, seq))
            header = line
            seq = ""
        else:
            seq += line.replace("\\s+", "")
    if header != "":
        records.append((header, seq))
    ins.close()
    return records


if __name__ == '__main__':
    main()
