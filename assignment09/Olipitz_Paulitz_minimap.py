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
            print("q:" + str(query[1][i_min:i_max]))
            print("t:" + str(sk(targets[t][1], ii_min, ii_max - ii_min, r)))


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

    # return reverse (!) complement
    
    if r == 1:  
        if i != 0:
            return ''.join([comp(c) for c in s[i + k - 1:i - 1:-1]])
        return ''.join([comp(c) for c in s[k - 1::-1]])
            

    # else return substring
    
    return s[i:i+k]
    
def comp(c: chr) -> chr:
    if (c == 'A' or c == 'a'):
        return 'T'
    if (c == 'T' or c == 't'):
        return 'A'
    if (c == 'G' or c == 'g'):
        return 'C'
    if (c == 'C' or c == 'c'):
        return 'G'
    return c

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

    return sum([h_char(c) * 4**i for i,c in enumerate(s[::-1])])

def h_char(c: chr) -> int:
    if (c == 'A' or c == 'a'):
        return 0
    if (c == 'C' or c == 'c'):
        return 1
    if (c == 'G' or c == 'g'):
        return 2
    if (c == 'T' or c == 't'):
        return 3
    return 0

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

    if len(s) >= w + k - 1:

        for i in range(len(s) - w - k + 1):

            m = math.inf

            for j in range(w):

                u = h(sk(s, i + j, k, 0))
                
                v = h(sk(s, i + j, k, 1))
                
                if u != v:
                    
                    m = min(m,u,v)

            for j in range(w):

                u = h(sk(s, i + j, k, 0))
                
                v = h(sk(s, i + j, k, 1))

                if u < v and u == m:

                    M.add((m, i + j, 0))

                elif u > v and v == m:

                    M.add((m, i + j, 1))

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

    for t,st in enumerate(fastA):

        M = minimizerSketch(st[1], w, k)

        for h,i,r in M:
            
            if (h in H.keys()):
                
                H[h] += [(t,i,r)]
                
            else:
                
                H[h] = [(t,i,r)]

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
    
    A = []

    M = minimizerSketch(q, w, k)
    
    for h, i, r in M:
        
        if (h in H.keys()):

            for t, ii, rr in H[h]:
                                    
                if (r == rr):
                    
                    A += [(t,0, i-ii, ii)]
                    
                else:
                    
                    A += [(t,1, i+ii, ii)]

    #(t,r,c,ii)
    A.sort()
    
    b = 0
            
    for e in range(len(A)):
        
        if (e == len(A) - 1 or A[e][0] != A[e+1][0] or A[e][1] != A[e+1][1] 
            or A[e+1][2] - A[e][2] >= eps):

            # in case of different diagonals, start indices are not 
            # monotonically in-/decreasing
            # i, ii are only start positions: add k for end positions
            min_i = min([extract_i(*A[ix]) for ix in range(b,e+1)])
            max_i = max([extract_i(*A[ix]) for ix in range(b,e+1)]) + k
            min_ii = min([A[ix][3] for ix in range(b,e+1)])
            max_ii = max([A[ix][3] for ix in range(b,e+1)]) + k
                
            # current t, current x (or r as it is called here)
            result += [(A[e][0], A[e][1], min_i, max_i, min_ii, max_ii)]
            
            b = e + 1
    
    return result

def extract_i(t,r,c,ii):
    if (r == 1):
        return c - ii
        
    if (r == 0):
        return c + ii

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
