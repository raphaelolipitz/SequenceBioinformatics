#!/usr/bin/python
import sys, getopt
from cmath import inf
import pandas as pd

"""assignment_02.py: Pairwise alignment."""

__author__      = "Emil und Raphael"

# This function parses commandline options, as follows:
#  -i <input-file>: if provided, read sequences from named file, not from standard input.
#  -o <output-file>: If provided, write aligned sequences to named file, not to standard
# output.
#  -a <score>: Set match score.
#  -b <score>: Set mismatch score.
#  -d <score>: Set gap open penalty.
#  -e <score>: Set gap extend penalty.

def parse_commandline(argv):
    help_str = 'test.py -i <inputfile> -o <outputfile> -a <match score> -b <mismatch score> -c <gap penalty>'
    input_parameters = {}
    try:
        opts, args = getopt.getopt(argv, "h:i:o:a:b:d:e:",
                                   ["ifile=",
                                    "ofile=",
                                    "match_score=",
                                    "mismatch_score=",
                                    "gap_open=",
                                    "gap_extend="])
    except getopt.GetoptError:
        print(help_str)
        sys.exit()

    for opt, arg in opts:
        if opt == '-h':
            print(help_str)
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_parameters["ifile"] = arg

        elif opt in ("-o", "--ofile"):
            input_parameters["ofile"] = arg

        elif opt in ("-a", "--match_score"):
            input_parameters["match_score"] = int(arg)

        elif opt in ("-b", "--mismatch_score"):
            input_parameters["mismatch_score"] = int(arg)

        elif opt in ("-d", "--gap_open"):
            input_parameters["gap_open"] = int(arg)

        elif opt in ("-e", "--gap_extend"):
            input_parameters["gap_extend"] = int(arg)

    print("Options:")
    for i, j in input_parameters.items():
        print("\t",i,": ", j,sep="")
    return input_parameters

# The program starts here:
print("Program: Assignment 02 - pairwise aligner with affine gap scoring")
print("Author:",__author__)

# parse options:
input_parameters=parse_commandline(sys.argv[1:])

# problem 2: check options PLEASE_IMPLEMENT
print("TODO: check parameters")

# read input:
#input= "sequences.txt" in input_parameters.keys() and open(input_parameters["sequences.txt"],mode="r") or sys.stdin # read from stdin if no input file given

input= "ifile" in input_parameters.keys() and open(input_parameters["ifile"],mode="r") or sys.stdin # read from stdin if no input file given


x_header=input.readline().strip()
x_sequence=input.readline().strip() # assuming whole sequence is on one line
y_header=input.readline().strip()
y_sequence=input.readline().strip() # assuming whole sequence is on one line

input==sys.stdin or input.close() # close if reading from a file

# problem 1: compute alignment score PLEASE_IMPLEMENT
print("TODO: perform dynamic programming to compute alignment score")
# set the gap-penalties

matchScore = 1

mismatchScore = -2

gap_open = 3

gap_extend = 2

def aff_gap_pen(k):
    return gap_open + (k - 1) * gap_extend

def match(x,y):
    if (x==y):
        return matchScore
    return mismatchScore

# initialization
rows, cols = (len(x_sequence), len(y_sequence))

# Used pandas because lists have a really dumb behaviour: 
# If you define the list of lists like this: [[-inf] * (cols + 1)] * (rows + 1)
# all the rows are actually one list and will always be the same. 
# To avoid behaviour like this, I used pandas, which has a kind of strange 
# cell selection syntax but works great once you have got the hang of it
M = pd.DataFrame([[-inf] * (cols + 1)] * (rows + 1))
M.loc[0, 0] = 0

I = pd.DataFrame([[-inf] * (cols + 1)] * (rows + 1))
I.loc[1:,0] = [gap_open + (i - 1) * gap_extend for i in range(1, rows + 1)]
I.loc[0,1:] = [gap_open + (i - 1) * gap_extend for i in range(1, cols + 1)]
    
#print(pd.DataFrame(M))
#print(pd.DataFrame(I))

# recursion

for i in range(1, rows + 1):
    for j in range(1, cols + 1):

        M.loc[i, j] = max(M.loc[i - 1, j - 1] + match(x_sequence[i-1], y_sequence[j-1]), 
                      I.loc[i - 1, j - 1] + match(x_sequence[i-1], y_sequence[j-1]))

        I.loc[i, j] = max(M.loc[i-1, j] - gap_open,
                      I.loc[i-1, j] - gap_extend,
                      M.loc[i, j-1] - gap_open,
                      I.loc[i, j-1] - gap_extend)


optimalScore = max( M.loc[rows, cols], I.loc[rows, cols])

print("optimal alignment score: {}".format(optimalScore))
#print(pd.DataFrame(M))
#print(pd.DataFrame(I))
# problem 3: perform trace-back PLEASE_IMPLEMENT
print("TODO: perform trace-back to compute alignment")

# Making two empty strings for the optimal alingment.

alignedSequenceX = ""

alignedSequenceY = ""

# Loop for the Traceback. Goes from last cell to first and concatinates the two sequences.
while rows > 0 or cols > 0:

    # if anweisung für den Fall eines Matches.
    if rows > 0 and j > 0 and Ix[i - 1] == Iy[j - 1]:
        align_X = Ix[rows - 1] + align_X
        align_Y = Iy[cols - 1] + align_Y
        rows = rows - 1
        cols = cols - 1

    # Fall das die X Sequence einen gap bekommt.
    elif cols > 0 and current_score == left_score + penalty:
        align_X = Ix[rows - 1] + align_X
        alignedSequenceY = "-" + alignedSequenceY
        rows = rows - 1

    # Fall das die Y Sequence einen gap bekommt.
    else:
        alignedSequenceX = "-" + alignedSequenceX
        align_Y = Iy[cols - 1] + align_Y
        cols = cols - 1


output= "ofile" in input_parameters.keys() and open(input_parameters["ofile"],mode="w") or sys.stdout # write to stdout if no output file given

# This currently writes the two sequences to the given output file... replace this by the output that your algorithm computes
print("TODO: replace output of sequences by output of score and output of alignment")

print(alignedSequenceX,alignedSequenceY,optimal, sep="\n",file=output)

output == sys.stdout or output.close() # close if writing to a file




