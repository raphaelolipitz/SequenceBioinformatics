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
    help_str = 'test.py -i <inputfile> -o <outputfile> -a <match score> -b <mismatch score> -d <gap open penalty> -e <gap extension penalty>'
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

# problem 2: check options 
# check input file
if not "ifile" in input_parameters.keys():
    print("Please provide an input file. For help use flag -h")
    sys.exit()
elif not (".txt" in input_parameters["ifile"] or ".fa" in input_parameters["ifile"]):
    print("Input path has to point to a txt or fa file. For help use flag -h")
    sys.exit()
    
# check output file
if not "ofile" in input_parameters.keys():
    print("Please provide an output file. For help use flag -h")
    sys.exit()    
elif not (".txt" in input_parameters['ofile']):
    print("Output file has to be a txt file. For help use flag -h")
    sys.exit()

# check entered values for their correct sign
if ("match_score" in input_parameters.keys()
    and input_parameters["match_score"] < 0):
        print("Warning: Match score should be positive."
              "For help use flag -h")

if ("mismatch_score" in input_parameters.keys()
    and input_parameters["match_score"] > 0):
        print("Warning: Mismatch score should be negative."
              "For help use flag -h")

if ("gap_open" in input_parameters.keys()
    and input_parameters['gap_open'] < 0):
        print("gap penalties are defined positive. Please enter a positive value."
              "For help use flag -h")
        sys.exit()

if ('gap_extend' in input_parameters.keys()
    and input_parameters['gap_extend'] < 0):
        print("gap penalties are defined positive. Please enter a positive value."
              "For help use flag -h")
        sys.exit()

# check mismatch score vs extension penalty
if ("mismatch_score" in input_parameters.keys()
    and 'gap_extend' in input_parameters.keys()
    and input_parameters["mismatch_score"] <= -2 * input_parameters['gap_extend']):
        print("Mismatch score has to be larger than -2 times gap extension penalty."
              "For help use flag -h")
        sys.exit()

# check open penalty vs extension penalty        
if ('gap_open' in input_parameters.keys()
    and 'gap_extend' in input_parameters.keys()
    and input_parameters['gap_open'] < input_parameters['gap_extend']):
        print("Warning: gap open penalty is smaller than gap extension penalty."
              "This will lead to more and small gaps."
              "For help use flag -h")

# read input:

input= "ifile" in input_parameters.keys() and open(input_parameters["ifile"],mode="r") or sys.stdin # read from stdin if no input file given


x_header=input.readline().strip()
x_sequence=input.readline().strip() # assuming whole sequence is on one line
y_header=input.readline().strip()
y_sequence=input.readline().strip() # assuming whole sequence is on one line

input==sys.stdin or input.close() # close if reading from a file

# problem 1: compute alignment score 

if not 'match_score' in input_parameters.keys():
    print("Match score not provided. Will use default value of 1.")
    input_parameters["match_score"] = 1

if not 'mismatch_score' in input_parameters.keys():
    print("Mismatch score not provided. Will use default value of -2.")
    input_parameters["mismatch_score"] = -2

if not 'gap_open' in input_parameters.keys():    
    print("Gap open penalty not provided. Will use default value of 3.")
    input_parameters["gap_open"] = 3

if not 'gap_extend' in input_parameters.keys():
    print("Gap extension penalty not provided. Will use default value of 2.")
    input_parameters["gap_extend"] = 2

d = input_parameters["gap_open"]
e = input_parameters["gap_extend"]

def aff_gap_pen(k):
    return input_parameters["gap_open"] + (k - 1) * input_parameters["gap_extend"]

def match(x,y):
    if (x==y):
        return input_parameters["match_score"]
    return input_parameters["mismatch_score"]

# initialization
rows, cols = (len(x_sequence), len(y_sequence))

# Used pandas because lists have an annoying behaviour: 
# If you define the list of lists like this: [[-inf] * (cols + 1)] * (rows + 1)
# all the rows are actually one list and will always be the same. 
# Also, one can easily pretty print matrices
M = pd.DataFrame([[-inf] * (cols + 1)] * (rows + 1))
M.loc[0, 0] = 0

I = pd.DataFrame([[-inf] * (cols + 1)] * (rows + 1))
I.loc[1:,0] = [-d - (i - 1) * e for i in range(1, rows + 1)]
I.loc[0,1:] = [-d - (i - 1) * e for i in range(1, cols + 1)]
    
#print(pd.DataFrame(M))
#print(pd.DataFrame(I))

# recursion

for i in range(1, rows + 1):
    for j in range(1, cols + 1):

        M.loc[i, j] = max(M.loc[i - 1, j - 1] + match(x_sequence[i-1], y_sequence[j-1]), 
                      I.loc[i - 1, j - 1] + match(x_sequence[i-1], y_sequence[j-1]))

        I.loc[i, j] = max(M.loc[i-1, j] - d,
                      I.loc[i-1, j] - e,
                      M.loc[i, j-1] - d,
                      I.loc[i, j-1] - e)


optimalScore = max( M.loc[rows, cols], I.loc[rows, cols])

print("optimal alignment score: {}".format(optimalScore))
#print(pd.DataFrame(M))
#print(pd.DataFrame(I))

# problem 3: perform trace-back

# Making two empty strings for the optimal alingment.

alignedSequenceX = ""

alignedSequenceY = ""

i, j = (rows, cols)
are_at_M = M.loc[i, j] >= I.loc[i, j]
while (i > 0 or j > 0):
    
    if are_at_M:
        # we have a match or substitution event
        alignedSequenceX = x_sequence[i-1] + alignedSequenceX
        alignedSequenceY = y_sequence[j-1] + alignedSequenceY
        
        # wont give out of bounds because we are never at M if one of i and j are 0
        are_at_M = M.loc[i-1, j-1] >= I.loc[i-1, j-1]
        i -= 1
        j -= 1
        
    else:
        # are we coming from M or I?
        are_at_M = (max(M.loc[i-1, j], M.loc[i, j-1]) - d >= 
                    max(I.loc[i-1, j], I.loc[i, j-1]) - e)
        
        # is the gap in X or Y?
        if (are_at_M):
            gap_in_x = M.loc[i-1, j] < M.loc[i, j-1]
        else:
            gap_in_x = I.loc[i-1, j] < I.loc[i, j-1]
            
        if (gap_in_x):
            alignedSequenceX = '-' + alignedSequenceX
            alignedSequenceY = y_sequence[j-1] + alignedSequenceY
            j -= 1
        else:
            alignedSequenceX = x_sequence[i-1] + alignedSequenceX
            alignedSequenceY = '-' + alignedSequenceY
            i -= 1
        
output= "ofile" in input_parameters.keys() and open(input_parameters["ofile"],mode="w") or sys.stdout # write to stdout if no output file given

print(alignedSequenceX,alignedSequenceY,"Final Score: " + str(optimalScore), sep="\n",file=output)

output == sys.stdout or output.close() # close if writing to a file




