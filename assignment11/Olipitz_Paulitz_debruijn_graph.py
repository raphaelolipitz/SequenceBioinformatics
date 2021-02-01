# Sequence bioinformatics, WS 20/21, Daniel Huson

import sys
from optparse import OptionParser
import networkx as nx
import matplotlib.pyplot as plt
from typing import Tuple, Dict, Set

__author__ = "Your name here"


def main():
    """
        computes a de-Bruijn graph using the naive algorithm

        Usage: debruin_graph.py [options] input-file

        Options:
            -h, --help              show this help message and exit
            -k K                    k-mer size

     """

    parser = OptionParser("%prog infile", description="Naive de-Bruijn graph", epilog="Author(s): " + __author__)

    parser.add_option("-k", default=4, action="store", dest="k", help="k-mer size")
    parser.add_option("-a", default=False, action="store_true", dest="all",
                      help="All-against-all k-mer overlaps (not just per read)")
    parser.add_option("-s", default=False, action="store_true", dest="simplify", help="Simply graph")

    (options, args) = parser.parse_args()

    if len(args) != 1:
        sys.exit("Must specify fastA file, got: " + str(len(args)) + " arguments")

    reads_file = args[0]
    k = int(options.k)

    reads = read_fastA(reads_file)

    children = dict()  # maps each k-mer to its children in the graph

    for read in reads:
        process_read(read[1], k, children)

    if not options.all:
        print("read-based:")
        print_graph(children)

    if options.all:
        add_all_others(k, children)

    if options.all:
        print("read-based and others:")
        print_graph(children)

    if options.simplify:
        simplify_graph(k, children)

    if options.simplify:
        print("simplified:")
        print_graph(children)

    # stupid work around:
    # nx.spring_layout() requires an undirected graph to produce a good layout
    # nx.draw_networkx_edges() requires a directed graph if you want to display edges as arrows
    # so need to construct both an undirected and directed graph and the same time,
    # compute node positions for undirected graph and then use them for directed graph:
    graph = nx.Graph()
    di_graph = nx.DiGraph()

    # 1. construct graph and di_graph (using graph methods add_node() and add_edge()
    for u in children.keys():
        graph.add_node(u)
        di_graph.add_node(u)
        for v in children[u]:
            graph.add_edge(u,v)
            di_graph.add_edge(u,v)
            
    # 2. compute node positions for graph (using nx layout method)
    node_pos = nx.spring_layout(graph)
    #edge_pos = nx.draw_networkx_edges(di_graph, node_pos)
    
    # 3. draw edges, nodes and labels for di_graph using positions computed for graph (using draw_networkx methods)
    nx.draw_networkx(di_graph, node_pos, font_color = 'red', node_size = 20, node_color = 'white', font_size = 6)
    
    # 4. show the plot (use the method plt.waitforbuttonpress() for window to remain open)
    plt.waitforbuttonpress()

def process_read(seq: str, k: int, children: Dict[str, Set[str]]) -> None:
    """processes a read

        Parameters
        ----------
        seq : str
           The sequence
        k : int
            K-mer size
        children: Dict[str, Set[str]]
            Maps each k-mer to the set of all k-mers that direct successor in graph

        Returns
        -------
        None
       """
       
    for i in range(len(seq) - k):
        if (seq[i:i+k] in children):
            children[seq[i:i+k]].add(seq[i + 1:i + k + 1])
        else: 
            children[seq[i:i+k]] = set([seq[i + 1:i + k + 1]])
            
    if (not seq[-k:] in children):
        children[seq[-k:]] = set()


def add_all_others(k: int, children: Dict[str, Set[str]]) -> None:
    """Finds all other children for a given k-mer (that are not necessarily from the same read):

        Parameters
        ----------
        k : int
            K-mer size
        children: Dict[str, Set[str]]
            Maps each k-mer to the set of all its children in graph

        Returns
        -------
        None
       """
    
    alphabet = ['A', 'G', 'C', 'T']
    
    for key in children:
        for base in alphabet:
            if key[1:] + base in children:
                children[key].add(key[1:] + base)

def simplify_graph(k: int, children: Dict[str, Set[str]]) -> None:
    """Simplifies the graph

        Parameters
        ----------
        k : int
            K-mer size
        children: Dict[str, Set[str]]
            Maps each k-mer to the set of all its children in graph

        Returns
        -------
        None
       """
       
    # repeat merging until no more nodes can be merged
    has_changed = True
    while(has_changed):
        has_changed = False
        
        # iterate over all keys to look for possible merges
        for key in children:
            
            # only bother with keys having a single child
            if (len(children[key]) != 1):
                continue
        
            
            [child] = children[key]
            
            # do not merge with itself to prevent infinite loop
            if (child == key):
                continue
            
            only_parent = True
            
            # make sure the child does not have any other parents
            for key2 in children:
                if (key != key2 and child in children[key2]):
                    only_parent = False
                    break
            
            # merge key and its child, then restart for-loop
            if (only_parent):
                new_key = key + child[k-1:]

                # merge key and its child into new entry
                children[new_key] = children[child]
                
                # make sure to update occurrences of the key in other keys
                for key2 in children:
                    if key in children[key2]:
                        children[key2].add(new_key)
                        children[key2].remove(key)
                
                # remove entry of key and its child
                children.pop(child)
                children.pop(key)      
                
                has_changed = True
                break
     

def print_graph(children: Dict[str,Set[str]]) -> None:
    """ print the graph to the console in this format:
        kmer -> children in alphabetical order
        ...
        nodes: <number of nodes>
        edges: <number of edges>

                        Parameters
                        ----------
                         children: Dict[str, Set[str]]
                            Maps each k-mer to the set of all its children in graph

                        Returns
                        -------
                        None
                       """
                       
    keys = list(children.keys())
    keys.sort()
    for key in keys:
        ls = list(children[key])
        ls.sort()
        print("{} -> {}".format(key, ls))
    
    print("nodes: {}".format(len(children.keys())))
    print("edges: {}".format(sum([len(children[key]) for key in children.keys()])))


def read_fastA(filename: str) -> [Tuple[str, str]]:
    """Gets list of headers and sequences in fastA format

        Parameters
        ----------
        filename : str
            The file location of the spreadsheet

        Returns
        -------
        list
            a list of tuples, each containing a header and sequence
        """
    ins = open(filename)

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

