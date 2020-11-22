import pandas as pn

import numpy as np

import networkx as nx

import matplotlib.pyplot as plt

# The given sequences

s1 = "CATG"

s2 = "CAGT"

s3 = "AGTT"

# reprecenting the sequences as an graph

#representing the tree sequences as directed graphs.
Seq1Graph = nx.Graph()
Seq2Graph = nx.Graph()
Seq3Graph = nx.Graph()

for x in s1:
 Seq1Graph.add_node(x)


for x in s2:
 Seq2Graph.add_node(x)

for x in s3:
 Seq3Graph.add_node(x)


# init a graph for concat all graphs with each other.
G = nx.Graph()
H = nx.Graph()
# add edges








print(Seq1Graph.nodes())
print(Seq2Graph.nodes())
print(Seq3Graph.nodes())
print("Nodes of graph: ")
print(G.nodes())
print("Edges of graph: ")
print(G.edges())

#nx.draw(G)
#plt.savefig("path_graph1.png")
#plt.show()
