import networkx as nx
import glob

topo_files = glob.glob("*.co")

for tf in topo_files:
  g = nx.Graph()
  with open(tf) as f:
    node_count, edge_count = [int(x) for x in f.readline().split()]
    for n in range(node_count):
      f.readline()
    for e in range(edge_count*2):
      u, v, lat, cap = [int(x) for x in f.readline().split()]
      g.add_edge(u, v)
  if not nx.is_connected(g):
    print tf.split('.')[0], 'is not connected', node_count
