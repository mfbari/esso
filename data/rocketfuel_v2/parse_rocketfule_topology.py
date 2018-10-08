import argparse
import re
import os
import random

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dataset_root", type = str, required = True)
    args = parser.parse_args()
    current_directory = current_directory = os.getcwd()
    directory = os.path.join(current_directory, args.dataset_root)
    alldirs = os.listdir(directory)
    asnos = set()
    for d in alldirs:
        as_no = d.split(":")[0]
        if as_no == "README":
            continue
        asnos.add(as_no)
    for asno in asnos:
        topo_dir = os.path.join(directory, asno + ":" + asno)
        node_map = dict()
        node_count = 0
        edge_count = 0
        edges = []
        max_latency = 0
        if os.path.isdir(topo_dir):
            with open(topo_dir + "/edges.lat") as f:
                lines = f.readlines()
                for line in lines:
                    # 3967:Sunnyvale, CA -> 3967:Santa Clara, CA 0.0288213047287407
                    matches = re.match("\d+:(.*) -> \d+:(.*) (\d+\.*\d*)", line.strip("\r\n"))
                    src_city = matches.group(1)
                    dst_city = matches.group(2)
                    latency = int(round(float(matches.group(3)), 0))
                    max_latency = max(max_latency, latency)
                    if src_city not in node_map:
                        node_map[src_city] = node_count
                        node_count = node_count + 1
                    if dst_city not in node_map:
                        node_map[dst_city] = node_count
                        node_count = node_count + 1
                    u, v = node_map[src_city], node_map[dst_city]
                    if u <= v:
                        edges.append((u, v, latency))
            with open(str(asno) + ".co", "w") as f:
                edge_count = len(edges)
                f.write(str(node_count) + " " + str(edge_count)+ "\n")
                for i in range(0, node_count):
                    f.write(str(i) + " " + str(round(random.uniform(0.8, 1.5), 2)) + " " + str(random.randint(0, 1)) + "\n")
                for edge in edges:
                    f.write(str(edge[0]) + " "+ str(edge[1]) + " "+ str(edge[2]) + " 100000\n")
                    f.write(str(edge[1]) + " "+ str(edge[0]) + " "+ str(edge[2]) + " 100000\n")
        print "\t".join(str(x) for x in [asno, node_count, edge_count, max_latency])



if __name__ == "__main__":
    main()
