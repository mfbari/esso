import argparse
import fnss
import random
import networkx as nx

middleboxes = ['0', '1', '2', '3', '4', '5']
def GenerateMiddleboxSequence(n):
  random.shuffle(middleboxes)
  return middleboxes[0:n]
  
def read_co_topology(co_file):
    g = nx.Graph()
    with open(co_file, "r") as f:
        line = f.readline()
        tokens = map(int, line.strip("\n\r").split(","))
        num_nodes = tokens[0]
        for i in range(0, num_nodes):
            line = f.readline()
        num_edges = tokens[1]
        for line in f:
            tokens = line.strip("\n\r").split(",")
            u, v, lat, bw = int(tokens[0]), int(tokens[1]), float(tokens[2]), int(tokens[3])
            g.add_edge(u, v, {'bw':int(bw), 'latency':int(lat)})
    return g

def get_timeslot_index(ts, slot_width):
    return (ts / slot_width)

def generate_traffic_matrix():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--as_no", type = str, required = True)
    parser.add_argument('--arrival_rate', type=float, default=0.003) # 10 sfc per hour
    parser.add_argument('--max_simulation_time', type=int, default=86400) # 24 hour
    parser.add_argument('--sfc_lifetime', type=int, default=5400) # 1 hour 30 minute
    parser.add_argument('--num_time_slots', type=int, default=24)

    args = parser.parse_args()
    co_file = args.as_no + ".co"
    co_topology = read_co_topology(co_file)
    slot_width = args.max_simulation_time / args.num_time_slots

    times = []
    current_time = 0
    while current_time < args.max_simulation_time:
        wait_time = int(random.expovariate(args.arrival_rate)) + 1
        duration = int(random.expovariate(1.0 / args.sfc_lifetime)) + 1
        times.append((current_time, current_time + duration))
        current_time += wait_time
    timeslot_sfc_count = dict()
    flow_time_slot = dict()
    sources = random.sample(range(0, len(co_topology.nodes()) - 1), len(co_topology.nodes())/2)
    destinations = random.sample(range(0, len(co_topology.nodes()) - 1), len(co_topology.nodes())/2)
    topology = fnss.Topology(co_topology)
    fnss.set_capacities_constant(topology, 10000, 'Mbps')
    tmc = fnss.sin_cyclostationary_traffic_matrix(
           topology, 
           mean=0.1, # average flow in TM is 100 Mbps 
           stddev=0.03, # this is the std among average flows of different OD pairs 
           gamma=0.8,     # gamma and log_psi are parameters for fitting the std of 
           log_psi=-0.033, # volume fluctuations over time. Look at Nucci et al. paper
           delta=0.2, # traffic variation from period max and avg as fraction of average
           n=1, # number of samples per each period
           periods=args.num_time_slots, # number of periods
           max_u=0.50, # max link utilization desired 
           origin_nodes=sources, 
           destination_nodes=destinations)
    flow_id = 0
    all_flows = []
    for i in range(len(tmc)):
      tm = tmc.get(i)
      flows = tm.flows()
      for (endpoint, traffic) in flows.items():
        [source, destination] = endpoint
        traffic *= 1000
        traffic /= round((24*1.0/args.num_time_slots) * 3600.0, 0)
        if source == destination:
          continue
        if traffic < 10:
          continue
        if traffic > 8000:
          traffic = 8000
        mbox_seq = GenerateMiddleboxSequence(random.randint(3, 6))
        flow = dict()
        flow["flow_id"] = flow_id
        flow["src"] = source
        flow["dst"] = destination
        flow["bw"] = int(round(traffic, 0))
        flow["mbox_seq"] = mbox_seq
        flow["delay"] = int(round(random.uniform(1000, 2000), 0))
        all_flows.append(flow)
        flow_id += 1 
    # traffic_file.write(str(current_time) + "," + str(source) + "," +
    # str(destination) + "," + str(int(round(traffic, 0))) + "," + str(max_latency) + 
    # ",0.00000010," + mbox_seq[0] + "," + mbox_seq[1] + "," + mbox_seq[2] +
    # "\n")
    num_flows = len(all_flows)
    with open(args.as_no + ".timeslots", "w") as f:
        f.write(str(len(times)) + "\n")
        f.write(str(args.num_time_slots) + "\n")
        flow_index = 0
        for time_point in times:
            current_time, ends_at = time_point
            cur_slot_idx, end_slot_idx = get_timeslot_index(current_time, slot_width), get_timeslot_index(ends_at, slot_width)
            if cur_slot_idx not in timeslot_sfc_count.keys():
                timeslot_sfc_count[cur_slot_idx] = 0
                flow_time_slot[cur_slot_idx] = []
            timeslot_sfc_count[cur_slot_idx] += 1
            ttl = end_slot_idx - cur_slot_idx + 1
            flow_idx = random.randint(0, num_flows - 1)
            flow = all_flows[flow_idx]
            flow['ttl'] = ttl
            flow_str = " ".join([str(flow_index), str(flow["src"]),str(flow["dst"]), str(flow["ttl"]), str(len(flow["mbox_seq"])), " ".join(flow["mbox_seq"]), str(flow["bw"]), str(flow["delay"])])
            flow_time_slot[cur_slot_idx].append(flow_str)
            flow_index += 1
        for (ts, flows) in flow_time_slot.iteritems():
            f.write(str(len(flows)) + "\n")
            for flow in flows:
                f.write(flow + "\n")

if __name__ == "__main__":
    generate_traffic_matrix()
    



