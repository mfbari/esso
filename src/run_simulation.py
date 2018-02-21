#!/usr/bin/python

import sys
import os
import argparse
import subprocess
import shutil
from timeit import default_timer as timer
from collections import defaultdict

def int_or_float(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

class essobject:
    pass

class esso_sfc:
    def __init__(self):
        self.cec = -1.0 # current embedding cost
        pass

    # reads data for an SFC from a stream that 
    # provides a readline() function
    def read(self, strm):
        self.sfc_data = [int(x) for x in strm.readline().split()]
        self.sfc_data[5:-2] = [int(vnf_flavor_to_cpu[int(x)]['cpu_count'])
                for x in self.sfc_data[5:-2]]

    def ttl(self):
        return self.sfc_data[3]

    def id(self):
        return self.sfc_data[0]

    def dec_ttl(self):
        self.sfc_data[3] -= 1
        if self.sfc_data[3] < 0:
            raise ValueError('Negative TTL')

    def __str__(self):
        return " ".join([str(x) for x in self.sfc_data])

class esso_object:
    def __init__(self, data):
        self.data = data.strip()

    def __str__(self):
        return self.data

class esso_server:
    cpu_idx = 5
    def __init__(self, data):
        self.data_val = data.split()
        self.data_val[esso_server.cpu_idx] = int(
                self.data_val[esso_server.cpu_idx])

    def dec_cpu_count(self, val):
        self.data_val[esso_server.cpu_idx] -= val

    def __str__(self):
        return " ".join([str(x) for x in self.data_val])

class esso_edge:
    bw_idx = 4
    def __init__(self, data):
        self.data_val = data.split()
        self.data_val[esso_edge.bw_idx] = int(self.data_val[esso_edge.bw_idx])

    def dec_bandwidth(self, val):
        self.data_val[esso_edge.bw_idx] -= val

    def __str__(self):
        return " ".join(str(x) for x in self.data_val)

vnf_flavor_to_cpu = {}
sfcs = []
sfc_in = [set()]
sfc_out = [set()]
sfc_count = 0
timeslot_count = 0
co_list = []
node_list = []
edge_list = []
edge_dir = defaultdict(lambda: defaultdict(int))

def update_write_topology_file(t, mapping):
    if mapping:
        mv = [int_or_float(x) for x in mapping.split()]
        if mv[0] == 200:
            vc = mv[5]
            cpus = mv[6:6+vc]
            servers = mv[10+vc:10+2*vc]
            bw = mv[6+vc]
            paths_data = mv[10+2*vc:]
            for i in range(vc):
                node_list[servers[i]].dec_cpu_count(cpus[i])
            itr = iter(paths_data)
            path_count = itr.next()
            for p in range(path_count):
                edge_count = itr.next()
                u = itr.next()
                for e in range(edge_count-1):
                    v = itr.next()
                    edge_list[edge_dir[u][v]].dec_bandwidth(bw)
                    u = v

    with open('res_topology_'+str(t)+'.dat', 'w') as f:
        f.write(str(len(co_list)) + '\n')
        for co in co_list:
            f.write(str(co) + '\n')
        f.write(str(len(node_list)) + ' ' + str(len(edge_list)) + '\n')
        for node in node_list:
            f.write(str(node) + '\n')
        for edge in edge_list:
            f.write(str(edge) + '\n')

def read_topology_file(dataset_path):
    with open(os.path.join(dataset_path, 'init_topology.dat')) as f:
        co_count = int(f.readline())
        for c in range(co_count):
            co = esso_object(f.readline())
            co_list.append(co)
        node_count, edge_count = [int(x) for x in f.readline().split()]
        for n in range(node_count):
            line = f.readline()
            values = line.split()
            if values[1] == "c":
                node_list.append(esso_server(line))
            else:
                node_list.append(esso_object(line))
        for e in range(edge_count):
            line = f.readline()
            values = line.split()
            u, v = values[1:3]
            edge_dir[u][v] = len(edge_list)
            edge_list.append(esso_edge(line))

def read_vnf_types_file(dataset_path):
    with open(os.path.join(dataset_path, "vnf_types.dat")) as f:
        type_count = int(f.readline())
        # skip over the types
        for i in range(type_count): f.readline()
        # read the flavors
        flavor_count = int(f.readline())
        for line in f.readlines():
            fid, tid, cpu, pd = [int(x) for x in line.split()]
            vnf_flavor_to_cpu[fid] = {'type_id':tid,
                    'cpu_count':cpu,
                    'proc_delay': pd}

def read_timeslots_file(dataset_path):
    global sfc_count
    global timeslot_count
    global sfc_in
    global sfc_out
    with open(os.path.join(dataset_path, "timeslots.dat")) as f:
        sfc_count = int(f.readline())
        timeslot_count = int(f.readline())
        sfc_in = [set() for t in range(timeslot_count)]
        sfc_out = [set() for t in range(timeslot_count)]
        for t in range(timeslot_count):
            t_sfc_count = int(f.readline())
            for s in range(t_sfc_count):

                ###############################
                #values = [x for x in f.readline().split()]
                # convert from flavor-id to cpu-count
                #values[5:-2] = [str(vnf_flavor_to_cpu[int(x)]['cpu_count'])
                #                    for x in values[5:-2]]
                #ttl = int(values[3])
                #sfc_out[t+ttl].add(sfc_id)
                ###############################

                sfc = esso_sfc()
                sfc.read(f)
                sfc_in[t].add(sfc.id())
                sfc_out[t+sfc.ttl()].add(sfc.id())

                sfcs.append(sfc)

if __name__ == '__main__':

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help = "path to dataset folder")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--cplex', action='store_true',
            help="run CPLEX code")
    group.add_argument('-t', '--tabusearch', action='store_true',
            help="run tabu search code")
    parser.add_argument('-i', '--id', required=True,
            help="unique id for this run")
    parser.add_argument('-r', '--replace', action='store_true',
            help="replace existing run dir")
    parser.add_argument('-d', '--dryrun', action='store_true',
            help="only generate data")
    parser.add_argument('-m', '--migthr', type=float, default=0.3,
            help='migration threshold')
    args = parser.parse_args()

    # set the migration_threshold
    migration_threshold = args.migthr

    # path to the dataset and run folder
    dataset_path = args.path
    run_path = '../runs/'
    if args.cplex:
        run_path += 'c_'
        executable = 'esso_cplex.o'
    else:
        run_path += 'h_'
        executable = 'esso_heuristic.o'
    run_path += args.id

    # check whether dataset folder exists
    if not os.path.isdir(dataset_path):
        print "error: dataset folder " + dataset_path + " not found"
        exit()

    # check whether run folder exists
    if not args.replace and os.path.isdir(run_path):
        print "error: run folder " + run_path + " already exists"
        exit()

    # create run directory
    try:
        os.mkdir(run_path)
    except OSError:
        pass

    # copy data files to the run folder
    shutil.copy(os.path.join(dataset_path, 'paths.dat'), run_path)
    #shutil.copyfile(os.path.join(dataset_path, 'init_topology.dat'),
    #                os.path.join(run_path,'res_topology_0.dat'))

    # copy the executable to the run folder
    if args.cplex:
        shutil.copy('esso_cplex.o', run_path)
    else:
        shutil.copy('esso_heuristic.o', run_path)

    # read input files
    read_vnf_types_file(dataset_path)
    #print vnf_flavor_to_cpu

    read_timeslots_file(dataset_path)

    read_topology_file(dataset_path)

    # main simulation loop
    # x_sfcs: sfcs that can be considered for migration at a time instance
    x_sfcs = set()
    # change path to the run folder
    cdir = os.getcwd()
    os.chdir(run_path)
    update_write_topology_file(0, None)
    for t in range(timeslot_count):
        x_sfcs = x_sfcs.difference(sfc_out[t])
        # x_sfcs represent alive sfcs that arrived between [0,t)
        # --- --- start code for simulation
        print "timeslot:", t, "n_sfc:", sfc_in[t], "x_sfc:", x_sfcs
        # write files for cplex/heuristic code
        n_sfc_filename = 'n_sfc_t' + str(t)
        with open(n_sfc_filename, 'w') as f:
            f.write(str(len(sfc_in[t])) + "\n")
            for s in sfc_in[t]:
                f.write(str(sfcs[s]) + "\n")
        x_sfc_filename = 'x_sfc_t' + str(t)
        with open(x_sfc_filename, 'w') as f:
            f.write(str(len(x_sfcs)) + "\n")
            for s in x_sfcs:
                f.write(str(sfcs[s]) + "\n")
        # invoke cplex/heuristic code if not dryrun
        if not args.dryrun:
            with open('run.log', 'w') as exe_log:
                topo_filename = 'res_topology_' + str(t) + '.dat'
                exe_path = './' + executable + ' ' + \
                        topo_filename + ' paths.dat'
                sys.stdout.flush()
                #start = timer()
                ts_sfcs = list(sfc_in[t])
                ts_sfcs.extend(list(x_sfcs))
                for s in ts_sfcs:
                    exe_proc = subprocess.Popen(exe_path, shell=True,
                            stdout=subprocess.PIPE,
                            stdin=subprocess.PIPE,
                            stderr=exe_log)
                    if s in sfc_in[t]:
                        mt = 0.0
                    else:
                        mt = migration_threshold
                    stdin_str = str(sfcs[s]) + ' ' + \
                            str(sfcs[s].cec) + ' ' + \
                            str(mt) + '\n'
                    #print stdin_str
                    exe_proc.stdin.write(stdin_str)
                    mapping = exe_proc.communicate()[0]
                    #print mapping.strip()
                    mapping_values = [int_or_float(x) for x in mapping.split()]
                    if mapping_values[0] == 200:
                        sfcs[s].cec = mapping_values[8+mapping_values[5]]
                    exe_proc.stdin.close()
                    print mapping.strip()
                    if exe_proc.wait() != 0:
                        print 'failed to execute code'
                        exit()
                    update_write_topology_file(t+1, mapping)
                #print "took", (timer()-start), "sec."
        # end of cplex/heuristic code execution
        x_sfcs = x_sfcs.union(sfc_in[t])
        # --- --- end code for simulation

    # change pwd to esso/src
    os.chdir(cdir)
