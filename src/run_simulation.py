#!/usr/bin/python

import sys
import os
import argparse
import subprocess
import shutil
from timeit import default_timer as timer
from collections import defaultdict
import json
import numpy as np
from scipy import stats

def int_or_float(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

class EssoSfc:
    def __init__(self):
        self.curr_emb_cost = -1.0 # current embedding cost
        pass

    # global settings
    id_idx = 0
    ttl_idx = 3

    # reads data for an SFC from a stream that 
    # provides a readline() function
    def read(self, strm):
        # save all sfc data as ints
        self.data = [int(x) for x in strm.readline().split()]
        # convert the vnf_types into their cpu requirements
        self.data[5:-2] = [int(vnf_flavor_to_cpu[int(x)]['cpu_count'])
                for x in self.data[5:-2]]

    def ttl(self):
        return self.data[EssoSfc.ttl_idx]

    def id(self):
        return self.data[EssoSfc.id_idx]

    def dec_ttl(self):
        self.data[EssoSfc.ttl_idx] -= 1
        if self.data[EssoSfc.ttl_idx] < 0:
            raise ValueError('Negative TTL')

    def __str__(self):
        return " ".join([str(x) for x in self.data])

class EssoObject:
    def __init__(self, data):
        self.data = data.strip()

    def __str__(self):
        return self.data

class EssoServer:
    # global settings
    cpu_idx = 5

    def __init__(self, data):
        self.data = data.split()
        self.data[EssoServer.cpu_idx] = int(
                self.data[EssoServer.cpu_idx])

    def dec_cpu_count(self, val):
        self.data[EssoServer.cpu_idx] -= val

    def inc_cpu_count(self, val):
        self.data[EssoServer.cpu_idx] += val

    def __str__(self):
        return " ".join([str(x) for x in self.data])

class EssoEdge:
    bw_idx = 5
    def __init__(self, data):
        self.data = data.split()
        self.data[EssoEdge.bw_idx] = int(self.data[EssoEdge.bw_idx])

    def dec_bandwidth(self, val):
        assert(self.data[EssoEdge.bw_idx] >= val);
        self.data[EssoEdge.bw_idx] -= val

    def inc_bandwidth(self, val):
        self.data[EssoEdge.bw_idx] += val

    def __str__(self):
        return " ".join(str(x) for x in self.data)

class SfcMapping:
    def __init__(self, data):
        data = [int_or_float(x) for x in data.split()]
        self.code = data[0]
        if self.code == 200:
            self.vnf_count = data[5]
            self.cpu_counts = data[6:6+self.vnf_count]
            self.bandwidth = data[6+self.vnf_count]
            self.emb_cost = data[8+self.vnf_count]
            self.emb_servers = data[10+self.vnf_count:10+2*self.vnf_count]
            paths_data = data[10+2*self.vnf_count:-1] 
            self.path_stretch = 0
            itr = iter(paths_data)
            path_count = itr.next()
            self.emb_paths = [[] for i in range(path_count)]
            for p in range(path_count):
                edge_count = itr.next()
                self.path_stretch += edge_count
                u = itr.next()
                for e in range(edge_count - 1):
                    v = itr.next()
                    self.emb_paths[p].append((u, v))
                    u = v
            self.run_time = data[-3]
            self.brown_energy = data[-2]
            self.green_energy = data[-1]
    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__)

vnf_flavor_to_cpu = {}
sfcs = []
sfc_mappings = []
sfc_in = [set()]
sfc_out = [set()]
sfc_count = 0
timeslot_count = 0
co_list = []
node_list = []
edge_list = []
#edge_dir = defaultdict(lambda: defaultdict(int))
edge_dir = {}
carbon_fp = 0
brown_energy = 0
green_energy = 0

def allocate_resource(sfc_id):
    global node_list
    global edge_list
    global edge_dir
    global carbon_fp
    global brown_energy
    global green_energy
    sfc_map = sfc_mappings[sfc_id]
    if sfc_map:
        for i in range(sfc_map.vnf_count):
            node_list[sfc_map.emb_servers[i]].dec_cpu_count(
                    sfc_map.cpu_counts[i])
        for path in sfc_map.emb_paths:
            for (u, v) in path:
                if u != v:
                    str_u, str_v = str(u), str(v)
                    edge_list[edge_dir[str_u][str_v]].dec_bandwidth(
                            sfc_map.bandwidth)
        carbon_fp += sfc_map.emb_cost
        brown_energy += sfc_map.brown_energy
        green_energy += sfc_map.green_energy
    #else:
        #print "ERROR: failed to find sfc", sfc_id, "mapping"

def release_resource(sfc_id):
    global node_list
    global edge_list
    global edge_dir
    global carbon_fp
    global brown_energy
    global green_energy
    sfc_map = sfc_mappings[sfc_id]
    if sfc_map:
        for i in range(sfc_map.vnf_count):
            node_list[sfc_map.emb_servers[i]].inc_cpu_count(
                    sfc_map.cpu_counts[i])
        for path in sfc_map.emb_paths:
            for (u, v) in path:
                if u != v:
                    str_u, str_v = str(u), str(v)
                    edge_list[edge_dir[str_u][str_v]].inc_bandwidth(
                            sfc_map.bandwidth)
        carbon_fp -= sfc_map.emb_cost
        brown_energy -= sfc_map.brown_energy
        green_energy -= sfc_map.green_energy
    #else:
        #print "ERROR: failed to find sfc", sfc_id, "mapping"

def update_write_topology_file():
    global co_list
    global node_list
    global edge_dir
    global edge_list

    #if mapping:
    #    m = SfcMapping(mapping)
    #    print m.toJSON()
    #    mv = [int_or_float(x) for x in mapping.split()]
    #    if mv[0] == 200:
    #        vc = mv[5]
    #        cpus = mv[6:6+vc]
    #        servers = mv[10+vc:10+2*vc]
    #        bw = mv[6+vc]
    #        paths_data = mv[10+2*vc:]
    #        for i in range(vc):
    #            node_list[servers[i]].dec_cpu_count(cpus[i])
    #        itr = iter(paths_data)
    #        path_count = itr.next()
    #        for p in range(path_count):
    #            edge_count = itr.next()
    #            u = itr.next()
    #            for e in range(edge_count-1):
    #                v = itr.next()
    #                print u, v
    #                if u != v:
    #                    str_u, str_v = str(u), str(v)
    #                    edge_list[edge_dir[str_u][str_v]].dec_bandwidth(bw)
    #                u = v
    #        ms_time = itr.next()
    #        print 'emb time: ', ms_time

    with open('res_topology.dat', 'w') as f:
        f.write(str(len(co_list)) + '\n')
        for co in co_list:
            f.write(str(co) + '\n')
        f.write(str(len(node_list)) + ' ' + str(len(edge_list)) + '\n')
        for node in node_list:
            f.write(str(node) + '\n')
        for edge in edge_list:
            f.write(str(edge) + '\n')

def read_topology_file(dataset_path):
    global co_list
    global node_list
    global edge_dir
    global edge_list
    with open(os.path.join(dataset_path, 'init_topology.dat')) as f:
        co_count = int(f.readline())
        for c in range(co_count):
            co = EssoObject(f.readline())
            co_list.append(co)
        node_count, edge_count = [int(x) for x in f.readline().split()]
        for n in range(node_count):
            line = f.readline()
            values = line.split()
            if values[1] == "c":
                node_list.append(EssoServer(line))
            else:
                node_list.append(EssoObject(line))
        for e in range(edge_count):
            line = f.readline()
            values = line.split()
            u, v = values[1:3]
            if u not in edge_dir:
                edge_dir[u] = {}
            edge_dir[u][v] = len(edge_list)
            edge_list.append(EssoEdge(line))

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
    global sfcs
    global sfc_mappings
    with open(os.path.join(dataset_path, "timeslots.dat")) as f:
        sfc_count = int(f.readline())
        sfc_mappings = [None] * sfc_count
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

                sfc = EssoSfc()
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
    except OSError as e:
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
    topo_filename = 'res_topology.dat'
    carbon_fp = 0
    embed_sfc_count = 0
    prced_sfc_count = 0
    running_times = []
    for t in range(timeslot_count):
        x_sfcs = x_sfcs.difference(sfc_out[t])
        # x_sfcs represent alive sfcs that arrived between [0,t)
        # --- --- start code for simulation
        #print "timeslot:", t, "n_sfc:", list(sfc_in[t]), "x_sfc:", list(x_sfcs)


        # to store the path stretch
        path_stretches = []
        # to count migrations
        migration_count = 0

        # release resource for expired sfcs
        for s in sfc_out[t]:
            release_resource(s)
        update_write_topology_file()
  
        print t, round(carbon_fp, 3), round(brown_energy, 3), \
                round(green_energy, 3),

        ########################################
        # write files for cplex/heuristic code
        #n_sfc_filename = 'n_sfc_t' + str(t)
        #with open(n_sfc_filename, 'w') as f:
        #    f.write(str(len(sfc_in[t])) + "\n")
        #    for s in sfc_in[t]:
        #        f.write(str(sfcs[s]) + "\n")
        #x_sfc_filename = 'x_sfc_t' + str(t)
        #with open(x_sfc_filename, 'w') as f:
        #    f.write(str(len(x_sfcs)) + "\n")
        #    for s in x_sfcs:
        #        f.write(str(sfcs[s]) + "\n")
        ########################################        
        # invoke cplex/heuristic code if not dryrun
        if not args.dryrun:
            with open('run.log', 'w') as exe_log:
                if args.cplex:
                    exe_path = './' + executable + ' ' + \
                            topo_filename + ' paths.dat'
                else:
                    exe_path = './' + executable + ' ' + \
                            '../' + os.path.join(dataset_path, 
                                'co_topology.dat') + \
                            ' ' + topo_filename
                #print 'run_sim: exe_path:', exe_path 
                sys.stdout.flush()
                #start = timer()
                ts_sfcs = list(sfc_in[t])
                ts_sfcs.extend(list(x_sfcs))
                for s in ts_sfcs:
                    if s in sfc_in[t]:
                        prced_sfc_count += 1
                    exe_proc = subprocess.Popen(exe_path, shell=True,
                            stdout=subprocess.PIPE,
                            stdin=subprocess.PIPE,
                            stderr=exe_log)
                    if s in sfc_in[t]:
                        mt = 0.0
                    else:
                        mt = migration_threshold
                    stdin_str = str(t) + ' ' + str(sfcs[s]) + ' ' + \
                            str(sfcs[s].curr_emb_cost) + ' ' + \
                            str(mt) + '\n'
                    #print 'input', stdin_str.strip()
                    exe_proc.stdin.write(stdin_str)
                    mapping = exe_proc.communicate()[0]
                    #print mapping.strip()
                    mapping_values = [int_or_float(x) for x in mapping.split()]
                    if mapping_values[0] == 200:
                        smp = SfcMapping(mapping)
                        sfcs[s].curr_emb_cost = smp.emb_cost
                        running_times.append(smp.run_time)
                        if s in sfc_in[t]: # new sfc
                            embed_sfc_count += 1
                            path_stretches.append(smp.path_stretch)
                        if s in x_sfcs: # migration
                            migration_count += 1
                            release_resource(s)
                        sfc_mappings[s] = smp
                        allocate_resource(s)
                    exe_proc.stdin.close()
                    #print mapping.strip()
                    if exe_proc.wait() != 0:
                        print 'failed to execute code'
                        exit()
                    update_write_topology_file()
                    #print embed_sfc_count, '/', prced_sfc_count, \
                    #        '/', sfc_count
                #print "took", (timer()-start), "sec."
        # end of cplex/heuristic code execution
        print embed_sfc_count*100.0/prced_sfc_count
        # print migration count and path stretch stat
        print t, migration_count, 
        if path_stretches:
            print min(path_stretches), np.percentile(path_stretches, 5), \
                    np.mean(path_stretches), \
                    np.percentile(path_stretches, 95), max(path_stretches)
        else:
            print "0 0 0 0 0"
        x_sfcs = x_sfcs.union(sfc_in[t])
        #if t == 1: break
        # --- --- end code for simulation

    # print stat data
    print min(running_times), np.percentile(running_times, 5), \
            np.mean(running_times), np.median(running_times), \
            stats.mode(running_times)[0][0], \
            np.percentile(running_times, 95), max(running_times)

    # change pwd to esso/src
    os.chdir(cdir)
