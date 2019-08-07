#!/usr/bin/python

import logging
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

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

def int_or_float(s):
    """
    Converts `s` into either int or float and returns the value.
    First tries to convert `s` to int, if there is an ValueError
    exception, catches the error and tries to convert `s` to float.

    Args:
        s (string): A string representing numerical value

    Returns:
        int or float: The numerical value of the provided string

    Raises:
        ValueError: If fails to convert s into both int and float
    """
    try:
        return int(s)
    except ValueError:
        return float(s)


class EssoSfc:
    """
    EssoSfc represents an SFC in ESSO.
    This class is used to store the input SFC and
    the TTL value of an SFC during the simulation.
    """
    def __init__(self):
        self.curr_emb_cost = -1.0 # current embedding cost
        self.data = []
        pass

    # global settings
    # index of SFC-id and SFC-TTL in the input sting.split()
    id_idx = 0
    ttl_idx = 3

    # reads data for an SFC from a stream that 
    # supports a readline() function
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
    """
    EssoObject is a generic class just to store input data
    without any parsing (e.g., CO data)
    """
    def __init__(self, data):
        self.data = data.strip()

    def __str__(self):
        return self.data


class EssoServer:
    """
    Class to represent ESSO Server instances.
    Provides functions to increase and decrease available CPU counts
    during simulation.
    """
    # global settings
    # index for the CPU count in the input string
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
    """
    This class represents an ESSO edge in the simulation.
    Provides functions to increase and decrease the available
    bandwidth on an edge during simulation.
    """
    # global setting
    # index of bandwidth in the input string.split()
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
    """
    Class to store the SFC to infrastructure mapping.
    The constructor takes as input the output of an
    optimization algorithm and parses the data. It
    provided a function to retrieve the data as JSON.
    """
    def __init__(self, data):
        data = [int_or_float(x) for x in data.split()]
        self.code = data[0]
        if self.code == 200:
            self.vnf_count = data[5]
            self.cpu_counts = data[6:6+self.vnf_count]
            self.bandwidth = data[6+self.vnf_count]
            self.emb_cost = data[8+self.vnf_count]
            self.emb_servers = data[10+self.vnf_count:10+2*self.vnf_count]
            self.co_count = data[10+2*self.vnf_count]
            paths_data = data[11+2*self.vnf_count:-1]
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
            self.server_count = len(set(self.emb_servers))

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__)


# global data structures
# from vnf_types file
vnf_flavor_to_cpu = {} # stores a mapping between a VNF flavor and how many CPU it needs

# from timeslots file
sfcs = [] # list of all SFCs across timeslots, ordered according to their position in the input file
sfc_mappings = [] # list of `SfcMapping` instances represeting the mapping of SFC in `sfcs`
sfc_in = [set()] # a set to keep track of the inbound SFCs for a timeslot
sfc_out = [set()] # a set to keep track of the outbound SFCs for a timeslot
sfc_count = 0 # # number of total SFCs in the input file, counting across all timeslots
timeslot_count = 0 # number of timeslots in the input

# from co_topology/init_topology file
co_list = [] # list of COs in the topology
node_list = [] # list contains all nodes (switch and server)
edge_list = [] # all edges (inter and intra)
#edge_dir = defaultdict(lambda: defaultdict(int))
edge_dir = {} # a dict to store edge-id to `EssoEdge` instance mapping
carbon_fp = 0 # variable to keep track of carbon footprint during simulation
brown_energy = 0 # tracks brown energy during simulation
green_energy = 0 # tracks green energy during simulation


def allocate_resource(sfc_id):
    """
    This function allocates resources for an SFC.
    It assumes that the mapping for this SFC is already available through
    the list sfc_mappings. It reads the mapping from the list and allocates
    resources according to the mapping
    :param sfc_id:
    :return:
    """
    global node_list
    global edge_list
    global edge_dir
    global carbon_fp
    global brown_energy
    global green_energy
    sfc_map = sfc_mappings[sfc_id]
    if sfc_map:
        # allocate resource for VNF
        for i in range(sfc_map.vnf_count):
            node_list[sfc_map.emb_servers[i]].dec_cpu_count(
                    sfc_map.cpu_counts[i])
        # allocate resource for inter-VNF links on physical paths
        for path in sfc_map.emb_paths:
            for (u, v) in path:
                if u != v:
                    str_u, str_v = str(u), str(v)
                    edge_list[edge_dir[str_u][str_v]].dec_bandwidth(
                            sfc_map.bandwidth)
        # update the global variable related to
        # carbon footprint, brown and green energy
        carbon_fp += sfc_map.emb_cost
        brown_energy += sfc_map.brown_energy
        green_energy += sfc_map.green_energy
    #else:
        #print "ERROR: failed to find sfc", sfc_id, "mapping"


def release_resource(sfc_id):
    """
    This function releases the previously allocated resources for an SFC.
    There is no mechanism to ensure that resources for this SFC was actually
    allocated in an previous iteration. So, calling this function will either
    cause an assertion to fail somewhere in the code or give wrong results.
    TODO: add a check to ensure that resource was actually allocated for an SFC before freeing it.
    :param sfc_id:
    :return:
    """
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
    """
    This function updated the res_topology file based on the allocation
    and deallocation (release) of resources. The initial topology of the
    network is read from the init_topology file. After that the updated
    resource status is always writen and read from the res_topology file.
    :return:
    """
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

    # the class representation of objects greatly simplifies this code
    # the __str__ method of the EssoXYZ classes return a representation
    # that matches exactly with the format of the init_topology and
    # res_topology file hence the following code writes the updated resource
    # status of the network to the res_topology file.
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
    """
    Read the init_topology file, parses the data, and stores the information
    into EssoXYZ classes in different lists like
    co_list, node_list, edge_list and edge_dir
    :param dataset_path:
    :return:
    """
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
    """
    Reads in the vnf_type file
    :param dataset_path:
    :return:
    """
    global vnf_flavor_to_cpu
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
    """
    Reads the timeslots file
    :param dataset_path:
    :return:
    """
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
    """
    Main function for the simulation...
    """

    # parse command line arguments
    parser = argparse.ArgumentParser()
    # dataset path (must be relative to pwd) for the simulation
    parser.add_argument('dataset_path', help = "path to dataset folder (e.g., ../data/set0")

    # param group for selecting optimization algorithm
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--cplex', action='store_true',
            help="run CPLEX code")
    group.add_argument('-t', '--tabusearch', action='store_true',
            help="run tabu search code")
    group.add_argument('-f', '--firstfit', action='store_true',
                       help="run first fit code")

    # the `id` is used as the folder name for the current simulation
    # under the `runs` folder
    parser.add_argument('-i', '--id', required=True,
            help="unique id for this run")

    # if the `-r` option is provide then any existing directory under
    # the `runs` folder with the same run id will be replaced
    parser.add_argument('-r', '--replace', action='store_true',
            help="replace existing run dir")

    # the `-d` option is for generating the data only, no optimizer is run
    parser.add_argument('-d', '--dryrun', action='store_true',
            help="only generate data")

    # migration threshold determines the % reduction is cost to
    # trigger a migration
    parser.add_argument('-m', '--migthr', type=float, default=0.3,
            help='migration threshold (default=0.3)')
    args = parser.parse_args()

    # set the migration_threshold
    migration_threshold = args.migthr

    # path to the dataset and run folder
    dataset_path = args.dataset_path
    run_path = '../runs/'
    if args.cplex:
        run_path += 'c_'
        executable = 'esso_cplex.o'
    elif args.tabusearch:
        run_path += 't_'
        executable = 'esso_heuristic.o'
    else:
        run_path += 'f_'
        executable = 'esso_firstfit.o'
    run_path += args.id

    # check whether dataset folder exists
    if not os.path.isdir(dataset_path):
        logging.ERROR("dataset folder " + dataset_path + " not found")
        exit()

    # check whether run folder exists
    # if the 'replace' or '-r' option is provided then the following
    # check will not be done
    if not args.replace and os.path.isdir(run_path):
        logging.ERROR("run folder " + run_path + " already exists")
        exit()

    # create/replace run directory
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
    elif args.tabusearch:
        shutil.copy('esso_heuristic.o', run_path)
    else:
        shutil.copy('esso_firstfit.o', run_path)

    # read input files
    read_vnf_types_file(dataset_path)
    #print vnf_flavor_to_cpu

    # read the timeslots file
    read_timeslots_file(dataset_path)

    # read the init_topology file
    read_topology_file(dataset_path)

    ########################
    # main simulation loop #
    ########################
    # x_sfcs: sfcs that can be considered for migration at a time instance
    x_sfcs = set()
    # change path to the run folder, save the current pwd in cdir
    # will switch back to cdir at the end of the loop
    cdir = os.getcwd()
    os.chdir(run_path)
    # res_topology.dat is used to keep track of resources
    # during the simulation
    topo_filename = 'res_topology.dat'
    # to keep track of carbon footprint during the simulation
    carbon_fp = 0
    # counter to track how many SFCs were successfully embedded
    embed_sfc_count = 0
    # counter to track how many SFCs were attempted for embedding
    # the ratio of these two counters provide the acceptance ratio
    prced_sfc_count = 0
    # this list keeps track of the running times for each SFC
    running_times = []

    # open files for writing simulation output
    timeslot_data_file = open('timeslot_data.csv', 'w')
    sfc_data_file = open('sfc_data.csv', 'w')
    timeslot_data_file.write('timeslot,carbon_footprint,brown_energy,green_energy,' +
                        'acceptance_ratio,migration_count,' +
                        'ps_min,ps_5th,ps_mean,ps_95th,ps_max\n')
    timeslot_data_file.flush();
    sfc_data_file.write('timeslot,sfc_id,vnf_count,server_count,co_count,path_stretch\n')
    sfc_data_file.flush()
    # loop over the timeslots
    for t in range(timeslot_count):
        # remove the SFCs that are expiring at this timestamp
        # from the set x_sfcs
        # x_sfcs represent alive sfcs that arrived between [0,t)
        x_sfcs = x_sfcs.difference(sfc_out[t])
        logging.debug("timeslot: " + str(t) + " n_sfc: " + str(list(sfc_in[t])) + " x_sfc: " + str(list(x_sfcs)))

        #############################
        # start code for simulation #
        #############################

        # to store the path stretch
        path_stretches = []
        # to count migrations
        migration_count = 0

        # release resource for expired sfcs
        for s in sfc_out[t]:
            release_resource(s)
        # update the res_topology.dat file after releasing resources
        update_write_topology_file()

        # print <timeslot> <carbon-footprint> <brown-energy> <green-energy>
        # for each timeslot
        timeslot_data_file.write("{},{},{},{},".format(t, round(carbon_fp, 3),
                            round(brown_energy, 3), round(green_energy, 3)))

        # invoke cplex/heuristic code if not dryrun
        if not args.dryrun:
            # the `run.log` file within the `run` folder contains
            with open('run.log', 'a+') as exe_log:
                if args.cplex:
                    exe_path = './' + executable + ' ' + \
                            topo_filename + ' paths.dat'
                else:
                    exe_path = './' + executable + ' ' + \
                            '../' + os.path.join(dataset_path, 
                                'co_topology.dat') + \
                            ' ' + topo_filename
                logging.debug('run_sim: exe_path: %s', exe_path)
                # make sure that the stdout it written
                sys.stdout.flush()
                #start = timer()
                # ts_sfcs: sfcs (new and old) to be processed in this timestamp
                # first get the new ones
                ts_sfcs = list(sfc_in[t])
                # then get the old ones
                # comment out the following line to ignore migration
                ts_sfcs.extend(list(x_sfcs))
                # for each sfc s...
                for s in ts_sfcs:
                    # if s is a new sfc increment the prced_sfc_count
                    if s in sfc_in[t]:
                        prced_sfc_count += 1

                    # popen the cplex/heuristic program
                    exe_proc = subprocess.Popen(exe_path, shell=True,
                            stdout=subprocess.PIPE,
                            stdin=subprocess.PIPE,
                            stderr=exe_log)

                    # for new sfc migration threshold is zero
                    # otherwise we use the default/provided value during run
                    if s in sfc_in[t]:
                        mt = 0.0
                    else:
                        mt = migration_threshold

                    # build the stdin string to process the sfc
                    stdin_str = str(t) + ' ' + str(sfcs[s]) + ' ' + \
                            str(sfcs[s].curr_emb_cost) + ' ' + \
                            str(mt) + '\n'
                    logging.debug('input to optimizer: ' + stdin_str.strip())

                    # run the optimizer
                    exe_proc.stdin.write(stdin_str)
                    # get the output from the optimizer
                    mapping = exe_proc.communicate()[0]
                    logging.debug('output from optimizer: ' + mapping.strip())
                    # close the stdin of the exe_proc process
                    exe_proc.stdin.close()
                    # if failure to execute optimizer code then exit
                    if exe_proc.wait() != 0:
                        logging.error('failed to execute code')
                        exit()

                    # convert the string values in the mapping to int or float
                    mapping_values = [int_or_float(x) for x in mapping.split()]
                    # if the sfc was embedded, 200 indicates Okay.
                    if mapping_values[0] == 200:
                        smp = SfcMapping(mapping)
                        logging.debug('co_count: ' + str(smp.co_count))
                        sfcs[s].curr_emb_cost = smp.emb_cost
                        running_times.append(smp.run_time)
                        if s in sfc_in[t]: # new sfc
                            embed_sfc_count += 1
                            path_stretches.append(smp.path_stretch)
                            sfc_data_file.write('{},{},{},{},{},{}\n'.format(t, s,
                                                smp.vnf_count,
                                                smp.server_count,
                                                smp.co_count,
                                                smp.path_stretch))
                            #sfc_data_file.flush()
                        if s in x_sfcs: # migration
                            migration_count += 1
                            release_resource(s)
                        sfc_mappings[s] = smp
                        allocate_resource(s)
                    # update the res_topology.dat file after resource allocation
                    update_write_topology_file()

                    #print embed_sfc_count, '/', prced_sfc_count, \
                    #        '/', sfc_count
                #print "took", (timer()-start), "sec."
        # end of cplex/heuristic code execution
        # print the acceptance rate on the same line as carbon footprint,
        # brown & green energy
        timeslot_data_file.write('{},'.format(embed_sfc_count*100.0/prced_sfc_count))
        # print timestamp, migration count and path stretch stat
        timeslot_data_file.write('{},'.format(migration_count))
        if path_stretches:
            timeslot_data_file.write('{},{},{},{},{}'.format(min(path_stretches),
                                np.percentile(path_stretches, 5), np.mean(path_stretches),
                                np.percentile(path_stretches, 95), max(path_stretches)))
        else:
            timeslot_data_file.write('0.0,0.0,0.0,0.0,0.0')
        timeslot_data_file.write('\n')
        #timeslot_data_file.flush()

        # get the sfcs for the current timeslot and add them into x_sfcs
        # they will be considered for migration in the next timeslot
        x_sfcs = x_sfcs.union(sfc_in[t])

        # for testing purpose
        #if t == 1: break

        ###########################
        # end code for simulation #
        ###########################

    # print stat data
    timeslot_data_file.close()
    sfc_data_file.close()
    print min(running_times), np.percentile(running_times, 5), \
            np.mean(running_times), np.median(running_times), \
            stats.mode(running_times)[0][0], \
            np.percentile(running_times, 95), max(running_times)

    # change pwd to esso/src
    os.chdir(cdir)
