#!/usr/bin/python

import sys
import os

vnf_flavor_to_cpu = {}
sfcs = []
sfc_in_events = [[]]
sfc_out_events = [[]]

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
    global sfc_in_events
    global sfc_out_events
    with open(os.path.join(dataset_path, "timeslots.dat")) as f:
        sfc_count = int(f.readline())
        timeslot_count = int(f.readline())
        sfc_in_events = [[] for t in range(timeslot_count)]
        sfc_out_events = [[] for t in range(timeslot_count)]
        sfc_id = 0
        for t in range(timeslot_count):
            print 'timeslot', t
            t_sfc_count = int(f.readline())
            for s in range(t_sfc_count):
                print 'sfc', s
                sfc_in_events[t].append(sfc_id)
                values = [x for x in f.readline().split()]  
                ttl = int(values[2])
                sfc_out_events[t+ttl].append(sfc_id)
                sfcs.append(values[:2] + values[3:])
                sfc_id += 1

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "usage: ./run_simulation.py or python run_simulation.py",
        print "<dataset-folder-name>"
        print "name of the folder not the full path"
        exit()

    dataset_path = os.path.join("data", sys.argv[1])  

    if not os.path.isdir(dataset_path):
        print "dataset folder " + dataset_path + " not found"
        exit()

    #print "\n".join(os.listdir(dataset_path))

    read_vnf_types_file(dataset_path)
    print vnf_flavor_to_cpu
    
    read_timeslots_file(dataset_path)
    print sfcs
    print sfc_in_events
    print sfc_out_events

