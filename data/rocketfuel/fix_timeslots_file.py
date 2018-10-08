import sys

class SFC:
    def __init__(self, data):
        data = [int(x) for x in data.split()]
        self.id = data[0]
        self.ingress_co, self.egress_co = data[1:3]
        self.ttl = data[3]
        self.vnf_count = data[4]
        self.cpu_counts = data[5:5+self.vnf_count]
        self.bandwidth = data[5+self.vnf_count]
        self.max_delay = data[5+self.vnf_count+1]

if len(sys.argv) != 2:
    print 'input error'
    exit()

traffic_filename = sys.argv[1]

with open(traffic_filename, 'r') as f:
    sfc_count = int(f.readline())
    timeslot_count = int(f.readline())
    max_timepoint = 0
    timeslots = -1 # zero based indexing
    while True:
        line = f.readline()
        if line != '':
            timeslots += 1
            ts_sfc_count = int(line)
            #print 'timeslot:', timeslots, 'sfc_count:', ts_sfc_count,
            for s in range(ts_sfc_count):
                ttl = int(f.readline().split()[3])
                max_timepoint = max(max_timepoint, timeslots+ttl)
            #print 'max_timepoint:', max_timepoint
        else:
            break

timeslots += 1 # zero based indexing
max_timepoint += 1
print 'timeslots:', timeslots, 'max_timepoint:', max_timepoint

with open(traffic_filename, 'r') as f:
    file_data = f.readlines()

file_data[1] = str(max_timepoint) + '\n'
file_data.extend(['0\n'] * (max_timepoint - timeslots))

with open(traffic_filename, 'w') as f:
    f.writelines(file_data)
