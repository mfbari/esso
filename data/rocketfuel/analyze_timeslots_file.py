
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



as_no = input('AS number: ')

traffic_filename = str(as_no) + '.timeslots'

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

sfcs = []
sfc_in = [[] for t in range(max_timepoint)]
sfc_out = [[] for t in range(max_timepoint)]

min_delay = 999999999
max_bandwidth = 0
with open(traffic_filename, 'r') as f:
    sfc_count = int(f.readline())
    timeslot_count = int(f.readline())
    for t in range(timeslots):
        ts_sfc_count = int(f.readline())
        for s in range(ts_sfc_count):
            sfc = SFC(f.readline())
            sfcs.append(sfc)
            sfc_in[t].append(sfc.id)
            sfc_out[t+sfc.ttl].append(sfc.id)
            min_delay = min(sfc.max_delay, min_delay)
            max_bandwidth = max(sfc.bandwidth, max_bandwidth)

print 'min delay', min_delay, 'max_bandwidth', max_bandwidth

sfc_count = 0
cpu_count = 0
bandwidth = 0
for t in range(max_timepoint-1):
    sfc_count += (len(sfc_in[t]) - len(sfc_out[t]))
    cpu_count += ( sum([sum(sfcs[s].cpu_counts) for s in sfc_in[t]]) -
                   sum([sum(sfcs[s].cpu_counts) for s in sfc_out[t]]) )
    bandwidth += ( sum([sfcs[s].bandwidth for s in sfc_in[t]]) - 
                   sum([sfcs[s].bandwidth for s in sfc_out[t]]) )
    print t+1, len(sfc_in[t]), len(sfc_out[t]), sfc_count, cpu_count, bandwidth

