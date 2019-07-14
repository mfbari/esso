#ifndef PROBLEM_INSTANCE_HPP
#define PROBLEM_INSTANCE_HPP

#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>

#include "esso_topology.hpp"

using namespace std;

struct vnf_type {
  int type_id;
  string type_name;
  vnf_type(int type_id, const string& type_name) :
    type_id(type_id), type_name(type_name) {}
};

struct vnf_flavor {
  int flavor_id;
  int type_id;
  int cpu_cores;
  double proc_delay;
  vnf_flavor(int flavor_id, int type_id, 
            int cpu_cores, double proc_delay) :
    flavor_id(flavor_id), type_id(type_id), 
    cpu_cores(cpu_cores), proc_delay(proc_delay) {}
};

struct sfc_request {
  int id, vnf_count;
  int ingress_co, egress_co;
  int ttl;
  vector<int> vnf_flavors;
  vector<int> cpu_reqs;
  double latency;
  double bandwidth;
  friend istream& operator>>(istream& is, sfc_request& sfc_req);
};

istream& operator>>(istream& is, sfc_request& sfc_req) {
  is >> sfc_req.id >> sfc_req.ingress_co >> sfc_req.egress_co >>
      sfc_req.ttl >> sfc_req.vnf_count;
  for (int k = 0, cpu_count; k < sfc_req.vnf_count; ++k) {
    is >> cpu_count;
    sfc_req.cpu_reqs.push_back(cpu_count);
  }
  is >> sfc_req.bandwidth >> sfc_req.latency;
  return is;
}

using sfc_request_set = std::vector<sfc_request>;
//typedef vector<vnf_request> vnf_request_set;

struct problem_input {
  string topology_filename;
  string vnf_info_filename;
  string time_slot_filename;
};

struct sfc_mapping {
  struct node_mapping {
    int co_id, server_id;
    node_mapping(int co_id, int server_id) :
        co_id(co_id), server_id(server_id) {}
  };
  enum class edge_mapping_type {
    in_server,
    intra_co,
    inter_co
  };
  // if n vnfs then (n+1) edge_mapping
  struct edge_mapping {
    edge_mapping_type ee_type;
    std::vector<int> path;
    edge_mapping(edge_mapping_type ee_type = 
        edge_mapping_type::in_server, std::vector<int> path = 
        std::vector<int> {}) {}
  };
  std::vector<node_mapping> node_mappings; // location where vnf is embedded
  std::vector<edge_mapping> edge_mappings;
};

using sfc_mapping_set = std::vector<sfc_mapping>;

struct problem_instance {
  vector<vnf_type> vnf_types;
  vector<vnf_flavor> vnf_flavors;
  vector<sfc_request_set> time_slots;
  vector<sfc_mapping_set> solution;
  esso_topology topology;

  int time_slot_count() {return time_slots.size();}

  bool read_input(const problem_input& prob_input) {
    return 
    //read_vnf_info_file(prob_input.vnf_info_filename) &&
    //read_time_slot_file(prob_input.time_slot_filename) &&
    read_topology_file(prob_input.topology_filename);

  }

  bool read_topology_file(const string& filename) {
    fstream fin(filename.c_str());
    if (!fin) {
      cout << "ERROR: failed to open topology file" << endl;
      return false;
    }
    fstream fgc("../data/greencap.dat");
    if (!fgc) {
        fgc.open("../../data/greencap.dat");
        if (fgc == NULL) {
            cout << "ERROR: failed to open greencap.dat" << endl;
            return false;
        }
    }
    int node_count{0}, edge_count{0};
    fin >> node_count >> edge_count;
    topology.init(node_count);
    for (int n = 0; n < node_count; ++n) {
      int node_id{0}, co_id, has_green;
      double carbon{0.0};
      vector<double> green_cap(24, 0.0);
      fin >> node_id >> carbon >> has_green ;
      //has_green = 0;
      if (has_green) {
        for (double& gc : green_cap) {
          fgc >> gc;
        }
      }
      co_id = topology.add_co(green_cap, carbon);
      assert(node_id == co_id);
    }
    for (int e = 0; e < 2*edge_count; ++e) {
      int u{0}, v{0}, capacity{0};
      double latency{0.0};
      fin >> u >> v >> latency >> capacity;
      // create edges
      topology.add_edge(u, v , latency, capacity);
    }
    fin.close();
    fgc.close();
    return true;
  }

  bool read_time_slot_file(const string& filename) {
    fstream fin(filename.c_str());
    if (!fin) {
      cout << "ERROR: failed to open timeslot file" << endl;
      return false;
    }
    int time_slot_count{0};
    fin >> time_slot_count;
    for (int t = 0; t < time_slot_count; ++t) {
      sfc_request_set sfc_req_set;
      int sfc_request_count{0};
      fin >> sfc_request_count;
      for (int r = 0; r < sfc_request_count; ++r) {
        sfc_request sfc_req;
        fin >> sfc_req.ingress_co >> sfc_req.egress_co >> sfc_req.ttl >> 
            sfc_req.vnf_count;
        for (int v = 0; v < sfc_req.vnf_count; ++v) {
          int vnf_flavor{0};
          fin >> vnf_flavor;
          sfc_req.vnf_flavors.push_back(vnf_flavor);
          sfc_req.cpu_reqs.push_back(vnf_flavors[vnf_flavor].cpu_cores);
        }
        fin >> sfc_req.bandwidth >> sfc_req.latency;
        sfc_req_set.emplace_back(sfc_req);
      }
      time_slots.emplace_back(sfc_req_set);
    }
    fin.close();
    return true;
  }    

  bool read_vnf_info_file(const string& filename) {
    fstream fin(filename.c_str());
    if (!fin) {
      cout << "ERROR: failed to open vnfinfo file" << endl;
      return false;
    }
    int type_count{0};
    fin >> type_count;
    for (int i = 0; i < type_count; ++i) {
      int type_id{0};
      string type_name{""};
      fin >> type_id >> type_name;
      vnf_types.emplace_back(vnf_type(type_id, type_name));
    }
    int flavor_count{0};
    fin >> flavor_count;
    for (int i = 0; i < flavor_count; ++i) {
      int flavor_id{0}, type_id{0}, cpus{0}, delay{0};
      fin >> flavor_id >> type_id >> cpus >> delay;
      vnf_flavors.emplace_back(vnf_flavor(flavor_id, type_id, cpus, delay));
    }
    fin.close();
    return true;
  }
};

ostream& operator<<(ostream& out, const sfc_request& sfc) {
  out << sfc.id << " " << sfc.ingress_co << " " << sfc.egress_co << " ";
  out << sfc.ttl << " " << sfc.vnf_count << " ";
  for (auto cpu : sfc.cpu_reqs) out << cpu << " ";
  out << sfc.bandwidth << " " << sfc.latency;
  return out;
}

#endif // PROBLEM_INSTANCE_HPP
