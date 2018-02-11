#ifndef DATA_STORE_HPP
#define DATA_STORE_HPP

#include <vector>
#include <string>

#include "iz_topology.hpp"

constexpr int time_inst_count = 24;

struct data_store {
  //=========variables to represent the inputs=========//
  // variables representing physical infrastucture
  int co_count, node_count, edge_count, path_count;
  // these two vectors contain the node_ids for the
  // servers and switches
  std::vector<int> servers, switches;
  // this vector contains all info for the nodes
  std::vector<node_info> node_infos;
  // vector<vector> for 24 renewable energy data for each co
  std::vector<std::vector<double>> renewable_energy;
  // carbon/watt at the COs
  std::vector<double> carbon_per_watt;
  // topo represents the entire topology of the
  // physical infrastructure. It is used to track
  // delay and bandwidth usage
  izlib::iz_topology topo;
  // path to node and switch & edge to path mappings
  // these are used for the constraints in the cplex model
  std::vector<std::vector<int>> path_nodes;
  std::vector<std::vector<int>> path_switches;
  std::vector<std::vector<std::vector<int>>> edge_to_path;
  // contain the new and pre-existing sfcs
  sfc_request_set n_sfcs, x_sfcs;


  //=========primary function to process input=========//
  void read_input(char **argv);
  // functions to read specific files
  void read_res_topology_data(const std::string& filename,
      int& co_count, int& node_count, int& edge_count,
      vector<int>& servers, vector<int>& switches,
      vector<node_info>& node_infos,
      vector<vector<double>>& renewable_energy,
      vector<double>& carbon_per_watt,
      izlib::iz_topology& topo);
  void read_path_data(const std::string& paths_filename,
      const int node_count, int& path_count,
      vector<vector<int>>& path_nodes,
      vector<vector<int>>& path_switches,
      vector<vector<vector<int>>>& edge_to_path);
    void read_n_sfc_data(const std::string& n_sfc_filename,
        sfc_request_set& n_sfcs);
    void read_x_sfc_data(const std::string& x_sfc_filename,
        sfc_request_set& x_sfcs);
};

//=========data structures to hold input data=========//

// node_info represents both server and switch
struct node_info {
  int node_id, co_id;
  char node_category;
  double sleep_power, base_power;
  int cpu_capacity;
  double per_cpu_power;
  node_info(int node_id, int co_id, char node_category,
      double sleep_power, double base_power,
      int cpu_capacity = 0, double per_cpu_power = 0.0) :
    node_id(node_id), co_id(co_id), node_category(node_category),
    sleep_power(sleep_power), base_power(base_power),
    cpu_capacity(cpu_capacity), per_cpu_power(per_cpu_power) {}
};

// not used yet
struct edge {
  int u, v;
};

struct sfc_request {
  int vnf_count;
  int ingress_co, egress_co;
  int ttl;
  std::vector<int> vnf_flavors;
  std::vector<int> cpu_reqs;
  double latency;
  double bandwidth;
  int node_count() {return vnf_count + 2;}
  int edge_count() {return vnf_count + 1;}
};
using sfc_request_set = vector<sfc_request>;


#endif
