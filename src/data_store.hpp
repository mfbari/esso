#ifndef DATA_STORE_HPP
#define DATA_STORE_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <map>

using namespace std;

// assuming that all simulations run for 24 hours
constexpr int time_inst_count = 24;

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
  bool is_server() {return node_category == 'c';}
};

// edge_info.t == 'i' for intra-co links or internal links to a co
// edge_info.t == 'b' for inter-co links or backbone links
struct edge_info {
  int id, u, v;
  char type;
  int capacity, latency;
  vector<int> paths; // path-ids that include this edge
  edge_info(int id = -1, int u = -1, int v = -1, 
      char t = 'i', int c = -1, int l = -1) :
    id(id), u(u), v(v), type(t), capacity(c), latency(l) {
  }
  bool is_valid() {
    return id != -1; 
  }
};

struct sfc_request {
  int id, vnf_count;
  int ingress_co, egress_co;
  int ttl;
  std::vector<int> vnf_flavors;
  std::vector<int> cpu_reqs;
  double latency;
  double bandwidth;
  int node_count() {return vnf_count + 2;}
  int edge_count() {return vnf_count + 1;}
  friend istream& operator>>(istream& is, sfc_request& sfc_req);
};
using sfc_request_set = std::vector<sfc_request>;

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

ostream& operator<<(ostream& out, const sfc_request& sfc_req) {
  out << sfc_req.id << " " << sfc_req.ingress_co << " " << 
    sfc_req.egress_co << " " << sfc_req.ttl << " " << 
    sfc_req.vnf_count << " ";
  for (const auto& v : sfc_req.cpu_reqs) {
    out << v << " ";
  }
  out << sfc_req.bandwidth << " " << sfc_req.latency;
  return out;
}
//=========the main data_store=========//

struct data_store {
  //=========variables to represent the inputs=========//
  // variables representing physical infrastucture
  int co_count, node_count, server_count, switch_count, 
      edge_count, path_count;
  // these two vectors contain the node_ids for the
  // servers and switches
  std::vector<int> servers, switches;
  std::unordered_map<int, int> server_id_to_index;
  // this vector contains all info for the nodes
  std::vector<node_info> node_infos;
  // vector<vector> for 24 renewable energy data for each co
  std::vector<std::vector<double>> renewable_energy;
  // carbon/watt at the COs
  std::vector<double> carbon_per_watt;
  // edges represents all edges in the network 
  // it is used to track delay and bandwidth usage and
  // convert path -> (u, v) to path -> edge-id
  vector<edge_info> edges;
  map<pair<int, int>, int> edge_uv_to_id;
  // path to node and switch & edge to path mappings
  // these are used for the constraints in the cplex model
  std::vector<std::vector<int>> path_nodes;
  std::vector<std::vector<int>> path_switches;
  std::vector<std::vector<int>> path_edge_ids;
  vector<vector<int>> edge_id_to_paths;
  std::unordered_map<int, vector<int>> switch_to_paths;
  // co to server, switch, edge mapping
  std::vector<std::vector<int>> co_server_ids;
  std::vector<std::vector<int>> co_switch_ids;
  std::vector<std::vector<int>> co_edge_ids;
  std::vector<int> backbone_edge_ids;
  // contain the new and pre-existing sfcs
  sfc_request_set n_sfcs, x_sfcs;

  //=========primary function to process input=========//
  void read_input(int argc, char **argv);
  // functions to read specific files
  void read_res_topology_data(const string& filename);
  void read_path_data(const string& filename);
    void read_n_sfc_data(const std::string& n_sfc_filename,
        sfc_request_set& n_sfcs);
    void read_x_sfc_data(const std::string& x_sfc_filename,
        sfc_request_set& x_sfcs);
    int path_latency(int _p);
};

void data_store::read_input(int argc, char **argv) {
  // check args for correct format
  if (argc != 3) {
    cout << "usage: ./esso_cplex.o <path-to res_topology.dat> " <<
      "<path-to paths.dat>" << endl;
    exit(-1);
  }
  // filenames for reading in the inputs
  auto& res_topology_filename = argv[1];
  auto& paths_filename= argv[2];
  //auto& n_sfc_filename = argv[3];
  //auto& x_sfc_filename = argv[4];

  read_res_topology_data(res_topology_filename);

  read_path_data(paths_filename);

  //read_n_sfc_data(n_sfc_filename, n_sfcs);
  //read_x_sfc_data(x_sfc_filename, x_sfcs);
}

void data_store::read_res_topology_data(const string& filename) {
  fstream fin(filename);
  int co_id;
  fin >> co_count;
  renewable_energy.resize(co_count, vector<double>(time_inst_count));
  carbon_per_watt.resize(co_count);
  // initialize the vectors for co to server, edge, switch mapping
  // no need to resize backbone_edge_ids, used .push_back function
  co_server_ids.resize(co_count);
  co_switch_ids.resize(co_count);
  co_edge_ids.resize(co_count);
  // read renewable energy data from file
  double carbon;
  for (int i = 0; i < co_count; ++i) {
    fin >> co_id;
    fin >> carbon;
    carbon_per_watt[i] = carbon;
    for (auto& re : renewable_energy[co_id]) {
      fin >> re;
    }
  }

  fin >> node_count >> edge_count;
  // read the nodes
  int node_id, cpu_count;
  char type;
  double sleep_power, base_power, per_cpu_power;
  for (int i = 0; i < node_count; ++i) {
    fin >> node_id >> type >> co_id;
    switch(type) {
      case 'c':
        server_id_to_index[node_id] = servers.size();
        servers.push_back(node_id);
        fin >> sleep_power >> base_power >> cpu_count >> per_cpu_power;
        node_infos.emplace_back(node_id, co_id, type,
            sleep_power, base_power,
            cpu_count, per_cpu_power);
        co_server_ids[co_id].push_back(node_id);
        break;

      case 's':
        switches.push_back(node_id);
        fin >> sleep_power >> base_power;
        node_infos.emplace_back(node_id, co_id, type,
            sleep_power, base_power);
        co_switch_ids[co_id].push_back(node_id);
        break;

      default:
        cerr << "Error while parsing node list in phy_inf_clex.dat" << endl;

    }
  }
  server_count = servers.size();
  switch_count = switches.size();
  // read the links
  int id, node_u, node_v, capacity, latency;
  for (int i = 0; i < edge_count; ++i) {
    fin >> id >> node_u >> node_v >> type >> co_id >> capacity >> latency;
    if (type == 'i') co_edge_ids[co_id].push_back(id);
    else backbone_edge_ids.push_back(id);
    edge_uv_to_id.insert(make_pair(make_pair(node_u, node_v), edges.size()));
    edges.emplace_back(id, node_u, node_v, type, capacity, latency);
  }
  fin.close();
}

void data_store::read_path_data(const string& filename) {
  fstream fin(filename);
  if (!fin) {
    cerr << "Failed to open paths.dat" << endl;
    exit(-1);
  }
  fin >> path_count;
  // resize the vectors
  path_nodes.resize(path_count, vector<int>());
  path_switches.resize(path_count, vector<int>());
  path_edge_ids.resize(path_count);
  // read data
  int id, n, u, v;
  for (int i = 0; i < path_count; ++i) {
    fin >> id >> n >> u;
    path_nodes[i].push_back(u);
    for (int j = 1; j < n; ++j) {
      fin >> v;
      path_nodes[i].push_back(v);
      auto edge_itr = edge_uv_to_id.find(make_pair(u, v)); 
      // this check is needed for skipping server loops
      // if they are added to the res_topology file then 
      // this must be removed
      if (edge_itr != edge_uv_to_id.end()) {
        auto& edge = edges[edge_itr->second];
        edge.paths.push_back(i);
        path_edge_ids[i].push_back(edge.id);
      }
      u = v;
    }
    fin >> n; // switch-count
    path_switches[i].resize(n);
    for (auto& s : path_switches[i]) {
      fin >> s;
      switch_to_paths[s].push_back(i);
    }
  }
  fin.close();
}

void data_store::read_n_sfc_data(const string& n_sfc_filename,
    sfc_request_set& n_sfcs) {
  fstream fin(n_sfc_filename);
  int id, n, cpu_count;
  fin >> n;
  for (int i = 0; i < n; ++i) {
    sfc_request sfc_req;
    fin >> id >> sfc_req.ingress_co >> sfc_req.egress_co >>
        sfc_req.ttl >> sfc_req.vnf_count;
    for (int k = 0; k < sfc_req.vnf_count; ++k) {
      fin >> cpu_count;
      sfc_req.cpu_reqs.push_back(cpu_count);
    }
    fin >> sfc_req.bandwidth >> sfc_req.latency;
    n_sfcs.push_back(sfc_req);
  }
  fin.close();
}

void data_store::read_x_sfc_data(const string& x_sfc_filename,
    sfc_request_set& x_sfcs) {
  fstream fin(x_sfc_filename);
  int id, n, cpu_count;
  fin >> n;
  for (int i = 0; i < n; ++i) {
    sfc_request sfc_req;
    fin >> id >> sfc_req.ingress_co >> sfc_req.egress_co >>
        sfc_req.ttl >> sfc_req.vnf_count;
    for (int k = 0; k < sfc_req.vnf_count; ++k) {
      fin >> cpu_count;
      sfc_req.cpu_reqs.push_back(cpu_count);
    }
    fin >> sfc_req.bandwidth >> sfc_req.latency;
    x_sfcs.push_back(sfc_req);
  }
  fin.close();
}

int data_store::path_latency(int _p) {
  int latency{0}, u{path_nodes[_p].front()}, v;
  for (size_t i = 1; i < path_nodes[_p].size(); ++i) {
    v = path_nodes[_p].at(i);
    if (u == v) continue;
    auto edge_id = edge_uv_to_id[make_pair(u, v)];
    latency += edges[edge_id].latency;
    u = v;
  }
  return latency;
}


#endif
