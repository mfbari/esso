#include <iostream>

#include "data_store.hpp"

using namespace std;

void data_store::read_input(char **argv) {
  // check args for correct format
  if (argc != 5) {
    cout << "usage: ./esso_cplex.o <path-to res_topology.dat> " <<
      "<path-to paths.dat> <path-to n_sfc_tn.dat> " <<
      "<path-to x_sfc_tn.dat>" << endl;
    exit(-1);
  }
  // filenames for reading in the inputs
  auto& res_topology_filename = argv[1];
  auto& paths_filename= argv[2];
  auto& n_sfc_filename = argv[3];
  auto& x_sfc_filename = argv[4];

  read_res_topology_data(res_topology_filename,
      co_count, node_count, edge_count,
      servers, switches,
      node_infos,
      renewable_energy,
      carbon_per_watt,
      topo);

  read_path_data(paths_filename,
      node_count, path_count,
      path_nodes, path_switches, edge_to_path);

  read_n_sfc_data(n_sfc_filename, n_sfcs);
  read_x_sfc_data(x_sfc_filename, x_sfcs);
}

void data_store::read_res_topology_data(const string& filename,
    int& co_count, int& node_count, int& edge_count,
    vector<int>& servers, vector<int>& switches,
    vector<node_info>& node_infos,
    vector<vector<double>>& renewable_energy,
    vector<double>& carbon_per_watt,
    izlib::iz_topology& topo) {
  fstream fin(filename);
  int co_id;
  fin >> co_count;
  renewable_energy.resize(co_count, vector<double>(time_inst_count));
  carbon_per_watt.resize(time_inst_count);
  // read renewable energy data from file
  for (int i = 0; i < co_count; ++i) {
    fin >> co_id;
    for (auto& re : renewable_energy[co_id]) {
      fin >> re;
    }
  }
  // fill the co2/kWh data, random double between [1,2]
  //TODO: use random generator
  for (auto& cw : carbon_per_watt) {
    cw = 1.34;
  }

  fin >> node_count >> edge_count;
  cout << node_count << " " << edge_count << endl;
  // read the nodes
  int node_id, cpu_count;
  char type;
  double sleep_power, base_power, per_cpu_power;
  for (int i = 0; i < node_count; ++i) {
    fin >> node_id >> type >> co_id;
    switch(type) {
      case 'c':
        servers.push_back(node_id);
        fin >> sleep_power >> base_power >> cpu_count >> per_cpu_power;
        node_infos.emplace_back(node_id, co_id, type,
            sleep_power, base_power,
            cpu_count, per_cpu_power);
        break;

      case 's':
        switches.push_back(node_id);
        fin >> sleep_power >> base_power;
        node_infos.emplace_back(node_id, co_id, type,
            sleep_power, base_power);
        break;

      default:
        cerr << "Error while parsing node list in phy_inf_clex.dat" << endl;

    }
  }
  // read the links
  int node_u, node_v, bandwidth, latency;
  topo.init(node_count);
  for (int i = 0; i < edge_count; ++i) {
    fin >> node_u >> node_v >> type >> bandwidth >> latency;
    topo.add_edge(node_u, node_v, latency, bandwidth);
  }
  fin.close();
}

void data_store::read_path_data(const string& paths_filename,
    const int node_count, int& path_count,
    vector<vector<int>>& path_nodes,
    vector<vector<int>>& path_switches,
    vector<vector<vector<int>>>& edge_to_path) {
  fstream fin(paths_filename);
  if (!fin) {
    cerr << "Failed to open paths.dat" << endl;
    exit(-1);
  }
  fin >> path_count;
  // resize the vectors
  path_nodes.resize(path_count, vector<int>());
  path_switches.resize(path_count, vector<int>());
  edge_to_path.resize(node_count, vector<vector<int>>(node_count,
        vector<int>()));
  // read data
  int id, n, u, v;
  for (int i = 0; i < path_count; ++i) {
    fin >> id >> n >> u;
    path_nodes[i].push_back(u);
    for (int j = 1; j < n; ++j) {
      fin >> v;
      path_nodes[i].push_back(v);
      edge_to_path[u][v].push_back(i);
      edge_to_path[v][u].push_back(i);
      u = v;
    }
    fin >> n; // switch-count
    path_switches[i].resize(n);
    for (auto& s : path_switches[i]) fin >> s;
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
  ofstream fin(x_sfc_filename);
  // TODO: how to represent mapping


  fin.close();
}
