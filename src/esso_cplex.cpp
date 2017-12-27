#include <ilcplex/ilocplex.h>
#include <vector>
#include <string>

#include "iz_topology.hpp"

ILOSTLBEGIN

constexpr int time_inst_count = 24;

//typedef IloArray<IloIntVarArray> IloIntVarArray2D;
using IloIntVarArray2D = IloArray<IloIntVarArray>;
using IloIntVarArray3D = IloArray<IloIntVarArray2D>;
using IloIntVarArray4d = IloArray<IloIntVarArray3D>;

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

struct sfc_request {
  int vnf_count;
  int ingress_co, egress_co;
  int ttl;
  vector<int> vnf_flavors;
  vector<int> cpu_reqs;
  double latency;
  double bandwidth;
};
using sfc_request_set = vector<sfc_request>;

void read_physical_infrastucture_data(const string& filename, 
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
            per_cpu_power, cpu_count);
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

void read_vnf_info_data(const string& filename, vector<int>& flavor_cpu) {
  fstream fin("vnfinfo.dat");
  int type_count, flavor_count, dummy_int, cpu_count;
  string dummy_str;
  fin >> type_count;
  for (int k = 0; k < type_count; ++k) fin >> dummy_int >> dummy_str;
  fin >> flavor_count;
  for(int i = 0; i < flavor_count; ++i) {
    fin >> dummy_int >> dummy_int >> cpu_count >> dummy_int;
    flavor_cpu.push_back(cpu_count);
  }
  fin.close();
}

void read_time_instance_data(const string& filename, 
    const vector<int>& flavor_cpu,
    int& total_sfc_count, int& time_instance_count, 
    vector<vector<int>>& sfc_active, 
    vector<vector<int>>& sfc_arrival,
    vector<vector<int>>& sfc_departure,
    vector<sfc_request_set>& time_instances_sfcs) {
  // read sfc data
  fstream fin(filename);
  fin >> total_sfc_count >> time_instance_count;

  sfc_active.resize(total_sfc_count, vector<int>(time_instance_count, 0));
  sfc_arrival.resize(total_sfc_count, vector<int>(time_instance_count, 0));
  sfc_departure.resize(total_sfc_count, vector<int>(time_instance_count, 0));
  time_instances_sfcs.resize(time_instance_count);
  // read the timeslots.dat file for sfc requests and 
  // populate the active, arrival, departure events
  // and store the sfc requests in time_instances_sfcs 

  int sfc_count{0}, current_sfc{0};
  int flavor_id;
  for (int i = 0; i < time_instance_count; ++i) {
    fin >> sfc_count;
    time_instances_sfcs[i].resize(sfc_count);
    for (int j = 0; j < sfc_count; ++j) {
      sfc_request sfc_req;
      fin >> sfc_req.ingress_co >> sfc_req.egress_co >> sfc_req.ttl >> 
        sfc_req.vnf_count;
      sfc_arrival[current_sfc][i] = 1;
      sfc_departure[current_sfc][i+sfc_req.ttl] = 1;
      for (int k = 0; k < sfc_req.ttl; ++k) {
        sfc_active[current_sfc][i+k] = 1;
      }
      for (int k = 0; k < sfc_req.vnf_count; ++k) {
        fin >> flavor_id;
        sfc_req.vnf_flavors.push_back(flavor_id);
        sfc_req.cpu_reqs.push_back(flavor_cpu[flavor_id]);
      }
      fin >> sfc_req.bandwidth >> sfc_req.latency;
    }
  }

  fin.close();
}

void calculate_physical_paths(const int& node_count, 
    const vector<node_info>& node_infos, 
    const int& phy_k, const bool& use_one_path,
    izlib::iz_topology& topo, izlib::iz_path_list& phy_paths,
    vector<vector<int>>& path_to_switch, 
    vector<vector<vector<int>>>& path_to_edge) {
  for (int u = 0; u < node_count; ++u) {
    for (int v = 0; v < node_count; ++v) {
      if (u != v && node_infos[u].node_category == 'c' &&
          node_infos[v].node_category == 'c') {
        // both u and v are servers, and u != v
        if (use_one_path) {
          izlib::iz_path path;
          topo.shortest_path(u, v, path);
          phy_paths.push_back(path);
        }
        else {
          izlib::iz_path_list paths;
          topo.k_shortest_paths(u, v, phy_k, paths); 
          phy_paths.insert(phy_paths.end(), paths.begin(), paths.end());
        }
      }
    }
    cout << u << " done." << '\r' << flush;
  }
  cout << endl;
  cout << "total paths: " << phy_paths.size() << endl;
  path_to_switch.resize(phy_paths.size(), vector<int>(node_count, 0));
  path_to_edge.resize(phy_paths.size(), vector<vector<int>>(node_count, 
      vector<int>(node_count, 0)));
  for (size_t i = 0; i < phy_paths.size(); ++i) {
    for (auto& edge : topo.path_edges(phy_paths[i])) {
      path_to_edge[i][edge.u][edge.v] = 1;
      path_to_edge[i][edge.v][edge.u] = 1;
      if (node_infos[edge.u].node_category == 's') 
        path_to_switch[i][edge.u] = 1;
      if (node_infos[edge.v].node_category == 's') 
        path_to_switch[i][edge.v] = 1;
    }
  }
}

int main(int argc, char **argv) {

  // filenames for reading in input
  string phy_inf_filename = "phy_inf_cplex.dat";
  string vnf_info_filename = "vnfinfo.dat";
  string time_instance_filename = "timeslots.dat";

  // variables representing physical infrastucture
  int co_count, node_count, edge_count;
  vector<int> servers, switches;
  vector<node_info> node_infos;
  vector<vector<double>> renewable_energy;
  vector<double> carbon_per_watt;
  izlib::iz_topology topo;
  read_physical_infrastucture_data(phy_inf_filename, 
      co_count, node_count, edge_count,
      servers, switches,
      node_infos,
      renewable_energy, 
      carbon_per_watt,
      topo);
  // physical paths
  constexpr int phy_k = 3;
  bool use_one_path = true;
  izlib::iz_path_list phy_paths;
  //path to switch and path to edge mapping
  vector<vector<int>> path_to_switch;
  vector<vector<vector<int>>> path_to_edge;
  calculate_physical_paths(node_count, 
      node_infos,
      phy_k, use_one_path,
      topo, phy_paths,
      path_to_switch, path_to_edge);

  // read the vnfinfo.dat file for flavor_id to cpu_count
  vector<int> flavor_cpu;
  read_vnf_info_data(vnf_info_filename, flavor_cpu);

  // read sfc data
  int total_sfc_count{0}, time_instance_count{0};
  vector<vector<int>> sfc_active;
  vector<vector<int>> sfc_arrival;
  vector<vector<int>> sfc_departure;
  vector<sfc_request_set> time_instances_sfcs;
  read_time_instance_data(time_instance_filename, 
      flavor_cpu,
      total_sfc_count, time_instance_count,
      sfc_active, sfc_arrival, sfc_departure,
      time_instances_sfcs);

  IloEnv env;
  try {
    IloModel model(env);

    cout << "cplex code ..." << endl;

  }
  catch (IloException& e) {
    cerr << "Concert expection: " << e << endl;
    e.end();
  }
  catch (...) {
    cerr << "Unknown Expection: " << endl;
  }

  env.end();
  return 0;
}
