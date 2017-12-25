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

int main(int argc, char **argv) {

  int co_count, co_id, node_count, edge_count;
  vector<int> servers, switches;
  vector<node_info> node_infos;
  fstream fin("phy_inf_cplex.dat");
  fin >> co_count;
  vector<vector<double>> renewable_energy(co_count,
      vector<double>(time_inst_count));
  vector<double> carbon_per_watt(time_inst_count);
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
  izlib::iz_topology topo;
  topo.init(node_count);
  for (int i = 0; i < edge_count; ++i) {
    fin >> node_u >> node_v >> type >> bandwidth >> latency;
    topo.add_edge(node_u, node_v, latency, bandwidth);
  }

  fin.close();

  // physical paths
  constexpr int phy_k = 3;
  bool use_one_path = true;
  izlib::iz_path_list phy_paths;
  /*
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
  */

  //path to switch and path to edge mapping
  vector<vector<int>> path_to_switch(phy_paths.size(),
      vector<int>(node_count, 0));
  vector<vector<vector<int>>> path_to_edge(phy_paths.size(),
      vector<vector<int>>(node_count, vector<int>(node_count, 0)));
  /*
  for (int i = 0; i < phy_paths.size(); ++i) {
    for (auto& edge : topo.path_edges(phy_paths[i])) {
      path_to_edge[i][edge.u][edge.v] = 1;
      path_to_edge[i][edge.v][edge.u] = 1;
      if (node_infos[edge.u].node_category == 's') 
        path_to_switch[i][edge.u] = 1;
      if (node_infos[edge.v].node_category == 's') 
        path_to_switch[i][edge.v] = 1;
    }
  }
  */

  // read sfc data
  fstream ftin("timeslots.dat");
  int total_sfc_count{0}, time_instance_count{0};
  ftin >> total_sfc_count >> time_instance_count;

  vector<vector<int>> sfc_active(total_sfc_count, 
      vector<int>(time_instance_count, 0));
  vector<vector<int>> sfc_arrival(total_sfc_count, 
      vector<int>(time_instance_count, 0));
  vector<vector<int>> sfc_departure(total_sfc_count, 
      vector<int>(time_instance_count, 0));

  
  int sfc_count{0}, current_sfc{0};
  int ingress_co, egress_co, ttl, flavor_count, flavor_id, delay;
  for (int i = 0; i < time_instance_count; ++i) {
    ftin >> sfc_count;
    for (int j = 0; j < sfc_count; ++j) {
      ftin >> ingress_co >> egress_co >> ttl >> flavor_count;
      sfc_arrival[current_sfc][i] = 1;
      sfc_departure[current_sfc][i+ttl] = 1;
      for (int k = 0; k < ttl; ++k) {
        sfc_active[current_sfc][i+k] = 1;
      }
      for (int k = 0; k < flavor_count; ++k) {
        ftin >> flavor_id;
      }
      ftin >> bandwidth >> delay;
    }
  }

  ftin.close();

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
