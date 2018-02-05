#include <ilcplex/ilocplex.h>
#include <vector>
#include <string>

#include "iz_topology.hpp"
#include "stop_watch.hpp"

ILOSTLBEGIN

constexpr int time_inst_count = 24;

//typedef IloArray<IloIntVarArray> IloIntVarArray2D;
using IloIntVarArray2D = IloArray<IloIntVarArray>;
using IloIntVarArray3D = IloArray<IloIntVarArray2D>;
using IloIntVarArray4D = IloArray<IloIntVarArray3D>;

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

void read_res_topology_data(const string& filename, 
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

  sfc_active.resize(total_sfc_count);
  sfc_arrival.resize(total_sfc_count);
  sfc_departure.resize(total_sfc_count);
  time_instances_sfcs.resize(time_instance_count);
  // read the timeslots.dat file for sfc requests and 
  // populate the active, arrival, departure events
  // and store the sfc requests in time_instances_sfcs 

  int sfc_count{0}, current_sfc{0};
  int flavor_id;
  for (int i = 0; i < time_instance_count; ++i) {
    fin >> sfc_count;
    sfc_active[i].resize(sfc_count, 0);
    sfc_arrival[i].resize(sfc_count, 0);
    sfc_departure[i].resize(sfc_count, 0);
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

void read_path_data(const string& path_switch_filename,
    const string& path_edge_filename,
    const int node_count, 
    vector<vector<int>>& path_to_switch,
    vector<vector<vector<int>>>& path_to_edge,
    vector<int>& path_src, vector<int>& path_dst) {
  fstream fin(path_switch_filename);
  int path_count;
  fin >> path_count;
  // resize the vectors
  path_to_switch.resize(path_count, vector<int>(node_count, 0));
  path_to_edge.resize(path_count, vector<vector<int>>(node_count, 
        vector<int>(node_count, 0)));
  path_src.resize(path_count, -1);
  path_dst.resize(path_count, -1);
  // read data
  int s, d, n, sw;
  for (int i = 0; i < path_count; ++i) {
    fin >> s >> d >> n;
    path_src[i] = s;
    path_dst[i] = d;
    for (int j = 0; j < n; ++j) {
      fin >> sw;
      path_to_switch[i][sw] = 1;
    }
  }
  fin.close();
  fin.open(path_edge_filename);
  fin >> path_count;
  for (int i = 0; i < path_count; ++i) {
    fin >> s >> d >> n;
    int u, v;
    fin >> u;
    for (int j = 0; j < n; ++j) {
      fin >> v;
      path_to_edge[i][u][v] = 1;
      path_to_edge[i][v][u] = 1;
      u = v;
    }
  }
  fin.close();
}

void read_n_sfc_data(const string& n_sfc_filename, 
    sfc_request_set& n_sfcs) {
  fstream fin(n_sfc_filename);
  int n, cpu_count;
  fin >> n;
  for (int i = 0; i < n; ++i) {
    sfc_request sfc_req;
    fin >> sfc_req.ingress_co >> sfc_req.egress_co >> sfc_req.ttl >> 
      sfc_req.vnf_count;
    for (int k = 0; k < sfc_req.vnf_count; ++k) {
      fin >> cpu_count;
      sfc_req.cpu_reqs.push_back(cpu_count);
      fin >> sfc_req.bandwidth >> sfc_req.latency;
    }
    n_sfcs.push_back(sfc_req);
  }
  fin.close();
}

void read_x_sfc_data(const string& x_sfc_filename, 
    sfc_request_set& x_sfcs) {
  ofstream fin(x_sfc_filename);
  // TODO: how to represent mapping


  fin.close();
}

int main(int argc, char **argv) {

  // check args for correnct format
  if (argc != 6) {
    cout << "usage: ./esso_cplex.o <path-to res_topology.dat> " <<
      "<path-to path_switch.dat> <path-to path_edge.dat> " <<
      "<path-to n_sfc_tn.dat> <path-to x_sfc_tn.dat>" << endl;
    exit(-1);
  }
  // filenames for reading in the inputs
  auto& res_topology_filename = argv[1];
  auto& path_switch_filename= argv[2];
  auto& path_edge_filename = argv[3];
  auto& n_sfc_filename = argv[4];
  auto& x_sfc_filename = argv[5];
  //string phy_inf_filename = "phy_inf_cplex.dat";
  //string vnf_info_filename = "vnfinfo.dat";
  //string time_instance_filename = "timeslots.dat";

  stop_watch sw;

  // variables representing physical infrastucture
  int co_count, node_count, edge_count;
  // these two vectors contain the node_ids for the 
  // servers and switches
  vector<int> servers, switches;
  // this vector contains all info for the nodes
  vector<node_info> node_infos;
  vector<vector<double>> renewable_energy;
  // carbon/watt at the COs
  vector<double> carbon_per_watt;
  // topo represents the entire topology of the 
  // physical infrastructure. It is used to compute
  // the k-shortest paths between the servers
  izlib::iz_topology topo;
  sw.start();
  read_res_topology_data(res_topology_filename, 
      co_count, node_count, edge_count,
      servers, switches,
      node_infos,
      renewable_energy, 
      carbon_per_watt,
      topo);
  sw.stop();
  cout << "res_topology.dat read in " << sw << endl;
  // compute the physical paths
  // k
  //constexpr int phy_k = 3;
  //bool use_one_path = true;
  // phy_paths holds all the paths, this is used for embedding
  //izlib::iz_path_list phy_paths;
  // path to switch and path to edge mapping
  // these are used for the constraints in the model
  vector<vector<int>> path_to_switch;
  vector<vector<vector<int>>> path_to_edge;
  vector<int> path_src, path_dst;
  sw.start();
  read_path_data(path_switch_filename, path_edge_filename,
      node_count, path_to_switch, path_to_edge,
      path_src, path_dst);
  sw.stop();
  cout << "path_* read in " << sw << endl;
  // read the vnfinfo.dat file for flavor_id to cpu_count
  // this is the flavor -> cpu mapping
  //vector<int> flavor_cpu;
  //read_vnf_info_data(vnf_info_filename, flavor_cpu);

  // read n_sfc_tn file
  sfc_request_set n_sfcs, x_sfcs;
  sw.start();
  read_n_sfc_data(n_sfc_filename, n_sfcs);
  read_x_sfc_data(x_sfc_filename, x_sfcs);
  sw.stop();
  cout << "n_sfc read in " << sw << endl;


  // read sfc data
  // sfcs are numbered across time instances
  //int total_sfc_count{0}, time_instance_count{0};
  // the following vectors represent active, arrival, and departure
  // event for the sfc across time instances
  //vector<vector<int>> sfc_active;
  //vector<vector<int>> sfc_arrival;
  //vector<vector<int>> sfc_departure;
  // mapping between time instance and sfc
  //vector<sfc_request_set> time_instances_sfcs;
  //read_time_instance_data(time_instance_filename, 
  //    flavor_cpu,
  //    total_sfc_count, time_instance_count,
  //    sfc_active, sfc_arrival, sfc_departure,
  //    time_instances_sfcs);

  IloEnv env;
  try {
    IloModel model(env);
    IloCplex cplex(model);

    cout << "---------cplex---------" << endl;
    /*
    // decision variable x, indices: t, i, n, _n
    IloIntVarArray3D x(env, time_instance_count);
    for (int t = 0; t < time_instance_count; ++t) {
      x[t] = IloIntVarArray2D(env, time_instances_sfcs[t].size());
      for (int i = 0; i < time_instances_sfcs[t].size(); ++i) {
        x[t][i] = IloIntVarArray2D(env, 
            time_instances_sfcs[t][i].vnf_count+2);
        for (int n = 0; n < time_instances_sfcs[t][i].vnf_count+2; ++n) {
          x[t][i][n] = IloIntVarArray(env, phy_paths.size(), 0, 1);
        }
      }
    }
    // decision variable y, indices: t, i, l, _p
    IloIntVarArray4D y(env, time_instance_count);
    for (int t = 0; t < time_instance_count; ++t) {
      y[t] = IloIntVarArray3D(env, time_instances_sfcs[t].size());
      for (int i = 0; i < time_instances_sfcs[t].size(); ++i) {
        y[t][i] = IloIntVarArray2D(env, 
            time_instances_sfcs[t][i].vnf_count+1);
        for (int l = 0; l < time_instances_sfcs[t][i].vnf_count+1; ++l) {
          y[t][i][l] = IloIntVarArray(env, servers.size(), 0, 1);
        }
      }
    }
    */
    // x and y should be equal to the active when summed over all 
    // servers and vnf nodes
    

  }
  catch (IloException& e) {
    cerr << "Concert expection: " << e << endl;
    e.end();
  }
  catch (...) {
    cerr << "Unknown Expection: " << endl;
  }

  env.end();
  cout << "---------cplex---------" << endl;
  return 0;
}
