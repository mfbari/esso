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

struct edge {
  int u, v;
};

struct sfc_request {
  int vnf_count;
  int ingress_co, egress_co;
  int ttl;
  vector<int> vnf_flavors;
  vector<int> cpu_reqs;
  double latency;
  double bandwidth;
  int node_count() {return vnf_count + 2;}
  int edge_count() {return vnf_count + 1;}
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

void read_path_data(const string& paths_filename,
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

void read_n_sfc_data(const string& n_sfc_filename, 
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

void read_x_sfc_data(const string& x_sfc_filename, 
    sfc_request_set& x_sfcs) {
  ofstream fin(x_sfc_filename);
  // TODO: how to represent mapping


  fin.close();
}

int main(int argc, char **argv) {

  // check args for correnct format
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
  int path_count {0};
  vector<vector<int>> path_nodes;
  vector<vector<int>> path_switches;
  vector<vector<vector<int>>> edge_to_path;
  sw.start();
  read_path_data(paths_filename, 
      node_count, path_count,
      path_nodes, path_switches, edge_to_path);
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
    // decision variable x, indices: i, n, _n
    IloIntVarArray3D x(env, n_sfcs.size());
    for (int i = 0; i < n_sfcs.size(); ++i) {
      x[i] = IloIntVarArray2D(env, n_sfcs[i].node_count());
      for (int n = 0; n < n_sfcs[i].node_count(); ++n) {
        x[i][n] = IloIntVarArray(env, servers.size(), 0, 1);
      }
    } 
    // decision variable y, indices: i, l, _p
    IloIntVarArray3D y(env, n_sfcs.size());
    for (int i = 0; i < n_sfcs.size(); ++i) {
      y[i] = IloIntVarArray2D(env, n_sfcs[i].edge_count());
      for (int l = 0; l < n_sfcs[i].edge_count(); ++l) {
        y[i][l] = IloIntVarArray(env, path_count, 0, 1);
      }
    } 
    // derived vairable z
    IloIntVarArray z(env, servers.size(), 0, 1);

    //=========Constraint=========//
    // x and y should be equal to 1 when summed over all 
    // servers and physical paths
    for (int i = 0; i < n_sfcs.size(); ++i) {
      for (int n = 0; n < n_sfcs[i].node_count(); ++n) {
        IloExpr x__n(env);
        for (int _n = 0; _n < servers.size(); ++_n) {
          x__n += x[i][n][_n];
        }
        model.add(x__n == 1);
      }
    }
    for (int i = 0; i < n_sfcs.size(); ++i) {
      for (int l = 0; l < n_sfcs[i].edge_count(); ++l) {
        IloExpr y__p(env);
        for (int _p = 0; _p < path_count; ++_p) {
          y__p += y[i][l][_p];
        }
        model.add(y__p == 1);
      }
    }

    //=========Constraint=========//
    // z[_n] whether a server _n is active or not
    for (int _n = 0; _n < servers.size(); ++_n) {
      IloExpr sum(env);
      for (int i = 0; i < n_sfcs.size(); ++i) {
        for (int n = 0; n < n_sfcs[i].node_count(); ++n) {
          sum += x[i][n][_n];
        }
      }
      model.add(z[_n] <= sum);
      model.add(IloIfThen(env, sum > 0, z[_n] == 1));
    }

    //=========Constraint=========//
    // capacity constraint for server
    for (int _n = 0; _n < servers.size(); ++_n) {
      IloExpr allocated_cpu(env);
      for (int i = 0; i < n_sfcs.size(); ++i) {
        // n = 0 is the ingress co
        // n = n_sfcs[i].node_count()-1 is the egress co
        for (int n = 1; n < n_sfcs[i].node_count() - 1; ++n) {
          allocated_cpu += x[i][n][_n] * n_sfcs[i].cpu_reqs[n-1];
        }
      }
      model.add(allocated_cpu <= node_infos[servers[_n]].cpu_capacity);
    }
    // capacity constraint for edge
    for (int u = 0; u < node_count; ++u) {
      for (int v = 0; v < node_count; ++v) {
        if (u != v) {
          IloExpr allocated_capacity(env);
          for (int i = 0; i < n_sfcs.size(); ++i) {
            for (int l = 0; l < n_sfcs[i].edge_count(); ++l) {
              for (int _p: edge_to_path[u][v]) {
                allocated_capacity += y[i][l][_p] * n_sfcs[i].bandwidth;
              }
            }
          }
          model.add(allocated_capacity <= topo.residual(u, v));
        }
      }
    }
               
    //=========Constraint=========//
    // max delay constraint for sfc
    for (int i = 0; i < n_sfcs.size(); ++i) {
      IloExpr delay(env);
      for (int l = 0; l < n_sfcs[i].edge_count(); ++l) {
        for (int _p = 0; _p < path_count; ++_p) {
          int u = path_nodes[_p].front(), v;
          for (int j = 1; j < path_nodes[_p].size(); ++j) {
            v = path_nodes[_p][j];
            delay += y[i][l][_p] * topo.latency(u, v);
            u = v;
          }
        }
      }
      model.add(delay <= n_sfcs[i].latency);
    }
    
    //=========Objective=========//
    IloExpr objective(env);
    // minimize the number of active servers
    IloExpr cost(env);
    for (int _n = 0; _n < servers.size(); ++_n) {
      cost += z[_n];
    }

    objective += cost;

    /*Objective --> model*/
    model.add(objective >= 0);
    model.add(IloMinimize(env, objective));

    //=========Solver=========//
    IloTimer timer(env);
    timer.restart();

    const IloInt time_limit_seconds = 60*60; // 1 hour
    const IloNum relative_gap = 0.001; // 0.1% gap with optimal

    cplex.setParam(IloCplex::TiLim, time_limit_seconds);
    cplex.setParam(IloCplex::EpGap, relative_gap);
    cplex.setParam(IloCplex::PreDual, true);

    if(!cplex.solve()) {
      timer.stop();
      cout << "could not solve ILP!" << endl;
      cout << "status: " << cplex.getStatus() << endl;
      throw(-1);
    }

    timer.stop();
    cout << "ILP solved in " << timer.getTime() << endl;
    cout << "Objective value = " << cplex.getObjValue() << endl;

    // get value for x
    for (int i = 0; i < n_sfcs.size(); ++i) {
      for (int n = 0; n < n_sfcs[i].node_count(); ++n) {
        cout << "x[" << i << "][" << n << "] = ";
        for (int _n = 0; _n < servers.size(); ++_n) {
          if (cplex.getValue(x[i][n][_n] == 1)) {
            cout << _n << " ";
          }
        }
        cout << endl;
      }
    }
    cout << "---------" << endl;
    // get value for y
    for (int i = 0; i < n_sfcs.size(); ++i) {
      for (int l = 0; l < n_sfcs[i].edge_count(); ++l) {
        cout << "y[" << i << "][" << l << "] = ";
        for (int _p = 0; _p < path_count; ++_p) {
          if (cplex.getValue(y[i][l][_p] == 1)) {
            cout << _p << " ";
          }
        }
        cout << endl;
      }
    }  


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
