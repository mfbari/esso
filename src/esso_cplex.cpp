#include <ilcplex/ilocplex.h>
#include <iterator>

#include "data_store.hpp"

ILOSTLBEGIN

using IloIntVarArray2D = IloArray<IloIntVarArray>;
using IloIntVarArray3D = IloArray<IloIntVarArray2D>;
using IloIntVarArray4D = IloArray<IloIntVarArray3D>;

int main(int argc, char **argv) {

  data_store ds;
  ds.read_input(argc, argv);

  // define lambda function to access source and destination
  auto _s = [&ds](int _p) {return ds.path_nodes[_p].front();};
  auto _d = [&ds](int _p) {return ds.path_nodes[_p].back();};
  auto s = [](int l) {return l;};
  auto d = [](int l) {return l+1;};

  IloEnv env;
  try {
    IloModel model(env);
    IloCplex cplex(model);

    cout << "---------cplex---------" << endl;
    // decision variable x, indices: i, n, _n
    IloIntVarArray3D x(env, ds.n_sfcs.size());
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      x[i] = IloIntVarArray2D(env, ds.n_sfcs[i].node_count());
      for (int n = 0; n < ds.n_sfcs[i].node_count(); ++n) {
        x[i][n] = IloIntVarArray(env, ds.node_count, 0, 1);
      }
    }
    // decision variable y, indices: i, l, _p
    IloIntVarArray3D y(env, ds.n_sfcs.size());
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      y[i] = IloIntVarArray2D(env, ds.n_sfcs[i].edge_count());
      for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
        y[i][l] = IloIntVarArray(env, ds.path_count, 0, 1);
      }
    }
    // derived vairable z
    IloIntVarArray z(env, ds.node_count, 0, 1);
    // derived varaible q
    IloIntVarArray q(env, ds.edge_count, 0, 1);
    // derived variable w
    IloIntVarArray w(env, ds.switches.size(), 0, 1);

    //=========Constraint=========//
    // x and y should be equal to 1 when summed over all
    // servers and physical paths
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int n = 0; n < ds.n_sfcs[i].node_count(); ++n) {
        IloExpr x__n(env);
        for (int _n = 0; _n < ds.node_count; ++_n) {
          if (ds.node_infos[_n].is_server()) {
            x__n += x[i][n][_n];
          }
          else {
            model.add(x[i][n][_n] == 0);
          }
        }
        model.add(x__n == 1);
      }
    }
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
        IloExpr y__p(env);
        for (int _p = 0; _p < ds.path_count; ++_p) {
          y__p += y[i][l][_p];
        }
        model.add(y__p == 1);
      }
    }

    //=========Constraint=========//
    // z[_n] whether a server _n is active or not
    for (int _n = 0; _n < ds.node_count; ++_n) {
      IloExpr sum(env);
      for (int i = 0; i < ds.n_sfcs.size(); ++i) {
        for (int n = 0; n < ds.n_sfcs[i].node_count(); ++n) {
          sum += x[i][n][_n];
        }
      }
      model.add(z[_n] <= sum);
      model.add(IloIfThen(env, sum > 0, z[_n] == 1));
    }
    // q[_l] whether a physical link is active or not
    int _l{0}, u, v;
    for (auto& p : ds.edges) {
      tie(u, v) = p;
      IloExpr sum(env);
      for (int _p : ds.edge_to_path[u][v]) {
        for (int i = 0; i < ds.n_sfcs.size(); ++i) {
          for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
            sum += y[i][l][_p];
          }
        }
      }
      //model.add(q[_l] <= sum);
      model.add(IloIfThen(env, sum > 0, q[_l] == 1));
      ++_l;
    }
    // w[_s] whether a switch is active or not
    int si{0};
    for (auto& sw_paths : ds.switch_to_path) {
      IloExpr sum(env);
      for (int _p : sw_paths.second) {
        for (int i = 0; i < ds.n_sfcs.size(); ++i) {
          for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
            sum += y[i][l][_p];
          }
        }
      }
      model.add(IloIfThen(env, sum > 0, w[si] == 1));
      ++si;
    }

    //=========Constraint=========//
    // capacity constraint for server
    for (int _n = 0; _n < ds.node_count; ++_n) {
      if (ds.node_infos[_n].is_server()) {
        IloExpr allocated_cpu(env);
        for (int i = 0; i < ds.n_sfcs.size(); ++i) {
          // n = 0 is the ingress co
          // n = n_sfcs[i].node_count()-1 is the egress co
          for (int n = 1; n < ds.n_sfcs[i].node_count() - 1; ++n) {
            allocated_cpu += x[i][n][_n] * ds.n_sfcs[i].cpu_reqs[n-1];
          }
        }
        model.add(allocated_cpu <= ds.node_infos[_n].cpu_capacity);
      }
    }
    // capacity constraint for edge
    for (int u = 0; u < ds.node_count; ++u) {
      for (int v = 0; v < ds.node_count; ++v) {
        if (u != v) {
          IloExpr allocated_capacity(env);
          for (int i = 0; i < ds.n_sfcs.size(); ++i) {
            for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
              for (int _p: ds.edge_to_path[u][v]) {
                allocated_capacity += y[i][l][_p] * ds.n_sfcs[i].bandwidth;
              }
            }
          }
          if (ds.topo[u][v].is_valid()) {
            model.add(allocated_capacity <= ds.topo[u][v].capacity);
          }
        }
      }
    }

    //=========Constraint=========//
    // max delay constraint for sfc
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      IloExpr delay(env);
      for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
        for (int _p = 0; _p < ds.path_count; ++_p) {
          // skip intra-server self-loops
          if (ds.path_nodes[_p].size() == 2) continue; 
          // actual paths
          int u = ds.path_nodes[_p].front(), v;
          for (int j = 1; j < ds.path_nodes[_p].size(); ++j) {
            v = ds.path_nodes[_p][j];
            if (ds.topo[u][v].is_valid()){
              delay += y[i][l][_p] * ds.topo[u][v].latency;
            }
            else {
              cout << "invalid edge is path" << endl;
              exit(-1);
            }
            u = v;
          }
        }
      }
      model.add(delay <= ds.n_sfcs[i].latency);
    }

    //=========Constraint=========//
    // flow conservation 
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
        for (int _p = 0; _p < ds.path_count; ++_p) {
          model.add(y[i][l][_p] <= x[i][s(l)][_s(_p)]);
          model.add(y[i][l][_p] <= x[i][d(l)][_d(_p)]);
        }
      }
    }

    //=========Constraint=========//
    // placement constraint for CO
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int _n = 0; _n < ds.node_count; ++_n) {
        if (ds.node_infos[_n].is_server()) {
          if (ds.node_infos[_n].co_id == ds.n_sfcs[i].ingress_co) {
            model.add(x[i][0][_n] <= 1);
          }
          else {
            model.add(x[i][0][_n] == 0);
          }
          if (ds.node_infos[_n].co_id == ds.n_sfcs[i].egress_co) {
            model.add(x[i][ds.n_sfcs[i].node_count()-1][_n] <= 1);
          }
          else {
            model.add(x[i][ds.n_sfcs[i].node_count()-1][_n] == 0);
          }
        }
      }
    }

    

    //=========Objective=========//
    IloExpr objective(env);
    // minimize the number of active servers
    IloExpr cost(env);
    for (int _n = 0; _n < ds.node_count; ++_n) {
      cost += z[_n];
    }
    for (int _l = 0; _l < ds.edge_count; ++_l) {
      cost += q[_l];
    }
    for (int _s; _s < ds.switches.size(); ++_s) {
      cost += w[_s];
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
    cout << "ILP solved in " << timer.getTime() << " sec" << endl;
    cout << "Objective value = " << cplex.getObjValue() << endl;

    // get value for x
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int n = 0; n < ds.n_sfcs[i].node_count(); ++n) {
        cout << "x[" << i << "][" << n << "] = ";
        for (int _n = 0; _n < ds.node_count; ++_n) {
          if (cplex.getValue(x[i][n][_n] == 1)) {
            cout << ds.node_infos[_n].co_id << ":" << _n << " ";
          }
        }
        cout << endl;
      }
      cout << "---" << endl;
    }
    // get value for y
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
        cout << "y[" << i << "][" << l << "] = ";
        for (int _p = 0; _p < ds.path_count; ++_p) {
          if (cplex.getValue(y[i][l][_p] == 1)) {
            cout << _p << ": ";
            copy(ds.path_nodes[_p].begin(), ds.path_nodes[_p].end(),
                ostream_iterator<int>(cout, " "));
            cout << "l:" << ds.path_latency(_p) << " ";
          }
        }
        cout << endl;
      }
      cout << "---" << endl;
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
