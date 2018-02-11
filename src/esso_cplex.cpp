#include <ilcplex/ilocplex.h>

#include "data_store.hpp"

ILOSTLBEGIN

using IloIntVarArray2D = IloArray<IloIntVarArray>;
using IloIntVarArray3D = IloArray<IloIntVarArray2D>;
using IloIntVarArray4D = IloArray<IloIntVarArray3D>;

int main(int argc, char **argv) {

  data_store ds;
  ds.read_input(argv);

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
        x[i][n] = IloIntVarArray(env, ds.servers.size(), 0, 1);
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
    IloIntVarArray z(env, ds.servers.size(), 0, 1);

    //=========Constraint=========//
    // x and y should be equal to 1 when summed over all
    // servers and physical paths
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int n = 0; n < ds.n_sfcs[i].node_count(); ++n) {
        IloExpr x__n(env);
        for (int _n = 0; _n < ds.servers.size(); ++_n) {
          x__n += x[i][n][_n];
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
    for (int _n = 0; _n < ds.servers.size(); ++_n) {
      IloExpr sum(env);
      for (int i = 0; i < ds.n_sfcs.size(); ++i) {
        for (int n = 0; n < ds.n_sfcs[i].node_count(); ++n) {
          sum += x[i][n][_n];
        }
      }
      model.add(z[_n] <= sum);
      model.add(IloIfThen(env, sum > 0, z[_n] == 1));
    }

    //=========Constraint=========//
    // capacity constraint for server
    for (int _n = 0; _n < ds.servers.size(); ++_n) {
      IloExpr allocated_cpu(env);
      for (int i = 0; i < ds.n_sfcs.size(); ++i) {
        // n = 0 is the ingress co
        // n = n_sfcs[i].node_count()-1 is the egress co
        for (int n = 1; n < ds.n_sfcs[i].node_count() - 1; ++n) {
          allocated_cpu += x[i][n][_n] * ds.n_sfcs[i].cpu_reqs[n-1];
        }
      }
      model.add(allocated_cpu <= ds.node_infos[ds.servers[_n]].cpu_capacity);
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
          model.add(allocated_capacity <= ds.topo.residual(u, v));
        }
      }
    }

    //=========Constraint=========//
    // max delay constraint for sfc
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      IloExpr delay(env);
      for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
        for (int _p = 0; _p < ds.path_count; ++_p) {
          int u = ds.path_nodes[_p].front(), v;
          for (int j = 1; j < ds.path_nodes[_p].size(); ++j) {
            v = ds.path_nodes[_p][j];
            delay += y[i][l][_p] * ds.topo.latency(u, v);
            u = v;
          }
        }
      }
      model.add(delay <= ds.n_sfcs[i].latency);
    }

    //=========Objective=========//
    IloExpr objective(env);
    // minimize the number of active servers
    IloExpr cost(env);
    for (int _n = 0; _n < ds.servers.size(); ++_n) {
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
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int n = 0; n < ds.n_sfcs[i].node_count(); ++n) {
        cout << "x[" << i << "][" << n << "] = ";
        for (int _n = 0; _n < ds.servers.size(); ++_n) {
          if (cplex.getValue(x[i][n][_n] == 1)) {
            cout << _n << " ";
          }
        }
        cout << endl;
      }
    }
    cout << "---------" << endl;
    // get value for y
    for (int i = 0; i < ds.n_sfcs.size(); ++i) {
      for (int l = 0; l < ds.n_sfcs[i].edge_count(); ++l) {
        cout << "y[" << i << "][" << l << "] = ";
        for (int _p = 0; _p < ds.path_count; ++_p) {
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
