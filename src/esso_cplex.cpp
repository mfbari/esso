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

    sfc_request sfc;
    double current_cost;
    // migrate if there is migration_threshold * 100 % 
    // cost reduction due to migration
    double migration_threshold; 
    cin >> sfc >> current_cost >> migration_threshold;

    // decision variable x, indices: n, _n
    IloIntVarArray2D x(env, sfc.node_count());
    for (int n = 0; n < sfc.node_count(); ++n) {
      x[n] = IloIntVarArray(env, ds.node_count, 0, 1);
    }
    // decision variable y, indices: l, _p
    IloIntVarArray2D y(env, sfc.edge_count());
    for (int l = 0; l < sfc.edge_count(); ++l) {
      y[l] = IloIntVarArray(env, ds.path_count, 0, 1);
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
    for (int n = 0; n < sfc.node_count(); ++n) {
      IloExpr x__n(env);
      for (int _n = 0; _n < ds.node_count; ++_n) {
        if (ds.node_infos[_n].is_server()) {
          x__n += x[n][_n];
        }
        else {
          model.add(x[n][_n] == 0);
        }
      }
      model.add(x__n == 1);
    }
    for (int l = 0; l < sfc.edge_count(); ++l) {
      IloExpr y__p(env);
      for (int _p = 0; _p < ds.path_count; ++_p) {
        y__p += y[l][_p];
      }
      model.add(y__p == 1);
    }

    //=========Constraint=========//
    // z[_n] whether a server _n is active or not
    for (int _n = 0; _n < ds.node_count; ++_n) {
      IloExpr sum(env);
      for (int n = 1; n < sfc.node_count() - 1; ++n) {
        sum += x[n][_n];
      }
      model.add(z[_n] <= sum);
      model.add(IloIfThen(env, sum > 0, z[_n] == 1));
    }
    // q[_l] whether a physical link is active or not
    for (int _l = 0; _l < ds.edge_count; ++_l) {
      IloExpr sum(env);
      for (int _p : ds.edges[_l].paths) {
        for (int l = 0; l < sfc.edge_count(); ++l) {
          sum += y[l][_p];
        }
      }
      //model.add(q[_l] <= sum);
      model.add(IloIfThen(env, sum > 0, q[_l] == 1));
    }
    // w[_s] whether a switch is active or not
    int si{0};
    for (auto& sw_paths : ds.switch_to_paths) {
      IloExpr sum(env);
      for (int _p : sw_paths.second) {
        for (int l = 0; l < sfc.edge_count(); ++l) {
          sum += y[l][_p];
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
        // n = 0 is the ingress co
        // n = sfc.node_count()-1 is the egress co
        for (int n = 1; n < sfc.node_count() - 1; ++n) {
          allocated_cpu += x[n][_n] * sfc.cpu_reqs[n-1];
        }
        model.add(allocated_cpu <= ds.node_infos[_n].cpu_capacity);
      }
    }
    // capacity constraint for edge
    for (auto& edge : ds.edges) {
      IloExpr allocated_capacity(env);
      for (int l = 0; l < sfc.edge_count(); ++l) {
        for (int _p: edge.paths) {
          allocated_capacity += y[l][_p] * sfc.bandwidth;
        }
      }
      model.add(allocated_capacity <= edge.capacity);
    }

    //=========Constraint=========//
    // max delay constraint for sfc
    IloExpr delay(env);
    for (int l = 0; l < sfc.edge_count(); ++l) {
      for (int _p = 0; _p < ds.path_count; ++_p) {
        // skip intra-server self-loops
        if (ds.path_nodes[_p].size() == 2) continue; 
        // actual paths
        for (int _l : ds.path_edge_ids[_p]) {
            delay += y[l][_p] * ds.edges[_l].latency;
        }
      }
    }
    model.add(delay <= sfc.latency);

    //=========Constraint=========//
    // flow conservation 
    for (int l = 0; l < sfc.edge_count(); ++l) {
      for (int _p = 0; _p < ds.path_count; ++_p) {
        model.add(y[l][_p] <= x[s(l)][_s(_p)]);
        model.add(y[l][_p] <= x[d(l)][_d(_p)]);
      }
    }

    //=========Constraint=========//
    // placement constraint for CO
    for (int _n = 0; _n < ds.node_count; ++_n) {
      if (ds.node_infos[_n].is_server()) {
        if (ds.node_infos[_n].co_id == sfc.ingress_co) {
          model.add(x[0][_n] <= 1);
        }
        else {
          model.add(x[0][_n] == 0);
        }
        if (ds.node_infos[_n].co_id == sfc.egress_co) {
          model.add(x[sfc.node_count()-1][_n] <= 1);
        }
        else {
          model.add(x[sfc.node_count()-1][_n] == 0);
        }
      }
    }

    //=========Constraint=========//
    // minimize the number of active servers, switches, edges
    IloExpr cost(env);
    for (int _n = 0; _n < ds.node_count; ++_n) {
      cost += z[_n];
    }
    for (int _l = 0; _l < ds.edge_count; ++_l) {
      cost += q[_l];
    }
    for (int _s; _s < ds.switch_count; ++_s) {
      cost += w[_s];
    }

    // this constraint is only added for pre-existing sfcs
    // a negative current_cost indicate new sfc
    // zero and positive current_cost indicate pre-existing sfc
    if (current_cost >= 0) {
      model.add(cost <= (1 - migration_threshold) * current_cost);
    }


    //=========Objective=========//
    IloExpr objective(env);
    // server power
    for (int _n = 0; _n < ds.node_count; ++_n) {
      if (ds.node_infos[_n].is_server()) {
        objective += z[_n] * ds.node_infos[_n].base_power + 
          (1 - z[_n]) * ds.node_infos[_n].sleep_power;
      }
    }
    IloExpr cpu_power(env);
    for (int _n = 0; _n < ds.node_count; ++_n) {
      if (ds.node_infos[_n].is_server()) {
        for (int n = 1; n < sfc.node_count() - 1; ++n) {
          cpu_power += x[n][_n] * (sfc.cpu_reqs[n-1] * 
              ds.node_infos[_n].per_cpu_power);
        }
      }
    }
    objective += cpu_power;
    // switch power
    for (int _s; _s < ds.switches.size(); ++_s) {
      objective += w[_s] * ds.node_infos[ds.switches[_s]].base_power +
        (1 - w[_s]) * ds.node_infos[ds.switches[_s]].sleep_power;
    }
    // edge power
    for (auto& edge : ds.edges) {
      IloExpr allocated_capacity(env);
      for (int l = 0; l < sfc.edge_count(); ++l) {
        for (int _p: edge.paths) {
          if (ds.path_nodes[_p].size() > 2) { // skip server self-loops 
            allocated_capacity += y[l][_p] * sfc.bandwidth;
          }
        }
      }
      objective += allocated_capacity * 0.006;
      /*
      objective += 0.001 * IloPiecewiseLinear(allocated_capacity, 
        0.0,
        IloNumArray(env, 10, 
          0.1, 0.1,
          1.1, 1.1, 
          10.1 ,10.1,
          100.1, 100.1,
          1000.1, 1000.1), 
        IloNumArray(env, 10, 
          0.0, 1.0,
          1.0, 3.0,
          3.0, 6.0, 
          6.0, 10.0,
          10.0, 18.0),
        0.0);
      */
    }


    /*Objective --> model*/
    model.add(objective >= 0);
    model.add(IloMinimize(env, objective));

    //=========Solver=========//
    IloTimer timer(env);
    timer.restart();

    const IloInt time_limit_seconds = 60*60; // 1 hour
    const IloNum relative_gap = 0.001; // 0.1% gap with optimal

    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::TiLim, time_limit_seconds);
    cplex.setParam(IloCplex::EpGap, relative_gap);
    cplex.setParam(IloCplex::PreDual, true);

    if(!cplex.solve()) {
      timer.stop();
      if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        cout << "404 " << sfc << endl;
      }
      else { 
        cerr << "could not solve ILP!" << endl;
        cerr << "status: " << cplex.getStatus() << endl;
      }
      throw(-1);
    }

    timer.stop();
    //cout << "ILP solved in " << timer.getTime() << " sec" << endl;
    //cout << "Objective value = " << cplex.getObjValue() << endl;
    cout << "200 " << sfc << " " << cplex.getObjValue() << " ";

    //for (int _n = 0; _n < ds.node_count; ++_n) {
    //  cout << _n << ":" << cplex.getValue(z[_n]) << " ";
    //}
    //cout << endl;
    // get value for x
    cout << sfc.vnf_count << " ";
    for (int n = 1; n < sfc.node_count() - 1; ++n) {
      //cout << "x[" << n << "] = ";
      for (int _n = 0; _n < ds.node_count; ++_n) {
        if (cplex.getValue(x[n][_n]) == 1) {
          //cout << ds.node_infos[_n].co_id << ":" << _n << " ";
          cout << _n << " ";
        }
      }
      //cout << endl;
    }
    //cout << "---" << endl;
    // get value for y
    vector<int> alledges;
    for (int l = 0; l < sfc.edge_count(); ++l) {
      //cout << "y[" << l << "] = ";
      for (int _p = 0; _p < ds.path_count; ++_p) {
        if (cplex.getValue(y[l][_p]) == 1) {
          //cout << _p << ": ";
          copy(ds.path_nodes[_p].begin(), ds.path_nodes[_p].end(),
              back_inserter(alledges));
          //cout << "l:" << ds.path_latency(_p) << " ";
        }
      }
      //cout << endl;
    }
    cout << alledges.size() - 1 << " ";
    copy(alledges.begin(), alledges.end(), 
        ostream_iterator<int>(cout, " "));
    //cout << "---" << endl;
    cout << endl;
    /*
    for (auto& edge : ds.edges) {
      cout << edge.u << "->" << edge.v << " ";
      IloExpr cap(env);
      for (int l = 0; l < sfc.edge_count(); ++l) {
        for (int _p: edge.paths) {
          cap += sfc.bandwidth * y[l][_p];
        }
      }
      cout << cplex.getValue(cap) << " ";
      cout << cplex.getValue(IloPiecewiseLinear(cap, 
        0.0,
        IloNumArray(env, 10, 
          0.1, 0.1,
          1.1, 1.1, 
          10.1 ,10.1,
          100.1, 100.1,
          1000.1, 1000.1), 
        IloNumArray(env, 10, 
          0.0, 1.0,
          1.0, 3.0,
          3.0, 6.0, 
          6.0, 10.0,
          10.0, 18.0),
        0.0)) << endl;
    }
    */
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
