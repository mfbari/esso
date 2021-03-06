#include <ilcplex/ilocplex.h>
#include <iterator>
#include <set>

#include "data_store.hpp"
#include "iz_timer.hpp"

ILOSTLBEGIN

using IloIntVarArray2D = IloArray<IloIntVarArray>;
using IloIntVarArray3D = IloArray<IloIntVarArray2D>;
using IloIntVarArray4D = IloArray<IloIntVarArray3D>;

// if this flag is set to false then, the switch cost will not be
// part of the objective function. However, the switch cost will
// be calculated and added to the final cost output
constexpr bool use_switch_cost_in_objective = true;

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
    int timeslot;
    double current_cost;
    // migrate if there is migration_threshold * 100 % 
    // cost reduction due to migration
    double migration_threshold; 
    cin >> timeslot >> sfc >> current_cost >> migration_threshold;

    // the following code block is for no consolidation case
//    if (current_cost > 0) {
//      cout << "404 " << sfc << endl;
//      return 0;
//    }

    iz_timer ctimer;

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
    // derived variables cap_x_y
    IloIntVarArray cap_0_1g(env, ds.edge_count, 0, 1);
    IloIntVarArray cap_1g_inf(env, ds.edge_count, 0, 1);
    // derived variable w
    //IloIntVarArray w(env, ds.switches.size(), 0, 1);
    IloIntVarArray w(env, ds.node_count, 0, 1);

    //=========Constraint=========//
    // placement constraint for ingress and egress co
    auto co_gateway = [](int co_id) {return co_id * 9;};
    for (int _n = 0; _n < ds.node_count; ++_n) {
      if (_n == co_gateway(sfc.ingress_co))
        model.add(x[0][_n] == 1);
      else 
        model.add(x[0][_n] == 0);
      if (_n == co_gateway(sfc.egress_co)) 
        model.add(x[sfc.node_count()-1][_n] == 1);
      else
        model.add(x[sfc.node_count()-1][_n] == 0);
    }

    //=========Constraint=========//
    // x and y should be equal to 1 when summed over all
    // servers and physical paths
    for (int n = 1; n < sfc.node_count()-1; ++n) {
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
    /*
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
    */
    for (int _s = 0; _s < ds.node_count; ++_s) {
      // if not server, then its a switch
      if (!ds.node_infos[_s].is_server()) {
        auto sw_paths = ds.switch_to_paths.find(_s);
        if (sw_paths == ds.switch_to_paths.end()) {
          continue;
        }
        IloExpr sum(env);
        for (int _p : sw_paths->second) {
          for (int l = 0; l < sfc.edge_count(); ++l) {
            sum += y[l][_p];
          }
        }
        model.add(IloIfThen(env, sum > 0, w[_s] == 1));
      }
      else {
        model.add(w[_s] == 0); // for server w[_s] is zero
      }
    }

    //=========Constraint=========//
    // capacity constraint for server
    for (int _n = 0; _n < ds.node_count; ++_n) {
      if (ds.node_infos[_n].is_server()) {
        IloExpr allocated_cpu(env);
        // n = 0 is the ingress co
        // n = sfc.node_count()-1 is the egress co
        // they are just placeholders and assumed to consume zero resource
        for (int n = 1; n < sfc.node_count() - 1; ++n) {
          allocated_cpu += x[n][_n] * sfc.cpu_reqs[n-1];
        }
        model.add(allocated_cpu <= ds.node_infos[_n].cpu_capacity);
      }
    }
    // capacity constraint for edge
    for (int _l = 0; _l < ds.edge_count; ++_l) {
      IloExpr allocated_capacity(env);
      for (int l = 0; l < sfc.edge_count(); ++l) {
        for (int _p: ds.edges[_l].paths) {
          allocated_capacity += y[l][_p] * sfc.bandwidth;
        }
      }
      model.add(allocated_capacity <= ds.edges[_l].capacity);
      // the following constraints are to track allocated
      // bandwidth on edges. They are used to assing energy cost
      // according to allocated bandwidth
      model.add(IloIfThen(env, allocated_capacity > 0 && 
            allocated_capacity <= 1000, cap_0_1g[_l] == 1));
      model.add(IloIfThen(env, allocated_capacity > 1000, 
            cap_1g_inf[_l] == 1));
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
    /*
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
    */

    //=========Constraint=========//
    // minimize the number of active servers, switches, edges
    /*
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
    */

    // this constraint is only added for pre-existing sfcs
    // a negative current_cost indicate new sfc
    // zero cost is not considered for migration
    // positive current_cost indicate pre-existing sfc and
    // is considered for migration
//    if (current_cost > 0) {
//      model.add(cost <= (1 - migration_threshold) * current_cost);
//    }


    //=========Objective=========//
    IloExpr objective(env);
    // server power
    // -----------------------------------------------------
    /*
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
    for (int _l = 0; _l < ds.edge_count; ++_l) {
      objective += 0.0012 * cap_0_1g[_l] +
                   0.0043 * cap_1g_inf[_l]; 
    }
    */
    // -----------------------------------------------------
    // The above codes for server, switch, and edge power are 
    // commented out to calculate power on a per co basis, 
    // as done in the following code block

    IloExprArray co_powers(env, ds.co_count);
    IloExprArray switch_powers(env, ds.co_count);
    IloExpr backbone_power(env);
    IloExprArray co_brown_power(env, ds.co_count);
    IloExpr total_brown_power(env), total_green_power(env);
    for (int _c = 0; _c < ds.co_count; ++_c) {
      //IloExpr co_power(env);
      co_powers[_c] = IloExpr(env);
      switch_powers[_c] = IloExpr(env);
      co_brown_power[_c] = IloExpr(env);
      // 1. server power on a per co basis
      for (int _n : ds.co_server_ids[_c]) {
        // server base/sleep power
        co_powers[_c] += z[_n] * ds.node_infos[_n].base_power + 
          (1 - z[_n]) * ds.node_infos[_n].sleep_power;
        // server cpu power
        for (int n = 1; n < sfc.node_count() - 1; ++n) {
          co_powers[_c] += x[n][_n] * (sfc.cpu_reqs[n-1] *
              ds.node_infos[_n].per_cpu_power);
        }
      }
      // 2. switch power on a per co basis
      for (int _s : ds.co_switch_ids[_c]) {
        switch_powers[_c] += w[_s] * ds.node_infos[_s].base_power +
            (1 - w[_s]) * ds.node_infos[_s].sleep_power;
      }
      // 3. edge power on a per co basis
      for (int _l : ds.co_edge_ids[_c]) {
        // compute power for (u->v) only
        // not for both (u->v) and (v->u)
        if (ds.edges[_l].u <= ds.edges[_l].v) {
          switch_powers[_c] += 0.0012 * cap_0_1g[_l] +
                     0.0043 * cap_1g_inf[_l]; 
        }
      }
      //IloExpr brown_power(env);
      if (use_switch_cost_in_objective) {
        co_powers[_c] += switch_powers[_c];
      }

      co_brown_power[_c] = IloMax(0, co_powers[_c] -
          ds.renewable_energy[_c][timeslot]);
      // comment out the above eqn and use the following eqn
      // to ignore green energy
      //co_brown_power[_c] = co_powers[_c];
      total_brown_power += co_brown_power[_c];
      total_green_power += co_powers[_c] - co_brown_power[_c];
      objective += co_brown_power[_c] * ds.carbon_per_watt[_c];

      //cerr << _c << " " << ds.carbon_per_watt[_c] << endl;

    }
    // add power for backbone links
    // no green energy is considered here
    for (int _l : ds.backbone_edge_ids) {
      backbone_power += 0.0012 * cap_0_1g[_l] +
                 0.0043 * cap_1g_inf[_l]; 
    }
    total_brown_power += backbone_power;
    if (use_switch_cost_in_objective)
      objective += backbone_power * 1.12; //carbon_per_watt for backbone 1.12

    if (current_cost > 0) {
      model.add(objective <= (1 - migration_threshold) * current_cost);
    }
    
    /*Objective --> model*/
    model.add(objective >= 0);
    model.add(IloMinimize(env, objective));

    //=========Solver=========//
    //IloTimer timer(env);
    //timer.restart();

    const IloInt time_limit_seconds = 60*60; // 1 hour
    const IloNum relative_gap = 0.001; // 0.1% gap with optimal

    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::TiLim, time_limit_seconds);
    cplex.setParam(IloCplex::EpGap, relative_gap);
    cplex.setParam(IloCplex::PreDual, true);

    if(!cplex.solve()) {
      //timer.stop();
      /*
      if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        cout << "404 " << sfc << endl;
      }
      else { 
        cerr << "could not solve ILP!" << endl;
        cerr << "status: " << cplex.getStatus() << endl;
      }
      throw(-1);
      */
      cout << "404 " << sfc << endl;
      return 0;
    }

    double time = ctimer.time();

    //timer.stop();
    //cout << "ILP solved in " << timer.getTime() << " sec" << endl;
    //cout << "Objective value = " << cplex.getObjValue() << endl;
    double final_cost = cplex.getObjValue();
    if (!use_switch_cost_in_objective) {
      final_cost += cplex.getValue(backbone_power);
      for (int _c = 0; _c < ds.co_count; ++_c) {
        final_cost += cplex.getValue(switch_powers[_c]);
      }
    }
    cout << "200 " << sfc << " " << final_cost << " ";

    // get value for x and co stretch
    set<int> uniq_cos;
    cout << sfc.vnf_count << " ";
    for (int n = 1; n < sfc.node_count() - 1; ++n) {
      for (int _n = 0; _n < ds.node_count; ++_n) {
        //cout << "x[" << n << "][" << _n << "] = " <<
        //  cplex.getValue(x[n][_n]) << endl;
        if (IloRound(cplex.getValue(x[n][_n])) == 1) {
          cout << _n << " ";
          uniq_cos.insert(ds.node_infos[_n].co_id);
        }
      }
    }
    // output the co stretch
    cout << uniq_cos.size() << " ";
    // get value for y
    vector<vector<int>> all_paths;
    for (int l = 0; l < sfc.edge_count(); ++l) {
      //cout << "y[" << l << "] = ";
      for (int _p = 0; _p < ds.path_count; ++_p) {
        if (IloRound(cplex.getValue(y[l][_p])) == 1) {
          all_paths.push_back(ds.path_nodes[_p]);
        }
      }
    }
    cout << all_paths.size() << " ";
    for (auto& path_nodes : all_paths) {
      cout << path_nodes.size() << " ";
      copy(path_nodes.begin(), path_nodes.end(),
          ostream_iterator<int>(cout, " "));
    }
    cout << time << " " << cplex.getValue(total_brown_power) << 
      " " << cplex.getValue(total_green_power) << endl;

    // print per co cost
//    for (int _c = 0; _c < ds.co_count; ++_c) {
//      cerr << cplex.getValue(co_brown_power[_c]) << " " <<
//        cplex.getValue(co_brown_power[_c] * ds.carbon_per_watt[_c]) << " ";
//    }
//    cerr << cplex.getValue(backbone_power) << " " <<
//        cplex.getValue(backbone_power * 1.12) <<endl;


    // output to cerr for computing green energy utilization
//    cerr << timeslot << " " << sfc.id << " " << " " << ds.co_count << " ";
//    for (int _c = 0; _c < ds.co_count; ++_c) {
//      cerr << cplex.getValue(co_powers[_c]) << " ";
//    }
//    for (int _c = 0; _c < ds.co_count; ++_c) {
//      cerr << cplex.getValue(co_brown_power[_c]) << " ";
//    }
//    cerr << endl;
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
