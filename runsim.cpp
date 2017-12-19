#include <iostream>
#include <fstream>
#include <string>
#include <limits>

#include "problem_instance.hpp"
#include "iz_topology.hpp"

using namespace std;

void embed(const int time_slot, problem_instance& prob_inst, 
    sfc_request& sfc_req) {
  // find path between ingress and egress cos
  auto& inter_co_topo = prob_inst.topology.inter_co_topo;
  auto& cos = prob_inst.topology.cos;
  cout << sfc_req.ingress_co << " - " << sfc_req.egress_co << endl;
  izlib::iz_path_list paths;
  inter_co_topo.k_shortest_paths(sfc_req.ingress_co, sfc_req.egress_co, 
      3, paths, sfc_req.bandwidth);
  for (auto& path : paths) {
    cout << path << endl;
  }
  auto& path = paths[0];
  vector<vector<vector<double>>> cost_matrices(path.size());
  int cm_index{0};
  for (auto co_id : path.nodes) {
    cout << co_id << ":" << cos[co_id].green_residual[time_slot] << endl;
    cos[co_id].compute_embedding_cost(sfc_req.cpu_reqs, sfc_req.bandwidth, 
        time_slot, cost_matrices[cm_index++]);
  }

  cout << "---" << endl;
  for (auto& cost_matrix : cost_matrices) {
    for (auto& row : cost_matrix) {
      for (auto& val : row) {
        cout << val << " ";
      }
      cout << endl;
    }
    cout << "---" << endl;
  }
}


void write_topo_to_file(problem_instance& prob_inst) {
  cout << "Generating phy_inf_cplex.dat ... ";
  ofstream fout("phy_inf_cplex.dat");
  
  //---------------
  auto& inter_co_topo = prob_inst.topology.inter_co_topo;
  auto& cos = prob_inst.topology.cos;

  int total_node_count{0}, total_edge_count{0};
  int max_inodes{0};
  //total_node_count += cos.size();
  total_edge_count += inter_co_topo.edge_count;
  fout << cos.size() << endl;
  for (auto& co : cos) {
    total_node_count += co.intra_topo.node_count;
    total_edge_count += 2*co.intra_topo.edge_count;
    max_inodes = max(max_inodes, co.intra_topo.node_count);
    fout << co.id << " ";
    for(auto& ge : co.green_capacity) {
      fout << ge << " ";
    }
    fout << endl;
  }
  // ds to hold id mappings
  vector<vector<int>> id_map(cos.size(), vector<int>(max_inodes));
  int node_id{0}; // node_id for the merged graph
  fout << total_node_count << " " << total_edge_count << endl;
  //fout << cos.size() << " " << inter_co_topo.edge_count << endl;
  for (auto& co : cos) {
    // output the switches 
    for (int i = 0; i < co.server_ids[0]; ++i) {
      id_map[co.id][i] = node_id;
      fout << node_id++ << " s " << co.id << " " <<
          co.intra_nodes[i]->sleep_power << " " <<
          co.intra_nodes[i]->base_power << endl;
    }
    // now the servers
    for (auto i : co.server_ids) {
      id_map[co.id][i] = node_id;
      auto server_ptr = dynamic_pointer_cast<esso_server>(
          co.intra_nodes[i]);
      fout << node_id++ << " c " << co.id << " " <<
          co.intra_nodes[i]->sleep_power << " " <<
          co.intra_nodes[i]->base_power << " " << 
          server_ptr->cpu_capacity << " " <<
          server_ptr->per_cpu_power << endl;
    }
  }
  for (auto& e : inter_co_topo.edges()) {
    fout << id_map[e.u][0] << " " << id_map[e.v][0] << " b " << e.capacity << " " << e.latency << endl;
    fout << id_map[e.v][0] << " " << id_map[e.u][0] << " b " << e.capacity << " " << e.latency << endl;
  }
  for (auto& co : cos) {
    // output edges
    for (auto& e : co.intra_topo.edges()) {
      fout << id_map[co.id][e.u] << " " << id_map[co.id][e.v] << " i " << e.capacity << " " << e.latency << endl;
      fout << id_map[co.id][e.v] << " " << id_map[co.id][e.u] << " i " << e.capacity << " " << e.latency << endl;
    }
  }
  
  //---------------
  
  fout.close();
  cout << "done" << endl;
}

int main(int argc, char** argv) {

  problem_instance prob_inst;

  problem_input prob_input;
  prob_input.vnf_info_filename = "vnfinfo.dat";
  prob_input.time_slot_filename = "timeslots.dat";
  prob_input.topology_filename = "topo.dat";

  prob_inst.read_input(prob_input);

  if (argc > 1) {
    write_topo_to_file(prob_inst);
    return 0;
  }

  for (sfc_request_set& sfcs: prob_inst.time_slots) {
    for (sfc_request& sfc : sfcs) {
      cout << sfc << endl;
      embed(0, prob_inst, sfc);
      break;
    }
    break;
  }

  /*
  double carbon{0.0};
  for (auto& c: prob_inst.topology.cos) {
    for (int i = 0; i < 24; ++i) {
      carbon += c.get_carbon_fp(i);
    }
  }
  cout << carbon << endl;

  vector<vector<vector<double>>> cost_matrix(prob_inst.topology.cos.size());
  for (size_t c = 0; c < prob_inst.topology.cos.size(); ++c) {
    prob_inst.topology.cos[c].compute_embedding_cost(
        std::vector<int> {4, 4, 2}, 500, 0, cost_matrix[c]);
  }
  
  double min_cost {numeric_limits<double>::max()};
  for (auto& matrix : cost_matrix) {
    cout << matrix[0][2] << endl;
    min_cost = min(min_cost, matrix[0][2]);
  }
  cout << min_cost << endl;
  */

  return 0;
}
