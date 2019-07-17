#ifndef ESSO_HEURISTIC_HPP_
#define ESSO_HEURISTIC_HPP_

#include "iz_topology.hpp"
#include "problem_instance.hpp"

// Stage - 1: Find a path between the ingress and egress COs with maximum
// aggregate green energy.
std::unique_ptr<izlib::iz_path> GetPathWithMaxGreenEnergy(
    problem_instance& p_instance, sfc_request& chain, int k, int ts) {
  // the path that will be returned from this function
  std::unique_ptr<izlib::iz_path> ret_path(new izlib::iz_path());

  // Compute k-Shortest path between the ingress and egress CO of the chain.
  // iz_path_list paths will contain k paths, we will iterate over them
  // and select the path with the maximum green energy
  izlib::iz_path_list paths;
  auto& inter_co_topo = p_instance.topology.inter_co_topo;
  inter_co_topo.k_shortest_paths(chain.ingress_co, chain.egress_co, k, paths,
                                 chain.bandwidth);

  // For each computed path, compute the aggregate amount of green energy
  // available on the COs.
  auto& cos = p_instance.topology.cos;
  double max_avail_green = -1.0;
  int max_path_index = -1;
  for (int i = 0; i < paths.size(); ++i) {
    auto& path = paths[i];
    // if the path's letency is greater than the max-latency of the
    // chain, then we skip this path as it does not satisfy latency
    // constraint
    if (path.latency > chain.latency) continue;
    double total_avail_green = 0.0;
    for (auto co_id : path.nodes) {
      total_avail_green += cos[co_id].green_residual[ts];
    }
    if (total_avail_green > max_avail_green) {
      max_avail_green = total_avail_green;
      max_path_index = i;
    }
  }
  // if we have found a path and it satisfies the latency constraint,
  // then the max_path_index will be set to some other value than -1
  // we check this and reset the unique_ptr to the selcted path
  // if the max_path_index is still -1, then this function will return
  // an empty path
  if (max_path_index != -1) {
    ret_path.reset(new izlib::iz_path(paths[max_path_index]));
  }
  return std::move(ret_path);
}

// Stage - 2: Compute cost matrix containing cost of partially embedding a chain
// in a CO. Do that for all COs on the computed path.
std::unique_ptr<std::vector<std::vector<std::vector<double>>>>
ComputeCostMatrix(problem_instance& p_instance, const sfc_request& chain,
                  izlib::iz_path& path, int ts) {
  // cost_matrices is a 3D vector. Each entry of the first dimension represents
  // a cost_matrix corresponding to the CO on the path selected by Stage-1
  std::unique_ptr<std::vector<std::vector<std::vector<double>>>> cost_matrices(
      new std::vector<std::vector<std::vector<double>>>(path.size()));
  auto& cos = p_instance.topology.cos;
  for (int i = 0; i < path.size(); ++i) {
    auto co = path.nodes[i];
    // the co class provides a function to compute the cost_matrix
    std::vector<std::vector<double>> temp_matrix;
    cos[co].compute_embedding_cost(chain.cpu_reqs, chain.bandwidth, ts,
                                   temp_matrix);
    // individual cost_matrix is calculated and stored in the cost_matrices ds
    cost_matrices->at(i) = temp_matrix;
  }
  return std::move(cost_matrices);
}

// this function generates a first-fit embedding, which is used to
// initialize simulated annealing. The path selected by the stage-1 function
// is passed to this function and this function returns a vector of indexes
// into the path. If the chain has 3 vnfs then the returned vector will have
// 3 entries within the range [0-2], e.g., [0, 0, 1] indicates that the first
// two vnfs are embedded on the first node on the path and the last vnf is
// embedded on the second node on the path.
std::unique_ptr<std::vector<int>> GenerateFirstFitSolution(
    problem_instance& p_instance, sfc_request& chain, izlib::iz_path& path,
    std::vector<std::vector<std::vector<double>>>& cost_matrices) {
  int to_embed_start_idx = 0;
  // vnf_assignment stores the index of a path node
  // where a vnf is embedded
  std::unique_ptr<std::vector<int>> vnf_assignment(
      new std::vector<int>(chain.vnf_count, -1));
  // loop over the nodes on the path
  for (int i = 0; i < path.size(); ++i) {
    // if all the vnfs are embedded break out of the loop
    if (to_embed_start_idx == chain.vnf_count) break;

    auto& cost_matrix = cost_matrices[i];

    // find out the longest partial chain that can be embedded in the
    // current node
    int to_embed_end_idx = -1;
    for (int k = to_embed_start_idx; k < chain.vnf_count; ++k) {
      if (cost_matrix[to_embed_start_idx][k] == -1.0) break;
      to_embed_end_idx = k;
    }

    // if any part of the chain is embeddable, then update the index
    // values in the vector vnf_assignment
    if (to_embed_end_idx != -1) {
      for (int k = to_embed_start_idx; k <= to_embed_end_idx; ++k) {
        vnf_assignment->at(k) = i;
      }

      // update the start index for the next partial chain
      to_embed_start_idx = to_embed_end_idx + 1;
    }
  }
  return std::move(vnf_assignment);
}


// checks the validity of an vnf assignment
bool IsValidAssignment(
    const problem_instance& p_instance, const std::vector<int>& vnf_assignment,
    const std::vector<std::vector<std::vector<double>>>& cost_matrices,
    const izlib::iz_path& path) {
  // Each VNF should be assigned to a node on the path.
  // A -1 indicates that a VNF has not been assigned.
  for (int i = 0; i < vnf_assignment.size(); ++i)
    if (vnf_assignment[i] == -1) return false;

  // VNF i should be placed on a path node that precedes the placement
  // of VNF i + 1 on the path.
  for (int i = 1; i < vnf_assignment.size(); ++i) {
    if (vnf_assignment[i] < vnf_assignment[i - 1]) return false;
  }

  // If a segment of a chain (i ~> j) is placed in a path node k,
  // then the corresponding entry in cost matrix that was
  // computed part of Stage-2 should not be -1.
  int current_co = vnf_assignment[0];
  int seg_start = 0;
  int seg_end = 0;
  for (int i = 1; i < vnf_assignment.size(); ++i) {
    if (vnf_assignment[i] != current_co) {
      seg_end = i - 1;
      if (cost_matrices[current_co][seg_start][seg_end] == -1.0) return false;
      current_co = vnf_assignment[i];
      seg_start = seg_end = i;
    }
  }
  // check for the final leg of partial chain
  seg_end = vnf_assignment.size() - 1;
  if (cost_matrices[current_co][seg_start][seg_end] == -1.0) return false;

  return true;
}

double GetAssignmentCost(
    const std::vector<int>& vnf_assignment,
    const std::vector<std::vector<std::vector<double>>>& cost_matrices) {
  double ret = 0.0;
  int current_co = vnf_assignment[0];
  int seg_start = 0;
  int seg_end = 0;
  for (int i = 1; i < vnf_assignment.size(); ++i) {
    if (vnf_assignment[i] != current_co) {
      seg_end = i - 1;
      ret += cost_matrices[current_co][seg_start][seg_end];
      current_co = vnf_assignment[i];
      seg_start = seg_end = i;
    }
  }
  ret += cost_matrices[current_co][seg_start][seg_end];
  return ret;
}

// Stage - 3: Perform tabu search on VNF -> CO assignment matrix.
//
/*std::unique_ptr<std::vector<int>> TabuSearch(
    problem_instance& p_instance, sfc_request& chain, izlib::iz_path& path,
    std::vector<std::vector<std::vector<double>>>& cost_matrices,
    double& best_cost) {
  // Generate an initial solution by using first fit algorithm.
  auto initial_vnf_assignment =
      GenerateFirstFitSolution(p_instance, chain, path, cost_matrices);

  // First check if the initial solution is valid or not. If it is invalid then
  // declare embedding failed and return.
  if (!IsValidAssignment(p_instance, *initial_vnf_assignment.get(),
                         cost_matrices, path))
    return nullptr;

  // Initialize data structure to contain best solution. Also initialize best
  // solution cost.
  std::vector<int> current_vnf_assignment(*initial_vnf_assignment.get());
  std::unique_ptr<std::vector<int>> best_assignment(
      new std::vector<int>(current_vnf_assignment));
  best_cost = GetAssignmentCost(current_vnf_assignment, cost_matrices);

  std::set<std::pair<int, int>> tabu_set;
  std::vector<std::vector<int>> tabu_timer(
      path.size(), std::vector<int>(chain.vnf_count, 0));
  const int kTabuPeriod = 500;
  const int kMaxIterations = 100000;
  const int kMaxNonImprovedIterations = 1500;

  int best_cost_update_timestamp = 0;
  int iter = 0;
  for (; iter < kMaxIterations; ++iter) {
    double best_neighbor_cost = ~(1 << 31);
    std::unique_ptr<std::vector<int>> best_neighbor(nullptr);
    std::pair<int, int> potential_tabu_move(-1, -1);
    // Generate the neighbors of current solution. A neighbor is a valid
    // solution which differs with the current solution by one VNF placement.
    for (int j = 0; j < chain.vnf_count; ++j) {
      // Create a neighboring solution by changing the placement of one VNF.
      std::unique_ptr<vector<int>> neighbor(
          new std::vector<int>(current_vnf_assignment));
      neighbor->at(j) = (current_vnf_assignment[j] + 1) % path.size();

      // Now check the tabu set if the current move, i.e., vnf j to
      // neighbor->at(j) is in the tabu set or not. If yes, then do not make the
      // move.
      if (tabu_set.find(std::pair<int, int>(j, neighbor->at(j))) !=
          tabu_set.end())
        continue;

      // Check if making the move creates a valid solution or not.
      if (!IsValidAssignment(p_instance, *neighbor.get(), cost_matrices, path))
        continue;

      // Now updated best neighbor based on the cost and also create a potential
      // tabu move.
      double neighbor_cost = GetAssignmentCost(*neighbor.get(), cost_matrices);
      if (neighbor_cost < best_neighbor_cost) {
        best_neighbor_cost = neighbor_cost;
        best_neighbor = std::move(neighbor);
        potential_tabu_move.first = j;
        potential_tabu_move.second = best_neighbor->at(j);
      }
    }

    if (best_neighbor.get() != nullptr) {
      if (best_neighbor_cost < best_cost) {
        best_assignment = std::move(best_neighbor);
        best_cost = best_neighbor_cost;
        tabu_set.insert(potential_tabu_move);
        tabu_timer[potential_tabu_move.first][potential_tabu_move.second] =
            kTabuPeriod;
        best_cost_update_timestamp = iter;
      }
    }

    // Update tabu list.
    for (auto tabu_it = tabu_set.begin(); tabu_it != tabu_set.end();
         ++tabu_it) {
      if (tabu_timer[tabu_it->first][tabu_it->second] > 0) {
        tabu_timer[tabu_it->first][tabu_it->second]--;
        if (!tabu_timer[tabu_it->first][tabu_it->second]) {
          tabu_set.erase(tabu_it);
        }
      }
    }

    // If best cost is not updated in last kMaxNonImprovedIterations then break.
    if (iter - best_cost_update_timestamp > kMaxNonImprovedIterations) break;

    // If the move is in tabu list then continue;

    // Otherwise apply the move to the current solution and generate a
    // neighboring solution.

    ++iter;
  }
  return std::move(best_assignment);
}*/

// Heuristic driver funciton. This function calls the individual functions
// responsible for each of the stages and constructs the solution accordingly.
/*std::unique_ptr<std::vector<int>> EssoHeuristic(problem_instance& p_instance,
                                                sfc_request& chain, int ts,
                                                double& best_cost) {
  const int kNumShortestPaths = 3;
  auto embedding_path =
      GetPathWithMaxGreenEnergy(p_instance, chain, kNumShortestPaths, ts);
  auto cost_matrices =
      ComputeCostMatrix(p_instance, chain, *embedding_path.get(), ts);
  auto vnf_assignment = TabuSearch(p_instance, chain, *embedding_path.get(),
                                   *cost_matrices.get(), best_cost);
  return std::move(vnf_assignment);
}*/
#endif // ESSO_HEURISTIC_HPP_
