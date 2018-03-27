#include "iz_topology.hpp"
#include "problem_instance.hpp"

// Stage - 1: Find a path between the ingress and egress COs with maximum
// aggregate green energy.
std::unique_ptr<izlib::iz_path> GetPathWithMaxGreenEnergy(
    const problem_instance& p_instance, const sfc_request& chain, int k, int timeslot) {
  std::unique_ptr<izlib::iz_path> ret_path(new izlib::iz_path());
  
  // Compute k-Shortest path between the ingress and egress CO of the chain.
  izlib::iz_path_list paths;
  auto& inter_co_topo = p_instance.topology.inter_co_topo;
  inter_co_topo.k_shortest_paths(chain.ingerss_co, chain.egress_co, k, paths, chain.bandwidth);

  // For each computed path, compute the aggregate amount of green energy
  // available on the COs.
  auto& cos = prob_inst.topology.cos;
  double max_avail_green = -1.0;
  int max_path_index = -1;
  for (int i = 0; i < paths.size(); ++i) {
    auto& path = paths[i];
    if (path.latency > chain.latency) continue;
    double total_avail_green = 0.0;
    for (auto co_id : path.nodes) {
      total_avail_green += cos[co_id].green_residual[timeslot]
    }
    if (total_avail_green > max_avail_green) {
      max_avail_green = total_avail_green;
      max_path_index = i;
    }
  }
  if (max_path_index != -1) {
    ret_path.reset(new iz_path(paths[max_path_index]));
  }
  return std::move(ret_path);
}

