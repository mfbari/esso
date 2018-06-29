#include <limits>

//#include "esso_heuristic.hpp"
#include "iz_topology.hpp"
#include "problem_instance.hpp"

using namespace std;
using namespace izlib;

template <typename T>
using matrix = vector<vector<T>>;

const char* join_path(string& dir, const string& file) {
  if (dir.back() != '/') dir += '/';
  return (dir+file).c_str();
}

iz_path stage_one(const sfc_request& sfc, const int k, const int timeslot,
    problem_instance& prob_inst) {
  // find path set
  iz_path_list paths;
  auto& inter_co_topo = prob_inst.topology.inter_co_topo;
  inter_co_topo.k_shortest_paths(sfc.ingress_co, sfc.egress_co, k, paths,
      sfc.bandwidth);
  
  // find the path with max energy
  int max_energy_path_index{-1};
  double max_energy{-1.0};
  const auto& cos = prob_inst.topology.cos;
  for (int i = 0; i < paths.size(); ++i) {
    const auto& path = paths[i];
    if (path.latency > sfc.latency) continue;
    double path_energy{0.0};
    for (auto& co_id : path.nodes) {
      path_energy += cos[co_id].green_residual[timeslot];
    }
    cout << path << " " << path_energy << endl;
    if (path_energy > max_energy) {
      max_energy_path_index = i;
      max_energy = path_energy;
    }
  }

  // return either the path with max energy or
  // return an empty path is all paths' latency
  // is greater than the sfc's latency bound
  iz_path path;
  if (max_energy_path_index != -1)
    path = move(paths[max_energy_path_index]);
  return path;
}

void stage_two(const int co_id, const sfc_request& sfc, const iz_path& path, 
    const int timeslot, 
    vector<vector<double>>& cost_matrix, vector<vector<int>>& node_matrix,
    problem_instance& prob_inst) {
  auto& cos = prob_inst.topology.cos;
  cout << "calling compute cost for co: " << co_id << endl;
  cos[co_id].compute_embedding_cost(sfc.cpu_reqs, sfc.bandwidth,
      timeslot, cost_matrix, node_matrix);
}

bool first_fit(const sfc_request& sfc, const iz_path& path, 
    const vector<vector<vector<double>>>& cost_matrices,
    vector<vector<int>>& embedding_table) {
  int next_vnf{0}, curr_co_idx{0};
  while (next_vnf < sfc.vnf_count && curr_co_idx < path.size()) {
    for (int i = sfc.vnf_count-1; i >= 0; --i) {
      if (cost_matrices[curr_co_idx][next_vnf][i] != -1) {
        for (int j = next_vnf; j <= i; ++j) {
          embedding_table[curr_co_idx][j] = 1;
        }
        next_vnf = i+1;
        break;
      }
    }
    if (next_vnf == sfc.vnf_count) return true;
    ++curr_co_idx;
  }
  return false;
}

// compute the embedding cost of an embedding table
double embedding_cost(const sfc_request& sfc, const iz_path& path,
    const vector<vector<vector<double>>>& cost_matrices,
    const vector<vector<int>>& embedding_table) {
  double cost{0.0};
  for (int c = 0; c < path.size(); ++c) {
    // find starting 1 in embedding table
    int sfc_start{0}, sfc_end{sfc.vnf_count-1};
    while(sfc_start < sfc.vnf_count && 
        embedding_table[c][sfc_start] == 0) 
      ++sfc_start;
    while(sfc_end >= 0 && 
        embedding_table[c][sfc_end] == 0)
      --sfc_end;
    if (sfc_start <= sfc_end) {
      cost += cost_matrices[c][sfc_start][sfc_end];
    }
  }
  return cost;
}

// returns the nodes selected in an embedding table
vector<int> embedding_nodes(const sfc_request& sfc, const iz_path& path,
    const vector<vector<vector<int>>>& node_matrices,
    const vector<vector<int>>& embedding_table) {
  vector<int> nodes;
  for (int c = 0; c < path.size(); ++c) {
    // find starting 1 in embedding table
    int sfc_start{0}, sfc_end{sfc.vnf_count-1};
    while(sfc_start < sfc.vnf_count && 
        embedding_table[c][sfc_start] == 0) 
      ++sfc_start;
    while(sfc_end >= 0 && 
        embedding_table[c][sfc_end] == 0)
      --sfc_end;
    for (int i = sfc_start; i <= sfc_end; ++i) {
      nodes.push_back(path.nodes[c] * 9 + 
          node_matrices[c][sfc_start][i]);
    }
  }
  return nodes;
}

// checks whether an embedding table is valid
bool is_valid_embedding(const sfc_request& sfc,const iz_path& path,
    const vector<vector<vector<double>>>& cost_matrices,
    const vector<vector<int>>& embedding_table) {
  // vnf_co stores the co_ids for the vnfs in the sfc
  // according to the provided embedding table
  vector<int> vnf_co(sfc.vnf_count);
  // loop for each vnf
  for (int vnf_idx = 0; vnf_idx < sfc.vnf_count; ++vnf_idx) {
    int sum{0}, non_zero_count{0};
    // chech each co for vnf vnf_inx
    for (int co_idx = 0; co_idx < path.size(); ++co_idx) {
      sum += embedding_table[co_idx][vnf_idx];
      if (embedding_table[co_idx][vnf_idx]) {
        ++non_zero_count;
        // here the entries in the vnf_co are updated to store
        // the co_ids on the path
        vnf_co[vnf_idx] = path.nodes[co_idx];
      }
    }
    // both sum and non_zero_count must be equal to one
    if (sum != 1 || non_zero_count != 1) return false;
  }

  // check to ensure vnf embedding always moves in the 
  // forward direction
  for (int v = 1; v < sfc.vnf_count; ++v) {
    if (vnf_co[v-1] > vnf_co[v]) {
      return false;
    }
  }

  // check to see if the 1's in the embedding table
  // are really valid costs in the cost_matrices
  for (int c = 0; c < path.size(); ++c) {
    // find starting 1 in embedding table
    int sfc_start{0}, sfc_end{sfc.vnf_count-1};
    while(sfc_start < sfc.vnf_count && 
        embedding_table[c][sfc_start] == 0) 
      ++sfc_start;
    // find ending 1 in the embedding table
    while(sfc_end >= 0 && 
        embedding_table[c][sfc_end] == 0)
      --sfc_end;
    if (sfc_start <= sfc_end) {
      if (cost_matrices[c][sfc_start][sfc_end] == -1)
        return false;
    }
  }


  return true;
}

bool tabu_search(const sfc_request& sfc, const iz_path& path,
    const vector<vector<vector<double>>>& cost_matrices,
    vector<vector<int>>& best_solution) {
    // initial solution from first-fit
    vector<vector<int>> current_solution = vector<vector<int>>(
        path.size(), vector<int>(sfc.vnf_count, 0));
    auto res = first_fit(sfc, path, cost_matrices, 
        current_solution);
    // if no fist-fit solution, then return false
    if (!res) return false;

    // overall best solution and cost
    double best_solution_cost {0.0};
    // set best solution = current solution
    best_solution = current_solution;
    best_solution_cost = embedding_cost(sfc, path, cost_matrices, 
        best_solution);

    // tabu search specific data strutures
    set<pair<int, int>> tabu_set;
    vector<vector<int>> tabu_timers(path.size(), 
        vector<int>(sfc.vnf_count, 0));
    const int tabu_period = 500;
    const int max_iterations = 100000;
    const int max_no_improvement_iterations = 1500;
    int best_cost_update_timestamp = 0;

    // variables for random number
    default_random_engine rnd_engine;
    uniform_int_distribution<> uni_dist(0,99);
    // main loop for tabu search
    int iter = 0;
    for (; iter < max_iterations; ++iter) {
      // datastructure to keep track of best neighbour
      double best_nbr_cost = numeric_limits<double>::max();
      vector<vector<int>> best_nbr;
      pair<int, int> potential_tabu_move(-1, -1);

      // generate neighbors to find the best neighbor for this iteration
      for (int j = 0; j < sfc.vnf_count; ++j) {
        // initialize the nbr solution with current solution
        auto nbr(current_solution);
        // now change one vnf assignment 
        // find the current index of co for vnf j
        int curr_co_idx{0}, next_co_idx{0};
        while (nbr[curr_co_idx][j] == 0 && 
            curr_co_idx < path.size()) 
          ++curr_co_idx;

        int rnd_num = uni_dist(rnd_engine);
        if (rnd_num%2) { // if odd then move up (-1)
          next_co_idx = (curr_co_idx + path.size() - 1) % path.size();
        }
        else { // even, so move down (+1)
          next_co_idx = (curr_co_idx + 1) % path.size();
        }
        // update current and next co assignment to generate the neighbor
        nbr[curr_co_idx][j] = 0;
        nbr[next_co_idx][j] = 1;

        // if nbr_solution is not valid then continue
        if (!is_valid_embedding(sfc, path, cost_matrices, nbr))
          continue;

        // find cost of neighbor solution
        auto nbr_cost = embedding_cost(sfc, path, cost_matrices, nbr);

        // update best neighbor and potential tabu move 
        if (nbr_cost < best_nbr_cost) {
          best_nbr_cost = nbr_cost;
          best_nbr = nbr;
          potential_tabu_move.first = j;
          potential_tabu_move.second = next_co_idx;
        }
      } // end of for loop for generating neigbor solutions

      // now, update best_solution with the best neighbor so far
      if (best_nbr_cost < best_solution_cost) {
        best_solution_cost = best_nbr_cost;
        best_solution = best_nbr;
        tabu_set.insert(potential_tabu_move);
        tabu_timers[potential_tabu_move.first][potential_tabu_move.second] =
            tabu_period;
        best_cost_update_timestamp = iter;
      }

      // update tabu list
      for(auto itr = tabu_set.begin(); itr != tabu_set.end();) {
        if (tabu_timers[itr->first][itr->second] > 0)
          --tabu_timers[itr->first][itr->second];
        if (tabu_timers[itr->first][itr->second] == 0)
          itr = tabu_set.erase(itr);
        else
          ++itr;
      }

      // if best cost is not updated in the last 
      // max_no_improvement_iterations then break
      if (iter - best_cost_update_timestamp > 
          max_no_improvement_iterations)
        break;

    } // end of tabu search iterations
    return true;
}

// prints the 404... message when no embedding is found
void no_path_found(const sfc_request& sfc) {
  cout << "404 " << sfc << endl;
}


int main(int argc, char **argv) {

  if (argc != 2) {
    cerr << "usage: ./esso_heuristic.o <relative-path-to-dataset-dir>" 
        << endl;
    return -1;
  }
  // directory for the dataset
  string dataset_dir {argv[1]};

  // prob_inst contains all the input data
  problem_instance prob_inst;
  problem_input prob_input;
  prob_input.vnf_info_filename = join_path(dataset_dir, "vnf_types.dat");
  prob_input.time_slot_filename = join_path(dataset_dir, "timeslots.dat");
  prob_input.topology_filename = join_path(dataset_dir, "co_topology.dat");

  // if all input read successfully, then call the heuristic
  if (prob_inst.read_input(prob_input)) {
    // call heuristic algorithm
    cout << "calling tabu search ..." << endl;

    sfc_request sfc;
    double current_cost, migration_threshold;
    cin >> sfc >> current_cost >> migration_threshold;

    int timeslot = 0;

    // Stage-1: find a path for embedding sfc
    int k = 3; // number of candidate paths
    auto embedding_path = stage_one(sfc, k, timeslot, prob_inst);
    // if no path found then print 404 ...
    if (!embedding_path.is_valid()) no_path_found(sfc);
    cout << "embedding path: " << embedding_path << endl;
    
    // Stage-2: compute the cost matrix for all co's on the embedding path
    vector<vector<vector<double>>> cost_matrices(embedding_path.size());
    vector<vector<vector<int>>> node_matrices(embedding_path.size());
    for (int i = 0; i < embedding_path.size(); ++i) {
      int co_id = embedding_path.nodes[i];
      stage_two(co_id, sfc, embedding_path, timeslot, 
          cost_matrices[i], node_matrices[i], prob_inst);
    }

    cost_matrices[0][0][2] = -1;
    cost_matrices[0][0][3] = -1;
    cost_matrices[1][2][2] = -1;
    cost_matrices[1][2][3] = -1;
    cost_matrices[2][2][3] = -1;

    // PRINT all cost and node matrices
    for (int i = 0; i < embedding_path.size(); ++i) {
      cout << "cost matrix[" << i << "]" << endl;
      for (int j = 0; j < sfc.vnf_count; ++j) {
        for (int k = 0; k < sfc.vnf_count; ++k) {
          cout << cost_matrices[i][j][k] << " ";
        }
        cout << endl;
      }
      cout << "node matrix[" << i << "]" << endl;
      for (int j = 0; j < sfc.vnf_count; ++j) {
        for (int k = 0; k < sfc.vnf_count; ++k) {
          cout << node_matrices[i][j][k] << " ";
        }
        cout << endl;
      }
    }

    // 1 1 3 1 4 1 2 2 1 100 100 0 0.2
    // Stage-3: call tabu search
    vector<vector<int>> embedding_table = vector<vector<int>>(
        embedding_path.size(), vector<int>(sfc.vnf_count, 0));
    auto res = tabu_search(sfc, embedding_path, cost_matrices, 
        embedding_table);
    cout << "embedding table " << boolalpha << res << endl;
    for (auto& row : embedding_table) {
      for (int c : row) {
        cout << c << " ";
      }
      cout << endl;
    }
    cout << embedding_cost(sfc, embedding_path, cost_matrices, 
        embedding_table) << endl;

    auto emb_nodes = embedding_nodes(sfc, embedding_path, 
        node_matrices, embedding_table);
    cout << "embedding nodes " << endl;
    for (auto node : emb_nodes)
      cout << node << " ";
    cout << endl;

    //embedding_table[3][3] = 0;
    //embedding_table[1][3] = 1;
    cout << "embedding table " << boolalpha << res << endl;
    for (auto& row : embedding_table) {
      for (int c : row) {
        cout << c << " ";
      }
      cout << endl;
    }
    cout << "check valid " <<  boolalpha <<
      is_valid_embedding(sfc, embedding_path, cost_matrices,
          embedding_table) << endl;

    //auto vnf_assignment = EssoHeuristic(prob_inst, sfc, 0);
    //for (auto va : *vnf_assignment.get()) {
    //  cout << va << " ";
    //}
    //cout << endl;
  }
  else {
    cerr << "failed to read input files for porblem instance" << endl;
    return -1;
  }

  return 0;
}
