#include <iostream>
#include <fstream>
#include <iterator>
#include <sys/stat.h>

#include "problem_instance.hpp"
#include "iz_topology.hpp"
#include "stop_watch.hpp"

using namespace std;

const char* join_path(string& dir, const string& file) {
  if (dir.back() != '/') dir += "/";
  return (dir+file).c_str();
}

bool write_init_topology(string& dataset_dir, 
    problem_instance& prob_inst,
    izlib::iz_topology& topo, 
    vector<char>& node_info, vector<int>& server_ids) {
  cout << "Generating init_topology.dat ... ";
  ofstream fout(join_path(dataset_dir, "init_topology.dat"));
  if (!fout) {
    cerr << "failed to create init_topology.dat file in " +
      dataset_dir << endl;
    return false;
  }
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
  // resize the node_info vector
  node_info.resize(total_node_count);
  // ds to hold id mappings
  vector<vector<int>> id_map(cos.size(), vector<int>(max_inodes));
  int node_id{0}; // node_id for the merged graph
  fout << total_node_count << " " << total_edge_count << endl;
  //fout << cos.size() << " " << inter_co_topo.edge_count << endl;
  for (auto& co : cos) {
    // output the switches 
    for (int i = 0; i < co.server_ids[0]; ++i) {
      id_map[co.id][i] = node_id;
      node_info[node_id] = 's';
      fout << node_id++ << " s " << co.id << " " <<
        co.intra_nodes[i]->sleep_power << " " <<
        co.intra_nodes[i]->base_power << endl;
    }
    // now the servers
    for (auto i : co.server_ids) {
      id_map[co.id][i] = node_id;
      node_info[node_id] = 'c';
      server_ids.push_back(node_id);
      auto server_ptr = dynamic_pointer_cast<esso_server>(
          co.intra_nodes[i]);
      fout << node_id++ << " c " << co.id << " " <<
        co.intra_nodes[i]->sleep_power << " " <<
        co.intra_nodes[i]->base_power << " " << 
        server_ptr->cpu_capacity << " " <<
        server_ptr->per_cpu_power << endl;
    }
  }
  // initialize topo
  topo.init(total_node_count);
  // output inter co edges
  for (auto& e : inter_co_topo.edges()) {
    fout << id_map[e.u][0] << " " << id_map[e.v][0] << " b " << 
      e.capacity << " " << e.latency << endl;
    fout << id_map[e.v][0] << " " << id_map[e.u][0] << " b " << 
      e.capacity << " " << e.latency << endl;
    topo.add_edge(id_map[e.u][0], id_map[e.v][0], e.latency, e.capacity);
    topo.add_edge(id_map[e.v][0], id_map[e.u][0], e.latency, e.capacity);
  }
  // output intra co edges
  for (auto& co : cos) {
    // output edges
    for (auto& e : co.intra_topo.edges()) {
      fout << id_map[co.id][e.u] << " " << id_map[co.id][e.v] << 
        " i " << e.capacity << " " << e.latency << endl;
      fout << id_map[co.id][e.v] << " " << id_map[co.id][e.u] << 
        " i " << e.capacity << " " << e.latency << endl;
      topo.add_edge(id_map[co.id][e.u], id_map[co.id][e.v], 
          e.latency, e.capacity);
      topo.add_edge(id_map[co.id][e.v], id_map[co.id][e.u], 
          e.latency, e.capacity);
    }
  }

  //---------------

  fout.close();
  cout << "done" << endl;
  return true;
}

bool write_path_link(string& dataset_dir, izlib::iz_topology& topo, 
    const vector<char>& node_info, vector<int> server_ids,
    bool use_one_path = true, int phy_k = 3) {
  auto node_count = node_info.size();
  izlib::iz_path_list phy_paths;
  for (int u = 0; u < node_count; ++u) {
    for (int v = 0; v < node_count; ++v) {
      if (u != v && node_info[u] == 'c' &&
          node_info[v] == 'c') {
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

  // write to file
  int server_count = server_ids.size();
  ofstream fout(join_path(dataset_dir, "paths.dat"));
  fout << phy_paths.size() + server_count << endl;
  // insert dummy paths for in-server embedding
  for (int i = 0; i < server_ids.size(); ++i) {
    fout << i << " 2 " << server_ids[i] << " " 
        << server_ids[i] << " 0" << endl;
  }
  // now the actual paths
  for (size_t i = 0; i < phy_paths.size(); ++i) {
    // path_switch.dat file
    fout << server_count + i << " " << phy_paths[i].nodes.size() << " ";
    vector<int> switches;
    for (auto& node : phy_paths[i].nodes) {
      fout << node << " ";
      if (node_info[node] == 's') switches.push_back(node);
    }
    fout << switches.size() << " ";
    copy(switches.begin(), switches.end(), ostream_iterator<int>(
        fout, " "));
    fout<< endl;
  }
  fout.close();
  return true;
}

bool process_dataset(string& dataset_dir) {
  problem_instance prob_inst;

  problem_input prob_input;
  prob_input.vnf_info_filename = join_path(dataset_dir, "vnf_types.dat");
  prob_input.time_slot_filename = join_path(dataset_dir, "timeslots.dat");
  prob_input.topology_filename = join_path(dataset_dir, "co_topology.dat");

  if (prob_inst.read_input(prob_input)) {
    // if all input files read successfully, then ...
    izlib::iz_topology topo; // holds the entire topology (inter-co + intra-co)
    vector<char> node_info; // to differentiate between server & switch
    vector<int> server_ids;
    return write_init_topology(dataset_dir, prob_inst, topo, 
        node_info, server_ids) &&
      write_path_link(dataset_dir, topo, node_info, server_ids);
  }
  // error reading input file(s)
  return false;
}

int main(int argc, char**argv) {

  // check for the correct usage
  if (argc != 2) {
    cerr << "usage: ./generate_full_topology.o <dataset-dir-name>" << endl;
    return -1;
  }

  string dataset_dir {argv[1]};
  //try to open co_topology.dat file
  stop_watch sw;
  sw.start();
  auto result = process_dataset(dataset_dir);
  sw.stop();
  cout << "time: " << sw << endl;
  if (!result) cerr << "error" << endl;
  return 0;
}
