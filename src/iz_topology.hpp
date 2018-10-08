#ifndef IZ_TOPOLOGY
#define IZ_TOPOLOGY

#include <set>
#include <queue>
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <unordered_map>

#include <iostream>
//using namespace std;

#include "iz_priority_queue.hpp"

namespace izlib {

  typedef struct iz_edge iz_edge;
  using iz_edge_list = std::vector<iz_edge>;
  typedef struct iz_path iz_path;
  using iz_path_list = std::vector<iz_path>;
  using iz_node_list = std::vector<int>;

  class iz_topology {
      std::vector<std::unordered_map<int, iz_edge>> adj_matrix;
      void update_path_metrics(iz_path& path);
    public:
      int node_count, edge_count;
      explicit iz_topology();
      void init(int node_count);
      iz_node_list neighbors(int u);
      iz_edge add_edge(int u, int v, int latency, int capacity);
      iz_edge add_edge(int u, int v, int latency, int capacity, int residual);
      iz_edge_list remove_edge(int u, int v);
      iz_edge_list remove_node(int u);
      iz_edge_list edges() const;
      iz_edge_list edges(int u) const;
      iz_edge_list path_edges(const iz_path& path) const;
      int latency(int u, int v) const;
      int residual(int u, int v) const;
      void set_residual_bandwidth(int u, int v, int bandwidth);
      void allocate_bandwidth(int u, int v, int bandwidth);
      void release_bandwidth(int u, int v, int bandwidth);
      int consumed_bandwidth(int u, int v) const;
      void shortest_path(int s, int t, iz_path& path, 
          int min_capacity = 0);
      void k_shortest_paths(int s, int t, int K, iz_path_list& k_paths,
          int min_capacity = 0);
  };

  struct iz_edge {
    int u, v; // u is always less than v
    int latency;
    int capacity; // capacity capacity
    int residual; // residual capacity
    iz_edge(int u, int v, int latency, int capacity, int residual) :
        u(std::min(u,v)), v(std::max(u,v)), latency(latency), 
        capacity(capacity), residual(residual) {}
    int consumed_bandwidth() const {return capacity - residual;}
  };
  iz_edge_list operator+=(iz_edge_list& lhs, const iz_edge_list& rhs) {
    lhs.reserve(lhs.size() + rhs.size());
    lhs.insert(lhs.end(), rhs.begin(), rhs.end());
    return lhs;
  }

  struct iz_path {
    int latency;
    int capacity;
    iz_node_list nodes;
    iz_path(int latency = 0, int capacity = 0, 
        const iz_node_list nodes = std::vector<int>{}) :
        latency{latency}, capacity{capacity}, nodes{nodes} {}
    // used for priority queue
    bool operator<(const iz_path& rhs) const {
      return latency < rhs.latency;
    }
    bool operator>(const iz_path& rhs) const {
      return latency > rhs.latency;
    }
    iz_path operator+(const iz_path& rhs) {
      assert(nodes.back() == rhs.nodes.front());
      iz_path result;
      result.latency = latency + rhs.latency;
      result.capacity = std::min(capacity, rhs.capacity);
      result.nodes.reserve(nodes.size() + rhs.nodes.size());
      // assuming end node of this.nodes and first node of rhs.nodes 
      // are the same -- used in yen's k_shortest_path algorithm
      result.nodes.insert(result.nodes.end(), 
          nodes.begin(), nodes.end());
      result.nodes.insert(result.nodes.end(), 
          rhs.nodes.begin()+1, rhs.nodes.end());
      return result;
    }
    /*
    iz_path operator+=(const iz_path& rhs) {
      latency += rhs.latency;
      capacity = std::min(capacity, rhs.capacity);
      nodes.reserve(nodes.size() + rhs.nodes.size());
      // assuming end node of this.nodes and first node of rhs.nodes 
      // are the same -- used in yen's k_shortest_path algorithm
      nodes.insert(nodes.end(), rhs.nodes.begin()+1, rhs.nodes.end());
      return *this;
    }
    */
    bool is_equal_upto(const iz_path& rhs, size_t n) {
      if (size() <= n) return false;
      if (rhs.size() <= n) return false;
      for (size_t i = 0;i <= n; ++i) {
        if (nodes[i] != rhs.nodes[i]) return false;
      }
      return true;
    }
    size_t size() const {return nodes.size();}
    bool is_valid() const {return nodes.size() > 0;}
    int source() const {return nodes.front();}
    int target() const {return nodes.back();}
    void clear() {
      latency = 0;
      capacity = 0;
      nodes.clear();
    }
  };

  iz_topology::iz_topology() :
      node_count{0}, edge_count {0} {}

  void iz_topology::init(int node_count_) {
    node_count = node_count_;
    adj_matrix.resize(node_count, 
        std::unordered_map<int, iz_edge> {});
  }

  iz_node_list iz_topology::neighbors(int u) {
    assert(u >= 0 && u < node_count);
    iz_node_list node_list;
    for (int i = 0; i < u; ++i) {
      if (adj_matrix[i].find(u) != adj_matrix[i].end()) {
        node_list.push_back(i);
      }
    } 
    for (auto& v_edge : adj_matrix[u]) {
      node_list.push_back(v_edge.first);
    }
    return node_list;
  }

  iz_edge iz_topology::add_edge(int u,  int v, 
      int latency, int capacity, int residual) {
    assert(u >= 0 && u < node_count);
    assert(v >= 0 && v < node_count);
    assert(latency >= 0 && capacity >= 0);
    //std::tie(u, v) = std::minmax({u,v});
    if (u > v) std::swap(u, v);
    iz_edge edge(u, v, latency, capacity, residual);
    adj_matrix[u].insert(std::make_pair(v, edge));
    ++edge_count;
    return edge;
  }


  iz_edge iz_topology::add_edge(int u,  int v, 
      int latency, int capacity) {
    return add_edge(u, v, latency, capacity, capacity);
  }

  iz_edge_list iz_topology::remove_edge(int u, int v) {
    assert(u >= 0 && u < node_count);
    assert(v >= 0 && v < node_count);
    iz_edge_list edge_list;
    std::tie(u, v) = std::minmax({u,v});
    auto edge_itr = adj_matrix[u].find(v);
    if (edge_itr == adj_matrix[u].end()) {
      return edge_list;
    }
    edge_list.push_back(edge_itr->second);
    adj_matrix[u].erase(v);
    --edge_count;
    return edge_list;
  }

  iz_edge_list iz_topology::remove_node(int u) {
    assert(u >= 0 && u < node_count);

    iz_edge_list edge_list;
    for (int i = 0; i < u; ++i) {
      if (adj_matrix[i].find(u) != adj_matrix[i].end()) {
        edge_list.push_back(adj_matrix[i].find(u)->second);
        adj_matrix[i].erase(u);
        --edge_count;
      }
    }
    for (auto& v_edge : adj_matrix[u]) {
      edge_list.push_back(v_edge.second);
      --edge_count;
    }
    adj_matrix[u].clear();
    return edge_list;
  }

  iz_edge_list iz_topology::edges() const {
    iz_edge_list edge_list;
    for (int u = 0; u < node_count; ++u) {
      for (auto& v_edge : adj_matrix[u]) {
        edge_list.push_back(v_edge.second);
      }
    }
    return edge_list;
  }

  iz_edge_list iz_topology::edges(int u) const {
    iz_edge_list edge_list;
    for (int i = 0; i < u; ++i) {
      auto edge_itr = adj_matrix[i].find(u);
      if (edge_itr != adj_matrix[i].end()) {
        edge_list.push_back(edge_itr->second);
      }
    }
    for (auto& v_edge : adj_matrix[u]) {
      edge_list.push_back(v_edge.second);
    }
    return edge_list;
  }
    

  iz_edge_list iz_topology::path_edges(const iz_path& path) const {
    iz_edge_list edges;
    if (path.size() == 0) return edges;
    int u = path.nodes[0], v, _u, _v;
    assert(u >= 0 && u < node_count);
    for (size_t i = 1; i < path.size(); ++i) {
      v = path.nodes[i];
      assert(v >= 0 && v < node_count);
      _u = std::min(u, v);
      _v = std::max(u, v);
      edges.emplace_back(adj_matrix[_u].find(_v)->second);
      u = v; 
    }
    return edges;
  }

  int iz_topology::latency(int u, int v) const {
    assert(u >= 0 && u < node_count);
    assert(v >= 0 && v < node_count);
    //std::tie(u, v) = std::minmax({u,v});
    if (u > v) std::swap(u, v);
    auto edge_itr = adj_matrix[u].find(v);
    if (edge_itr == adj_matrix[u].end()) 
      return std::numeric_limits<int>::max();
    return edge_itr->second.latency;
  }

  int iz_topology::residual(int u, int v) const {
    assert(u >= 0 && u < node_count);
    assert(v >= 0 && v < node_count);
    std::tie(u, v) = std::minmax({u,v});
    auto edge_itr = adj_matrix[u].find(v);
    if (edge_itr == adj_matrix[u].end()) return 0;
    return edge_itr->second.residual;
  }

  void iz_topology::set_residual_bandwidth(int u, int v, int bandwidth) {
    assert(u >= 0 && u < node_count);
    assert(v >= 0 && v < node_count);
    assert(bandwidth >= 0);
    std::tie(u, v) = std::minmax({u,v});
    auto edge_itr = adj_matrix[u].find(v);
    if (edge_itr != adj_matrix[u].end()) {
      edge_itr->second.capacity = bandwidth;
      edge_itr->second.residual = bandwidth;
    }
  }

  void iz_topology::allocate_bandwidth(int u, int v, int bandwidth) {
    assert(u >= 0 && u < node_count);
    assert(v >= 0 && v < node_count);
    assert(bandwidth >= 0);
    std::tie(u, v) = std::minmax({u,v});
    auto edge_itr = adj_matrix[u].find(v);
    if (edge_itr != adj_matrix[u].end()) {
      edge_itr->second.residual -= bandwidth;
    }
  }

  void iz_topology::release_bandwidth(int u, int v, int bandwidth) {
    assert(u >= 0 && u < node_count);
    assert(v >= 0 && v < node_count);
    assert(bandwidth >= 0);
    std::tie(u, v) = std::minmax({u,v});
    auto edge_itr = adj_matrix[u].find(v);
    if (edge_itr != adj_matrix[u].end()) {
      edge_itr->second.residual += bandwidth;
    }
  }

  int iz_topology::consumed_bandwidth(int u, int v) const {
    assert(u >= 0 && u < node_count);
    assert(v >= 0 && v < node_count);
    std::tie(u, v) = std::minmax({u,v});
    auto edge_itr = adj_matrix[u].find(v);
    if (edge_itr == adj_matrix[u].end()) return 0;
    return edge_itr->second.consumed_bandwidth();
  }

  void iz_topology::shortest_path(int s, int t, iz_path& path, 
      int min_capacity) {
    std::vector<int> dist(node_count, std::numeric_limits<int>::max());
    std::vector<int> capacity(node_count, std::numeric_limits<int>::max());
    std::vector<int> parent(node_count, -1);
    std::vector<bool> visited(node_count, false);

    dist[s] = 0;
    parent[s] = s;

    iz_priority_queue<> pq;
    for (int u = 0; u < node_count; ++u) {
      pq.push(u, dist[u]);
    }

    bool found_t = false;
    while (!found_t && !pq.empty()) {
    //while (!pq.empty()) {
      int u = pq.top();
      pq.pop();
      for (int v : neighbors(u)) {
        if (visited[v]) continue;
        if (residual(u, v) < min_capacity) continue;
        if (parent[u] != -1 && dist[v] > dist[u] + latency(u, v)) {
          dist[v] = dist[u] + latency(u, v);
          capacity[v] = std::min(capacity[u], residual(u, v));
          parent[v] = u;
          pq.update_key(v, dist[v]);
        }
        if (v == t) {
          found_t = true;
          break;
        }
      }
      visited[u] = true;
    }

    path.clear();
    if (parent[t] != -1) {
      path.latency = dist[t];
      path.capacity = capacity[t];
      path.nodes.push_back(t);
      int _p = t;
      while (_p != s) {
        _p = parent[_p];
        path.nodes.push_back(_p);
      }
      std::reverse(path.nodes.begin(), path.nodes.end());
    }
  }

  void iz_topology::k_shortest_paths(int s, int t, int K, 
      iz_path_list& k_paths, int min_capacity) {
    assert(s >= 0 && s < node_count);
    assert(t >= 0 && s < node_count);
    assert(K > 0);

    k_paths.clear();

    iz_path s_path;
    shortest_path(s, t, s_path, min_capacity);
    if (!s_path.is_valid()) return;
    k_paths.push_back(s_path);
      
    //std::priority_queue<iz_path> pq;
    std::set<std::vector<int>> pq_set;
    std::priority_queue<iz_path, std::vector<iz_path>, 
        std::greater<iz_path>> pq;
   
    for (int k = 1; k < K; ++k) {
      for (size_t i = 0; i < k_paths[k-1].size()-1; ++i) {
        int spur_node = k_paths[k-1].nodes[i];
        iz_path root_path;
        root_path.clear();
        root_path.nodes.insert(root_path.nodes.end(), 
            k_paths[k-1].nodes.begin(), k_paths[k-1].nodes.begin()+i+1);

        iz_edge_list removed_edges;

        for (auto& path : k_paths) {
          if (root_path.is_equal_upto(path, i)) {
            removed_edges += remove_edge(path.nodes[i], 
                path.nodes[i+1]);
          }
        }

        for (int u : root_path.nodes) {
          if (u != spur_node) {
            removed_edges += remove_node(u);
          }
        }

        s_path.clear();
        shortest_path(spur_node, t, s_path, min_capacity);
        for (auto& edge : removed_edges) {
          add_edge(edge.u, edge.v, edge.latency, edge.capacity, edge.residual);
        }

        if (!s_path.is_valid()) continue;
        iz_path total_path = root_path + s_path;
        update_path_metrics(total_path);

        if (pq_set.count(total_path.nodes) == 0) {
          pq.push(total_path);
          pq_set.insert(total_path.nodes);
        }
      }
      if (pq.empty()) break;
      iz_path kth_path = pq.top();
      k_paths.push_back(kth_path);
      pq.pop();
      pq_set.erase(kth_path.nodes);
    }
  
  }

  void iz_topology::update_path_metrics(iz_path& path) {
    path.latency = 0;
    path.capacity = std::numeric_limits<int>::max();
    if (path.size() == 0) return;
    int u = path.nodes[0], v;
    for (size_t i = 1; i < path.size(); ++i) {
      v = path.nodes[i];
      path.latency += latency(u, v);
      path.capacity = std::min(path.capacity, residual(u, v));
      u = v; 
    }
  }

} // end of namespace izlib

std::ostream& operator<<(std::ostream& out, const izlib::iz_path& path) {
  for (auto& n : path.nodes) out << n << " ";
  out << "l:" << path.latency << " ";
  out << "c:" << path.capacity;
  return out;
}

#endif
