#ifndef IZ_PRIORITY_QUEUE_H
#define IZ_PRIORITY_QUEUE_H

#include <vector>
#include <unordered_map>

//#include <iostream>
//using namespace std;

namespace izlib {

  template <class item_type = int, 
      class key_type = int,
      class compare = std::less<key_type>,
      class item_type_hash = std::hash<item_type>>
  class iz_priority_queue {
    public:
      struct pq_item {
        item_type item;
        key_type key;
        pq_item(const item_type& item, const key_type& key) :
            item(item), key(key) {}
      }; 

      std::vector<pq_item> heap;
      std::unordered_map<item_type, int, item_type_hash> 
          item_to_heap_index;

      key_type key(const int heap_index) const {
        return heap[heap_index].key;
      }
      item_type item(const int heap_index) const {
        return heap[heap_index].item;
      }

      void heap_swap(const int x, const int y) {
        std::swap(heap[x], heap[y]);
        item_to_heap_index.erase(item(x));
        item_to_heap_index.erase(item(y));
        item_to_heap_index.insert(
            std::make_pair(item(x), x));
        item_to_heap_index.insert(
            std::make_pair(item(y), y));
      }
    
      void shift_up(const int heap_index) {
        int p{parent(heap_index)}, larger{heap_index};
        if (p >= 0 && compare()(key(larger), key(p))) {
          larger = p;
        }
        if (larger != heap_index) {
          heap_swap(larger, heap_index);
          shift_up(larger);
        }
      }
      void shift_down(const int heap_index) {
        int lc{left_child(heap_index)}, rc{right_child(heap_index)},
            smaller{heap_index};
        if (lc < size() && compare()(key(lc), key(smaller))) {
          smaller = lc;
        }
        if (rc < size() && compare()(key(rc), key(smaller))) {
          smaller = rc;
        }
        if (smaller != heap_index) {
          heap_swap(smaller, heap_index);
          shift_down(smaller);
        }
      }

      int left_child(const int heap_index) const {
        return 2 * heap_index + 1;
      }
      int right_child(const int heap_index) const {
        return 2 * heap_index + 2;
      }
      int parent(const int heap_index) const {
        return (heap_index - 1) / 2;
      }
    public:
      void push(const item_type& item, const key_type& key) {
        int item_index = heap.size();
        heap.emplace_back(pq_item(item, key));
        item_to_heap_index.insert(std::make_pair(item, item_index));
        shift_up(item_index);
      }
      item_type top() {
        return heap.front().item;
      }
      key_type top_key() {
        return heap.front().key;
      }
      void pop() {
        heap_swap(0, heap.size() - 1);
        item_to_heap_index.erase(heap.back().item);
        heap.pop_back();
        shift_down(0);
      }

      bool update_key(const item_type& item, const key_type& new_key) {
        auto item_index_itr = item_to_heap_index.find(item);
        if (item_index_itr == item_to_heap_index.end()) return false;
        int item_index = item_index_itr->second;
        key_type old_key = key(item_index);
        heap[item_index].key = new_key;
        if (compare()(old_key, new_key)) {
          shift_down(item_index);
        }
        else {
          shift_up(item_index);
        }
        return true;
      } 
    
      int size() const {return heap.size();}
      bool empty() const {return size() == 0;}
  };
}
#endif
