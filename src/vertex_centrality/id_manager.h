#ifndef ID_MANAGER_H
#define ID_MANAGER_H

#include <vector>
#include <cstdlib>
using std::vector;

class IDManager {
  std::vector<int> alive_ids;
  std::vector<int> dead_ids;
  std::vector<int> pos_in_alive;
  std::vector<int> pos_in_dead;
public:
  IDManager(size_t V);
  void Add(bool alive = true);
  bool MakeAlive(int u);
  bool MakeDead(int u);
  int SampleAlive() const ;
  int SampleDead() const ;
  bool Full() const { return dead_ids.empty(); }
  size_t Size() const { return pos_in_alive.size(); }
  size_t NumAlive() const { return alive_ids.size(); }
};


#endif /* ID_MANAGER_H */
