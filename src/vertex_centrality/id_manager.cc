#include "id_manager.h"
#include "common.h"
#include <numeric>
using namespace std;

IDManager::IDManager(size_t V) : alive_ids(V), pos_in_alive(V), pos_in_dead(V) {
  iota(alive_ids.begin(), alive_ids.end(), 0);
  iota(pos_in_alive.begin(), pos_in_alive.end(), 0);
  fill(pos_in_dead.begin(), pos_in_dead.end(), -1);
}

void IDManager::Add(bool alive) {
  size_t n = pos_in_alive.size();

  if (alive){
    alive_ids.push_back(n);
    pos_in_alive.push_back(alive_ids.size() - 1);
    pos_in_dead.push_back(-1);
  } else {
    dead_ids.push_back(n);
    pos_in_alive.push_back(-1);
    pos_in_dead.push_back(alive_ids.size() - 1);
  }
}

bool IDManager::MakeAlive(int u){
  while ((size_t)u >= Size()){
    Add(false);
  }

  if (pos_in_alive[u] == -1){
    
    if (dead_ids.size() > 1u && dead_ids.back() != u){
      int l = dead_ids.back();
      CHECK(l == dead_ids[pos_in_dead[l]]);
      swap(dead_ids.back(), dead_ids[pos_in_dead[u]]);
      pos_in_dead[l] = pos_in_dead[u];
      pos_in_dead[u] = dead_ids.size() - 1;
    }

    alive_ids.push_back(u);
    dead_ids.pop_back();
    pos_in_dead[u]  = -1;
    pos_in_alive[u] = alive_ids.size() - 1;
    return true;
  } else {
    return false;
  }
}

bool IDManager::MakeDead(int u){
  if (pos_in_dead[u] == -1){
  
    // pop u from alive_ids.
    if (alive_ids.size() > 1u && alive_ids.back() != u){
      int l = alive_ids.back();
      CHECK(l == alive_ids[pos_in_alive[l]]);
      swap(alive_ids.back(), alive_ids[pos_in_alive[u]]);
      pos_in_alive[l] = pos_in_alive[u];
      pos_in_alive[u] = alive_ids.size() - 1;
    }

    // cout << "POP" << endl;
    alive_ids.pop_back();
    dead_ids.push_back(u);
    pos_in_dead[u]  = dead_ids.size() - 1;
    pos_in_alive[u] = -1;
    return true;
  } else {
    return false;
  }
}

int IDManager::SampleAlive() const {
  return !alive_ids.empty() ? alive_ids[rand() % alive_ids.size()] : -1;
}

int IDManager::SampleDead() const {
  return !dead_ids.empty() ? dead_ids[rand() % dead_ids.size()] : -1;
}
