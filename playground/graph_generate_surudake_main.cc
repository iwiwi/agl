#include "easy_cui.h"
#include "box_cover.h"
using namespace std;
using namespace agl;

int main(int argc, char** argv) {
  auto es = generate_uv_flower(2732, 2, 2);
  cout << es.size() << endl;
  for (auto e : es) {
    cout << e.first << " " << e.second << endl;
  }
}