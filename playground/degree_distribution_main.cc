#include "easy_cui.h"
#include "picojson.h"
using namespace picojson;
DEFINE_string(jlog_file, "", "JLOG file path");

string load_json(string& f) {
  ifstream ifs(f);
  if (ifs.fail()) {
    cerr << "faled" << endl;
    return "";
  }
  int begin = static_cast<int>(ifs.tellg());
  ifs.seekg(0, ifs.end);
  int end = static_cast<int>(ifs.tellg());
  int size = end - begin;
  ifs.clear();
  ifs.seekg(0, ifs.beg);
  char* str = new char[size + 1];
  str[size] = '\0';
  ifs.read(str, size);
  return str;
}

int main(int argc, char** argv) {
  G g = easy_cui_init(argc, argv);
  CHECK_MSG(FLAGS_force_undirected, "undirected only!!!");

  string json = load_json(FLAGS_jlog_file);
  value v;
  string err;
  parse(v, json.begin(), json.end(), &err);
  picojson::array& centers_array =
      v.get<object>()["centers"].get<picojson::array>();
  for (auto it = centers_array.begin(); it != centers_array.end(); it++) {
    object& co = it->get<object>();
    for (auto co_it = co.begin(); co_it != co.end(); co_it++) {
      auto rad = co_it->first;
      picojson::array centers = co[rad].get<picojson::array>();
    }
  }

  return 0;
}
