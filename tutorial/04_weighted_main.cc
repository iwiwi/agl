#include "agl.h"
using namespace std;
using namespace agl;

int main() {
  //
  // 1. 重み付きグラフ
  //
  {
    // 重み付きグラフの辺の型は weighted_edge<重みの型> です．
    {
      weighted_edge<double> e{1, 2.5};
      cout << to(e) << ", " << weight(e) << endl;
    }

    // 重み付きグラフの辺リストの型は，|vector<pair<V, weighted_edge>>| です．
    vector<pair<V, weighted_edge<double>>> el = {
        {0, {1, 2.5}},
        {1, {2, 0.3}},
    };

    // 重み付きグラフは以下のように作ります．
    weighted_graph<double> g(el);
    pretty_print(g);
  }

  //
  // 2. 重みの操作
  //
  {
    /*
     重みの操作には，必ず |src/graph/weight_type.h| 内で定義される
     |is_zero|, |is_eq|, |is_lt| 等の関数を通じて行って下さい．

     これは，以下の理由によります．

     (1) double, float 等の浮動小数点数を用いた時，計算誤差を考慮した比較が必要にる．
         一方，int 等の整数型ではその必要がない．
         |weight_type.h| 内で定義した関数は型による場合分けにより
         それらの違いを吸収し自動的に適切な処理を行う．

     (2) コードの抽象性を上げ，ユーザが定義する重み付けを用いることができるようになる．
　       重みとしたい型に対しこれらのフリー関数を用意すれば，
         殆どの型を重みとして用いることができるようになる．
     */
  }

  return 0;
}
