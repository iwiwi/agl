#include "agl.h"
using namespace std;
using namespace agl;

//
// 1. AGL の設計思想
//
/*
* 全ての関数をテンプレートで記述し，全体に高い抽象度を実現することは，かなりの手間です．
* 一方，重みが有る場合・ない場合，異なるおもみ型の場合でも同じ処理がしたくなることも少なくありません．
*
* AGL では，異なる抽象度のコードの共存を許し，必要に応じた容易な抽象化を実現することにより，
* 上述の 2 つの問題を解決します．
*
*
* (1) プロトタイピング時，及び使用頻度の低いアルゴリズムの実装時
*
*   この場合，何も考えず，チュートリアル序盤のように，
*   型 G, V, E, W などを用いて，どうぞ書きなぐって下さい．
*
*
* (2) 使用頻度の高いコアアルゴリズムの実装時
*
*   テンプレートを用いて関数を記述して下さい．
*   既に (1) の状態で記述してある関数が有れば，書き換えて下さい．
*
*   チュートリアル中盤で紹介したグラフへのアクセスは重み無し・有りのグラフで共通です．
*   また，データ構造等もなるべく共通した処理を実現するように設計してあります．
*   運が良ければ，中身を何もいじらずに抽象化できることでしょう．
*/


//
// 2. 実例：連結性判定
//
namespace {
// 以下は，テンプレートを用いて連結性判定を行うコードの例です．
// 中身に注目して下さい．以前に紹介した幅優先探索のコードとほぼ全く同じです．
// しかし，このコードは，重み無しグラフ，重み有りグラフの両方で動作します．
template<typename GraphType>
bool my_super_connectivity_test(const GraphType &g) {
  if (g.num_vertices() == 0) return true;

  std::queue<V> que;
  std::vector<bool> vis(g.num_vertices());
  que.push(0);
  vis[0] = true;
  while (!que.empty()) {
    V v = que.front();
    que.pop();
    for (V tv : g.neighbors(v)) {
      if (vis[tv]) continue;
      que.push(tv);
      vis[tv] = true;
    }
  }
  return std::find(vis.begin(), vis.end(), false) == vis.end();
}
}  // namespace

int main() {
  // 重みなしグラフで上の関数を使ってみます．
  {
    // ちなみに |G| は |unweighted_graph| の alias です．
    unweighted_graph g({
      {0, 1}, {1, 2}, {3, 4}
    });
    cout << my_super_connectivity_test(g) << endl;
  }
  // 重み有りグラフで上の関数を使ってみます．
  {
    weighted_graph<double> g({
      {0, {1, 2.5}},
      {1, {2, 0.3}},
    });
    cout << my_super_connectivity_test(g) << endl;
  }

  return 0;
}
