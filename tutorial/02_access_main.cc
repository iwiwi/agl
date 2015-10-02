#include "agl.h"
using namespace std;
using namespace agl;

int main() {
  G g(generate_grid(3, 3));
  pretty_print(g);

  //
  // 1. 辺へのアクセス
  //
  {
    // 辺の型は |E| です（実際は int）．
    // 関数 |edge| を使ってアクセスします．
    // 以下では頂点 1 の 0 番目の辺へアクセスしています．
    const E &e = g.edge(1, 0);

    // 辺の情報はフリー関数 |to| と |weight| を用いてアクセスします．
    // 重み無しグラフでは |weight| は常に 1 です．
    // 重みの型は W（実際は int）です．
    {
      V v = to(e);
      W w = weight(e);
      cout << v << ", " << w << endl;
    }

    // 次数は degree で取得でき，以下のように辺は走査できます．
    for (size_t i = 0; i < g.degree(1); ++i) {
      const E &e = g.edge(1, i);
      cout << to(e) << ", " << weight(e) << endl;
    }

    // ただし，関数 |edges| を用いると range-based for で辺を走査できます．
    // 以下では頂点 1 から出ている辺を全て舐めています．
    for (const auto &e : g.edges(1)) {
      cout << e << endl;
    }

    // 連結成分分解や幅優先探索など，重みを用いないアルゴリズムの記述では，
    // |edge|, |edges| の代わりに |neighbor|, |neighbors| を必ず用いて下さい．
    V n1 = g.neighbor(1, 0);  // |to(g.edge(1, 0))| と同値
    assert(n1 == to(g.edge(1, 0)));

    for (V n : g.neighbors(1)) {
      cout << n << endl;
    }
  }

  //
  // 2. アルゴリズム実装の例：幅優先探索
  //
  {
    std::queue<V> que;
    std::vector<bool> vis(g.num_vertices());

    // 頂点 0 から幅優先探索する
    que.push(0);
    vis[0] = true;

    while (!que.empty()) {
      V v = que.front();
      que.pop();
      for (V u : g.neighbors(v)) {
        if (vis[u]) continue;
        cout << "BFS: " << u << endl;  // 頂点 u に到達
        que.push(u);
        vis[u] = true;
      }
    }
  }

  //
  // 3. AGL 内の幅優先探索・Dijkstra のアルゴリズム
  //
  {
    // 頂点 0 から全頂点への距離を求めます．
    vector<int> distances = single_source_distance(g, 0);
    cout << distances << endl;

    // 頂点 0 から近い頂点を順に訪問します．
    // ラムダ式が |false| を返した瞬間に探索を終了します．
    visit_by_distance(g, 0, [](V v, W w) {
      cout << v << " " << w << endl;
      return true;
    });

    // ただし，|visit_by_distance| 関数は初期化に O(|V|) 時間を要します．
    // 狭い範囲の探索を繰り返す場合，再初期化を必要な部分に限定する
    // |visitor_by_distance| クラスを使うことで計算量を抑えられます．
    visitor_by_distance<G> vis(g);
    vis.visit(0, [](V v, W w) {
      cout << v << " " << w << endl;
      return true;
    });
  }

  return 0;
}
