#include "agl.h"
using namespace std;
using namespace agl;

int main() {
  //
  // 1. グラフを作ってみる
  //
  {
    // AGL では |V| が頂点の型です (実際は int）．
    // グラフを作るためには，まず辺リストを用意します．
    vector<pair<V, V>> el = {
        {0, 1},
        {1, 2},
        {2, 3},
    };

    // AGL では |G| がグラフの型です．
    G g(el);

    // |pretty_print| で概要を標準エラー出力に書きます．
    pretty_print(g);
  }

  //
  // 2. 有向・無向グラフ
  //
  {
    // 上と同じグラフをもう一度作ってみます．
    // |unweighted_edge_list| は |vector<pair<V, V>>| の alias です．
    unweighted_edge_list el1 = {
        {0, 1},
        {1, 2},
        {2, 3},
    };

    // 先ほどの |pretty_print| で分かったとおり，このように作成すると，
    // グラフは有向グラフになります．
    G g1(el1);
    assert(is_adjacent(g1, 0, 1) == true);
    assert(is_adjacent(g1, 1, 0) == false);

    // AGL では，無向グラフは，逆辺の追加により有向グラフと同じ扱いをします．
    // 上のパスグラフを無向にしたい場合，以下のようになります．
    unweighted_edge_list el2 = {
        {0, 1}, {1, 0},
        {1, 2}, {2, 1},
        {2, 3}, {3, 2},
    };

    G g2(el2);
    pretty_print(g2);
    assert(is_adjacent(g2, 0, 1) == true);
    assert(is_adjacent(g2, 1, 0) == true);

    // 一方向のみの辺リストに逆辺を追加するには，|make_undirected| 関数を用います．
    unweighted_edge_list el3 = make_undirected(el2);
    G g3(el3);
    pretty_print(g3);  // g2 と全く同じ
  }

  //
  // 3. 組み込みグラフデータ・ジェネレータ
  //
  {
    // パスグラフの辺リストは |generate_path| 関数により生成できます．
    G g(generate_path(3));

    // 他にも色々なジェネレータがあります．
    // 詳しくは |src/graph/generator.h| を見て下さい．
    pretty_print(G(generate_erdos_renyi(10, 2)));
    pretty_print(G(generate_grid(3, 3)));

    // また，小規模な実データも埋め込んでありすぐ使えるようになっています．
    pretty_print(G(built_in_edge_list("karate_club")));
    pretty_print(G(built_in_edge_list("ca_grqc")));
  }

  //
  // 4. 読み込み
  //
  {
    // 以下のようにすることでファイルから読み込みを行います．
    // （ファイルが存在しない環境でエラーが起きないよう，コメントアウトしています）
    // unweighted_edge_list el = read_edge_list_tsv("hoge.tsv");
  }

  return 0;
}
