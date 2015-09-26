//
// 1. 早く自分のアイディアが試したい！待ちきれない！
//
/*
 * そんな人の願いを，たった 2 行で叶えるのが easy_cui です．
 * #include "easy_cui.h" して，main 内で |easy_cui_init| してグラフを受け取って下さい．
 *
 * この 2 行を書いておけば，あなたは受け取ったグラフを用いてアルゴリズムを実装することに集中できます．
 * そして，これだけで，様々なグラフを用いた実験ができるようになります．
 *
 * (1) 標準入力からグラフを読み込み
 *   $ bin/03_easy_cui <<< "0 1"
 *
 * (2) ファイルからグラフを読み込み
 *   $ bin/03_easy_cui --graph=hoge.tsv
 *
 * (3) 組み込みのグラフを読み込み
 *   $ bin/03_easy_cui --type=built_in --graph=karate_club
 *   $ bin/03_easy_cui --type=built_in --graph=ca_grqc
 *
 * (4) ジェネレータを使ってグラフを生成
 *   $ bin/03_easy_cui --type=gen --graph="grid 3 3"
 *
 * (オプション) グラフを無向化したいときは --force_undirected を足す
 *   $ bin/03_easy_cui --type=gen --graph="grid 3 3" --force_undirected
 *
 * (オプション) ジェネレータの乱数シードを変えたいときは，--random_seed を足す
 *   $ bin/03_easy_cui --type=gen --graph="erdos_renyi 10 2" --random_seed
 *
 * 特に，組み込みデータセットを使うことで，
 * データセットを手元に落としてきていない場所でも簡単な実験ができます．
 *
 * 詳しくは |easy_cui.h| の中身を見てみて下さい．
 */

#include "easy_cui.h"

int main(int argc, char **argv) {
  G g = easy_cui_init(argc, argv);
  pretty_print(g);

  //
  // 2. デバッグが辛い！グラフの構造を描いてみないと分からない！
  //
  {
    // たった 1 行でグラフを描画できます．
    // ジェネレータを用いた小さなグラフの生成と併せて，効率的なデバッグができます．
    graphviz_draw_graph(g, "out.png");

    // 注意：graphviz がインストールされている必要があります．
    //   $ sudo apt-get install graphviz
  }

  return 0;
}
