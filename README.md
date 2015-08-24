AGL: Akihabara Graph Library (仮)
=================================

![ミナリンスキー](http://sifblog.net/wp-content/uploads/2015/01/1422485254-c9d4989e5ee5afcf888fa69bcd69cafd.jpg)

AGL は C++ のグラフ処理ライブラリです．ユニークな特徴として，高速アルゴリズムの研究と大規模データの解析をスムーズに接続するために，以下の 2 点を同時に達成することを目標としています．

1. **アルゴリズムの迅速なプロトタイピング** - - -アルゴリズムの研究者が自身のアイディアを試したり論文に向けて性能測定を行ったりする際の労力を最小化するように設計されています．
2. **最先端アルゴリズムを用いた大規模ネットワーク解析** - - - AGL を用いて開発されたアルゴリズムがそのまま取り込まれており，発表されたばかりの最先端の手法を用いて大規模グラフを解析することができます *(予定)* ．



## アルゴリズム研究者向けガイドライン


### プロトタイピング・ベンチマークの段階

AGL を private fork して自分のリポジトリを作り，そこで作業して下さい．rebase のために master にはコミットせず，ブランチを作って作業するようにだけして下さい． http://qiita.com/potix2/items/cddd7dd9cde6a9c6dde6

その中では何をやってもらっても構いません．好きにコードを書き殴って下さい．思う存分測定して下さい．

### AGL に取り込む段階

AGL を汚染されると困るので，以下の条件を満たすようにコードを整理して下さい．

* 他人が include し得るヘッダー内では using namespace をしない 
* インターフェース部では命名規則を守る
* グローバル変数は使わない
* テストをつける（そもそもプロトタイピングの時につけてると思いますが一応）
* 自分だけの namespace を区切る，namespace の中では何でもし放題（命名規則も好きにして下さい）

.cc の中や自分しか include しない header の中では using namespace していいです．#define もしてもいいです．

完成したら pull request を送って下さい．



## 利用条件

現在のところ，秋葉との共同研究にして頂ける場合に限ります．まだ未成熟のプロジェクトであるため，目が行き届きフィードバックを迅速に反映できるような範囲内での利用に制限したいと考えています．
