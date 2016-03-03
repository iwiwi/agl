# フラクタル次元絶対求めるマン

## Average Approximation Ratio
フラクタルなモデルのグラフの計算結果に対して、パラメータ k ごと、α ごとにそれぞれ分けて Average Approximation Ratio を計算してグラフを出力してくれる。ただし、モデルは解析解が埋め込みされているものしか計算できないので、予め埋め込んでおく必要がある。

<img src="average_approximation_ratio_with_a.png" width=320px>
<img src="average_approximation_ratio_with_k.png" width=320px>

```
./average-approximation-ratio.py [解析に使いたいJLOGファイルたち]
./average-approximation-ratio.py jlog/*
```
## フィッティング計算
各アルゴリズムの各グラフに対する結果について、exponentail と power-law な関数をそれぞれ非線形最小二乗法でフィッティングし、その二条誤差からどちらにより近いかを判定する。同じアルゴリズムかつ同じグラフの解析結果が複数含まれていた場合実行時間のより短かったほうが使用される。
```
./calc_residual.py [解析に使いたいJLOGファイルたち]
./calc_residual.py jlog/*
```

## jlog 統合
半径ごとに別々のジョブとして実行した場合、それぞれの半径について計算した結果がjlogとして出力されるが、combine_rads.pyによって統合できる。

```
./combine_rads.py [統合したjlogファイルたち]
./combine_rads.py jlog/sketch-hollywood-2011.agl-rad.*
```
