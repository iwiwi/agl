# フラクタル次元絶対求めるマン

## Average Approximation Ratio
フラクタルなモデルのグラフの計算結果に対して、パラメータ k ごと、α ごとにそれぞれ分けて Average Approximation Ratio を計算してグラフを出力してくれる。ただし、モデルは解析解が埋め込みされているものしか計算できないので、予め埋め込んでおく必要がある。

<img src="average_approximation_ratio_with_a.png" width=320px>
<img src="average_approximation_ratio_with_k.png" width=320px>

./average-approximation-ratio.py [解析に使いたいJLOGファイル]
```
./average-approximation-ratio.py jlog/*
```
