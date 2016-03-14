# Box-Cover Algorithms

## Building

```
./waf configure
./waf
```

## Running

```
./bin/box_cover -force_undirected -type=tsv -graph=/data/graph_edges.tsv \
 -method=sketch -alpha=1.0 -least_coverage=1.0 -sketch_k=128 -multipass=10000 -rad_min=1 -rad_max=30 --random_seed=114514 
```

### Options
|Options||Type|Default|
|:----------------|:-----------------------------------------------|:-----:|:----:|
|-force_undirected|Automatically add reverse edges?                | bool  |false  |
|-type            |Graph file type (auto, tsv, agl, built_in, gen) |string | "auto"|
|-graph           |Input graph                                     |string | "-"   |
|-method          |Using method                                    |string |"sketch"|
|-alpha           |Index size limit to use MEMB (alpha*n*k)        |double |1.0    |
|-least_coverage  |Least coverage.                   |double |1.0|
|-sketch_k        |sketch k                                        |int32|128|
|-multipass       |Number of multi-pass                            |int32|1000000000|
|-rad_analytical  |Use analytical diameters for radius?            | bool|false|
|-rad_min         |Minimum radius.                                  |int32 |  1|
|-rad_max         |Maximum radius.                                  |int32 |100000000|
|-random_seed     |Random seed.                                    |int64|922337203685477583|

### Available Mtehods

|Name|Method|
|:--|:--|
|sketch|Sketch (Akiba et al. 2016)|
|memb|MEMB (Song et al. 2007)|
|coloring|Greedy Coloring (Song et al. 2007)|
|cbb|CBB (Song et al. 2007)|
|burning|Burning (Schneider et al. 2012)|
|analytical|Optimal solutions of Box-Cover of (u,v)-flower or SHM-model|

## LICENSE

HOGE
