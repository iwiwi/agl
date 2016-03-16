Box-Cover Algorithms
========================

# Building

```
./waf configure
./waf
```

# Running

```
./bin/box_cover -type=tsv -graph=/data/graph_edges.tsv \
 -method=sketch -alpha=1.0 -least_coverage=1.0 -sketch_k=128 -multipass=10000 -rad_min=1 -rad_max=30 \
 -random_seed=114514 
```

## Options
|Options          |                                                |Type   |Default|
|:----------------|:-----------------------------------------------|:-----:|:----:|
|-type            |Graph file type (auto, tsv, agl, gen) |string | "auto"|
|-graph           |Input graph                                     |string | "-"   |
|-method          |Using method (see below)           |string |"sketch"|
|-alpha           |Index size limit to use MEMB (alpha * N * k)        |double |1.0    |
|-least_coverage  |Least coverage.                   |double |1.0|
|-sketch_k        |sketch k                                        |int32|128|
|-multipass       |Maximum number of multi-pass      |int32|1000000000|
|-rad_analytical  |Use analytical diameters for radius?            | bool|false|
|-rad_min         |Minimum radius.                                  |int32 |  1|
|-rad_max         |Maximum radius.                                  |int32 |100000000|
|-random_seed     |Random seed.                                    |int64|922337203685477583|

## Available Methods

|Name|Method|
|:--|:--|
|sketch|Sketch (Akiba et al. 2016)|
|memb|MEMB (Song et al. 2007)|
|coloring|Greedy Coloring (Song et al. 2007)|
|cbb|CBB (Song et al. 2007)|
|burning|Burning (Schneider et al. 2012)|
|analytical|Optimal solutions of Box-Cover of (u,v)-flower or SHM-model|


# Supported Formats

`bin/box_cover` works with TSV file, AGL file and User-Generated graph.

## TSV files 
* By using the option `-type=tsv`, you can specify .tsv file as the graph file with the option `-graph=...`.
* In a .tsv file, each line should contain two integers describing an edge (see below).
* Vertices should be described by integers starting from zero.

### TSV file example
```
0 1
0 2
1 2
```

## AGL files
__math is GOD.

## User-Generated graph
* By using the option `-type=gen`, you can generate and use famous network models.

### Barabási–Albert model
### (u, v)-flower
### Song-Havlin-Makse model

# Unit Test
Execute `bin/test` to run tests.

# LICENSE

This software is released under the MIT License, see LICENSE.txt.
