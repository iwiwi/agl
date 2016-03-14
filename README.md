# Box-Cover Algorithms

## Building

```
./waf configure
./waf
```

## Running

### Options
|Options||Type|Default|
|:--|:--|:--|--:|
|-random_seed |Random seed. |int64|922337203685477583|
|-alpha       | Index size limit to use MEMB (alpha*n*k)|double|1.0|
|-least_coverage |coverage|double|1.0|
|-method|using method|string|"sketch"|
|-multipass |Number of multi-pass|int32|1000000000|
|-rad_analytical|Using analytical diameters for rads| bool|false|
|-rad_max |maximum radius|int32 |100000000|
|-rad_min |minimum radius|int32 |  1|
|-sketch_k|sketch k|int32|128|
|-force_undirected|Automatically add reverse edges?| boolfalse|
|-graph |input graph|string| "-"|
|-type  |Graph file type (auto, tsv, agl, built_in, gen) |string | "auto"|

## LICENSE

HOGE