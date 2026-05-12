# Colorful Path Timing

Use `bench_colorful_path.py` before and after colorful shortest-path changes.
It builds one synthetic directed graph, then repeats the same retry and gate-limit
loop shape used by `example/colab_oligos_design.ipynb`.

Baseline command:

```sh
python benchmarks/bench_colorful_path.py --runs 5 --retries 20
```

Larger run:

```sh
python benchmarks/bench_colorful_path.py \
  --nodes 800 \
  --edge-probability 0.02 \
  --runs 5 \
  --retries 100
```

Compare `median_seconds` before and after a change. Keep `nodes`,
`edge-probability`, `colors`, gate limits, retries, and seed fixed when comparing
two implementations.
