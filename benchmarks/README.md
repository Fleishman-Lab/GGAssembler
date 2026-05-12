# Colorful Path Timing

Use `bench_colorful_path.py` before and after colorful shortest-path changes.
It builds one synthetic directed graph, then repeats the same retry and gate-limit
loop shape used by `example/colab_oligos_design.ipynb`.

Baseline command:

```sh
python benchmarks/bench_colorful_path.py --runs 5 --retries 20 --seed 1
```

Larger run:

```sh
python benchmarks/bench_colorful_path.py \
  --nodes 800 \
  --edge-probability 0.02 \
  --runs 5 \
  --retries 100 \
  --seed 1
```

Compare `median_seconds` before and after a change. Keep `nodes`,
`edge-probability`, `colors`, gate limits, retries, and seed fixed when comparing
two implementations. The benchmark uses `--seed` for synthetic graph
construction and reseeds NumPy as `seed + run_index` before each timed run, so
separate invocations with the same arguments execute the same sequence of
random recolorings.

Real-data graph setup smoke run:

```sh
python benchmarks/bench_colorful_real_data.py \
  --runs 1 \
  --retries 0 \
  --min-gates 26 \
  --max-gates 26 \
  --seed 1
```

The real-data benchmark uses `example/wt_dna.fasta`,
`example/input.resfile`, and the checked-in `example/deg_table.csv`. It avoids
regenerating the degenerate table because that path can fetch codon usage data
from the network. Use nonzero `--retries` when you want to time the actual
colorful search. The matching real-data pytest always verifies graph
construction from the checked-in inputs; the expensive colorful search test is
separately opt-in:

```sh
python -m pytest dawdlib/dijkstra/tests/test_colorful_real_data.py -q
GGASSEMBLER_RUN_REAL_DATA_SEARCH=1 python -m pytest dawdlib/dijkstra/tests/test_colorful_real_data.py -q
```
