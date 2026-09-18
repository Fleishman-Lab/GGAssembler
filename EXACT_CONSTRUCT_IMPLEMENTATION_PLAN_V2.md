# Exact Construct Implementation Plan — Review Draft

## Status

Approved for implementation by the user.

## Goal

Complete the existing exact-construct prototype so GGAssembler can design one exact protein variant per well while minimizing the acquisition cost of the distinct physical DNA parts shared across those variants.

The implementation will preserve the existing Golden Gate graph formulation:

- node: candidate Golden Gate junction;
- edge: candidate fragment interval;
- edge weight: acquisition cost of the distinct physical parts required for that interval;
- path: one segmentation shared by every requested variant;
- path validity: the existing Golden Gate compatibility constraints.

The implementation will extend the existing exact-construct code. It will not redesign the ordinary combinatorial-library workflow.

## Confirmed scope

### Included

- Exact variant DNA generation using the existing fixed-codon behavior.
- Explicit sequence-invariance checking for candidate overhangs.
- Per-interval deduplication of identical variant payloads.
- Independent sourcing of WT-equivalent and non-WT payloads within the same interval.
- Individual-oligo pricing through a caller-supplied length-cost function.
- Infeasibility reported by the pricing function rather than an arbitrary fallback cost.
- Optional comparison of WT PCR cost with direct WT-fragment synthesis.
- Integer eurocent edge weights.
- A unique physical-part registry.
- A variant-to-part bill of materials.
- Summary statistics derived from those tables.
- A documented notebook example for excluding user-specified overhangs.
- Unit, integration, and MthUPO end-to-end tests.

### Explicitly excluded

- Synonymous-codon optimization.
- Joint codon and segmentation optimization.
- Any change to the Rust colorful shortest-path implementation.
- Any attempt to fix or replace color-mask state handling.
- The `performance_improvments` branch.
- DAG dynamic programming.
- Optimality guarantees or near-optimal/Pareto search.
- Oligopool pricing or retrieval-handle implementation.
- Retrieval handles in individual-oligo mode.
- Monetary penalties for Echo transfers or fill-in chemistry.
- Changes to the existing non-exact cost model.

These topics may be revisited after the exact individual-oligo workflow is complete.

## Existing machinery to reuse

### Exact-variant input and sequence generation

`dawdlib/exact_construct.py` already provides:

- resfile parsing;
- variant CSV parsing;
- mutation records;
- deterministic construction of exact DNA variants from selected codons;
- restriction-site validation of generated variant DNA.

The existing codon-selection behavior will be preserved. No synonymous alternatives will be generated or optimized.

### Golden Gate graph

`dawdlib/golden_gate/graph_maker.py` already provides:

- candidate gate generation;
- graph construction;
- edge generation;
- default overhang efficiency filtering;
- local fidelity checks;
- a custom edge-weight callback;
- a custom node-validity callback through `build_custom_graph()`.

### Golden Gate compatibility

The implementation will retain:

- `GGData.filter_self_binding_gates()` for the standard acceptable-overhang set;
- pairwise fidelity checks in `GraphMaker`;
- the existing restriction graph and path validation;
- existing reaction/fidelity reporting.

### Sequence slicing

The interval slicing convention already used by `get_wt()` and `_gate_oligo_slice()` will remain the single definition of the payload between two gates.

### Pricing callback

The current prototype already passes `segment_length_cost` into the exact edge-weight function. That mechanism will be retained and corrected rather than replaced by a class hierarchy.

## Pricing-function contract

The exact individual-oligo workflow will accept a callable with this conceptual type:

```python
Callable[[int], Optional[int]]
```

Its contract is:

- input: complete physical ordered oligo length in nucleotides;
- output: price in integer eurocents;
- `None`: this length is unavailable and therefore physically infeasible.

The function receives the complete ordered length, including the fixed prefix and suffix actually present on that particular oligo. Terminal intervals may contain only one or neither addition.

The default individual-oligo function will be constructed with the existing `make_length_tier_cost()` helper using:

```python
individual_oligo_cost = make_length_tier_cost(
    [
        (13, 60, 13),  # cents per nt
        (61, 80, 17),
    ]
)
```

`make_length_tier_cost()` will change as follows:

- tier prices are integer cents per nucleotide;
- if multiple configured tiers match, use the cheapest matching price;
- if no tier matches and no explicit fallback was supplied, return `None`;
- never silently return `float(length)`;
- an explicitly supplied fallback remains allowed because it is then a deliberate caller-defined sourcing rule.

The following boundary behavior will be tested:

- 12 nt: infeasible;
- 13 nt: 169 cents;
- 60 nt: 780 cents;
- 61 nt: 1,037 cents;
- 80 nt: 1,360 cents;
- 81 nt: infeasible.

## WT PCR sourcing

The default WT PCR estimate will be configurable but initially use:

```text
2 primers × 20 nt × 13 cents/nt = 520 cents
```

Twenty nucleotides is an intentionally approximate experimental default. Experimentalists may configure another primer length or price when appropriate.

WT-equivalent payloads will support two modes:

1. `allow_wt_direct_synthesis=False`
   - use WT PCR for eligible internal WT-equivalent parts;
   - do not compare against direct synthesis.

2. `allow_wt_direct_synthesis=True`
   - price the complete directly ordered WT oligo with the supplied pricing function;
   - compare it with the configured PCR cost;
   - choose the cheaper feasible source.

Both modes will be covered by tests. The initial default will be `False`, preserving the existing preference for WT-template amplification while making the alternative explicit.

The universal fill-in primer will not add a per-part optimization cost.

## Terminal intervals

WT-equivalent material in an interval touching the source or target pseudogate will retain zero acquisition cost. A non-WT terminal payload still requires a valid individual-oligo quote, using only the fixed additions actually present on that terminal oligo.

Terminal WT-equivalent material will still be represented in the part registry and variant BOM so the assembly description remains complete, but its source type will indicate that it is terminal/free and its cost will be zero.

No attempt will be made in this task to model the case-specific cost of incorporating terminal sequence into a plasmid or another segment.

## Exact interval planning

A single helper will plan every candidate edge. Conceptually:

```python
plan_exact_interval(
    gate1,
    gate2,
    wt_dna,
    variant_dnas,
    prefix,
    suffix,
    oligo_length_cost,
    primer_length,
    primer_price_per_nt_cents,
    allow_wt_direct_synthesis,
) -> Optional[ExactIntervalPlan]
```

The exact name and lightweight record type may follow the local naming style, but the behavior will be:

1. Determine the interval payload slice using the existing gate slicing convention.
2. Extract that payload from every requested variant.
3. Group variant IDs by identical payload sequence.
4. Create one physical-part record per distinct payload in that interval.
5. For WT-equivalent payloads in terminal intervals, assign zero cost. Price non-WT terminal payloads as individual oligos.
6. For an internal WT-equivalent payload:
   - quote WT PCR;
   - optionally quote direct synthesis;
   - choose the configured/cheapest feasible source.
7. For an internal non-WT payload:
   - build the complete ordered sequence from the additions actually used on that edge;
   - pass its complete length to `oligo_length_cost`;
   - reject the interval when the function returns `None`.
8. Sum each unique physical part once to obtain `cost_cents`.
9. Retain the payload-to-variant mapping as edge metadata.

The same interval plan will drive both optimization and final output. Output generation will not independently reconstruct or re-price the interval.

## Edge feasibility

The exact path will not blindly reuse the legacy variable/constant classification.

An exact edge is feasible only if:

- its gates do not overlap;
- its internal overhang pair satisfies the existing minimum-fidelity rule;
- each required internal non-WT oligo receives a valid quote;
- a directly synthesized WT oligo, when selected, receives a valid quote;
- PCR-specific fragment constraints are satisfied when WT PCR is selected.

For individual oligos, the pricing callback defines the supported physical length range. The legacy `min_oligo_length`, `max_oligo_length`, and `oligo_addition` values will not be used as a second, inconsistent length model for exact individual-oligo sourcing.

`min_const_oligo_length` will apply only where it has a defined PCR-fragment meaning. It will not disqualify a valid synthesized part.

## Candidate overhang invariance

Exact mode will use a direct sequence-based check rather than the union of mutation positions.

For a candidate gate starting at `gate.idx` with length `gate.gatelength`, every requested variant must have exactly the same sequence:

```python
variant_dna[gate.idx : gate.idx + gate.gatelength] == gate.bps
```

This will be combined with the existing acceptable-overhang filter. It remains correct without relying on codon-position indexing and does not introduce synonymous recoding.

## Custom overhang filtering

`make_exact_variant_graph()` will accept one optional callable:

```python
overhang_filter: Optional[Callable[[str], bool]] = None
```

The standard GGAssembler filtering will run first. The custom callable may only remove otherwise-valid overhangs; it will not bypass efficiency, GC, sequence-invariance, or fidelity checks.

The exact notebook will include a documented example for excluding overhangs already used in a collaborator's plasmid assembly:

```python
from dawdlib.golden_gate.gate_data import reverse_complement

plasmid_overhangs = {"AATG", "GCTT"}
excluded_overhangs = plasmid_overhangs | {
    reverse_complement(overhang) for overhang in plasmid_overhangs
}

def overhang_is_available(overhang):
    return overhang not in excluded_overhangs

graph, source, target = make_exact_variant_graph(
    ...,
    overhang_filter=overhang_is_available,
)
```

The example will state that excluding reverse complements is intentional for this use case, while callers remain free to provide different filtering rules.

## Graph edge metadata and weights

Every accepted exact edge will store:

- `weight`: total acquisition cost in integer eurocents;
- its interval identifier/gates;
- its planned unique physical parts;
- the mapping from each requested variant to its part within the interval.

The existing path-search machinery will consume only `weight`. No path-search code or numeric representation will be changed in this task.

The graph-building change will be limited to allowing exact edge construction to attach this metadata while leaving ordinary graph construction behavior unchanged.

## Output model

The prototype `exact_variant_oligo_table()` currently creates one row for every variant occurrence and duplicates shared physical sequences. It will be replaced by an exact-output function returning at least two linked DataFrames and a summary.

Because the exact prototype is not used outside its tests or the untracked example notebook, no compatibility wrapper is planned unless implementation reveals another caller.

### Physical part registry

One row per distinct physical stock within an assembly interval.

Planned fields:

- `part_id`;
- `interval_id`;
- `gate1`;
- `gate2`;
- `payload_sequence`;
- `ordered_sequence`;
- `source_type`:
  - `wt_pcr`;
  - `individual_oligo`;
  - `terminal_free`;
- `is_wt`;
- `ordered_length`;
- `cost_cents`;
- `cost_eur` as a presentation-only derived column;
- `usage_count`;
- primer metadata where applicable.

Part identity will include the assembly role. Two coincidentally identical sequences used in different intervals will not be merged unless they are genuinely interchangeable physical stocks.

### Variant bill of materials

One row per requested variant and assembly interval.

Planned fields:

- `variant_id`;
- `interval_id`;
- `part_id`;
- `gate1`;
- `gate2`;
- `is_wt`;
- `source_type`.

Variants sharing an interval payload will reference the same `part_id`.

This table is the future input for Echo transfer/pick-list generation, but Echo file generation itself is outside this task.

## Summary statistics

Summary values will be calculated from the selected path, part registry, and BOM rather than maintained separately.

The report will include:

- number of requested variants;
- number of unique physical parts;
- number of unique synthesized non-WT oligos;
- number of directly synthesized WT oligos;
- number of WT-PCR fragments;
- number of terminal/free parts;
- total ordered oligo nucleotides;
- total primer nucleotides;
- total acquisition cost in cents and euros;
- total fragment usages/transfers;
- average non-WT parts per variant;
- average non-common parts per variant;
- total Golden Gate fragments per construct;
- number of internal Golden Gate junctions;
- selected junction sequences and positions;
- existing fidelity metrics available for the selected path.

Definitions:

- non-WT part: `is_wt == False`;
- non-common part: a part whose `usage_count` is less than the number of requested variants;
- transfer count: the number of BOM rows, including zero-cost terminal parts.

Neither average nor transfer count will affect the monetary objective.

## Input validation

Before graph construction, exact mode will validate:

- at least one requested variant exists;
- variant IDs are unique;
- every variant DNA has the same length as WT DNA;
- sequences contain only supported unambiguous DNA bases after normalization;
- generated variant translation matches the requested amino-acid substitutions;
- configured restriction-enzyme checks still pass;
- the pricing callback returns either a non-negative integer number of cents or `None`.

Invalid pricing callback results such as negative values, floats, strings, or booleans will raise a clear error rather than enter the graph.

## File-by-file implementation plan

### `dawdlib/golden_gate/graph_maker.py`

- Correct `make_length_tier_cost()` so unmatched lengths return `None` unless the caller explicitly supplied a fallback.
- Use integer-cent tier values.
- Preserve the caller-supplied `segment_length_cost` mechanism.
- Replace the prototype exact weight calculation with the shared interval planner.
- Remove retrieval-handle length from individual exact mode.
- Add sequence-invariant node validation for exact variants.
- Add the optional `overhang_filter` argument to `make_exact_variant_graph()`.
- Attach interval plan metadata to exact graph edges.
- Leave `create_default_weight_func()` and the ordinary graph path unchanged.

### `dawdlib/exact_construct.py`

- Retain resfile/variant parsing and fixed-codon exact DNA generation.
- Do not add synonymous optimization.
- Add the minimal record(s) needed to represent an interval plan and physical part.
- Add WT-PCR and optional WT-direct-synthesis selection.
- Generate the physical-part registry and BOM from graph edge metadata on the selected path.
- Generate summary statistics.
- Replace the duplicating prototype exact oligo table API with the linked output model.

### `dawdlib/golden_gate/tests/test_exact_variant_weights.py`

- Replace tests that expect fallback pricing or retrieval handles.
- Add pricing boundary and infeasibility tests.
- Add complete-ordered-length tests.
- Add unique-fragment reuse, WT reuse, mixed interval, and optional WT synthesis tests.
- Add terminal/free interval tests.
- Add custom overhang-filter tests.

### `dawdlib/integration_tests/test_exact_construct.py`

- Retain parsing and fixed-codon tests.
- Replace the duplicate-row output expectation.
- Add linked registry/BOM tests.
- Add summary-statistic tests, including both requested averages.
- Add a toy global segmentation test using the existing path search unchanged.
- Add explicit invariant-junction tests.

### `example/MthUPO/`

- Add the existing example inputs to version control so tests do not depend on untracked local files.
- Do not modify the scientific input data.

### `example/colab_exact_oligos_design.ipynb`

- Replace the current copied degenerate-library workflow with the exact-variant workflow.
- Load the supplied MthUPO WT sequence, resfile, and exact variant CSV.
- Show the default individual-oligo pricing function.
- Show both WT sourcing modes.
- Show the custom plasmid-overhang exclusion callback.
- Build the exact graph and use the existing path search unchanged.
- Display selected gates, cost, fidelity, physical-part registry, BOM, and summary.
- Export the registry, BOM, selected path, and summary.

### Files intentionally unchanged

- `src/lib.rs`;
- `src/colourful_dijkstra_impl.rs`;
- `src/scored.rs`;
- the performance branch;
- ordinary combinatorial-library output code;
- oligopool/retrieval code, because none exists yet.

## Test matrix

### Pricing

- Every default tier boundary from 12 through 81 nt described above.
- Multiple matching tiers choose the cheapest result.
- Explicit fallback works only when supplied.
- Invalid callback results are rejected.
- Complete ordered length, not payload length, controls feasibility.

### Deduplication and sourcing

- Many variants sharing one payload create and charge one part.
- WT used by many variants creates one WT part.
- Mixed WT plus two mutant payloads creates exactly three physical parts.
- Usage frequency changes `usage_count`, not acquisition cost.
- Individual mode adds no retrieval handles.
- WT PCR costs 520 cents with default primer settings.
- WT direct synthesis disabled forces PCR.
- WT direct synthesis enabled chooses synthesis when cheaper.
- WT direct synthesis enabled keeps PCR when synthesis is unavailable or more expensive.
- Terminal/free WT parts cost zero but remain in registry and BOM; terminal mutant parts are priced normally.

### Graph and Golden Gate behavior

- A gate identical in every exact variant remains eligible.
- A gate differing in any exact variant is excluded.
- The custom overhang filter removes an otherwise-valid overhang.
- The custom filter cannot restore an overhang rejected by standard checks.
- Existing pairwise and global compatibility validation remains active.

### Output

- Registry has one row per unique interval part.
- BOM reuses a part ID across variants.
- Every variant has exactly one BOM row per path interval.
- Usage counts equal BOM references.
- Total cost equals the sum of registry costs.
- Total ordered and primer nucleotides match the registry.
- Average non-WT and average non-common values are independently correct.

### End-to-end

- A small hand-calculable exact example chooses the expected segmentation under the supplied pricing function.
- MthUPO inputs parse into 150 unique requested variants.
- At least one complete MthUPO exact path is generated with the unchanged search machinery.
- The registry and BOM reconstruct every requested exact variant.
- Printed total cost agrees with an independent calculation from the registry.
- The existing non-exact test suite remains unchanged and passes.

## Validation procedure

After implementation approval:

1. Run the exact pricing/unit tests.
2. Run exact integration tests.
3. Run the existing full Python test suite.
4. Build the existing Rust extension without modifying its source.
5. Execute the small hand-calculable workflow.
6. Execute the MthUPO notebook workflow.
7. Manually verify one mixed interval containing WT and mutant parts.
8. Confirm `git diff` contains no Rust/search or synonymous-optimization changes.
9. Report files changed, tests run, results, assumptions, and any observed limitations.

## Implementation order

1. Correct and test the callable length-pricing helper.
2. Add and test exact interval planning and sourcing.
3. Use interval plans for exact graph feasibility, weight, and metadata.
4. Add invariant and custom overhang filtering.
5. Build the registry, BOM, and summary from selected-edge metadata.
6. Add the small end-to-end test.
7. Track and run the MthUPO fixtures.
8. Update and execute the exact notebook.
9. Run regression tests and produce the implementation report.

Each stage should leave a runnable test before moving to the next.

## Acceptance criteria

The implementation is complete when:

- each distinct internal physical part is charged once per assembly role;
- WT and mutant payloads on the same interval can use different source types;
- individual-oligo cost and feasibility come only from the supplied pricing function applied to complete ordered length;
- unsupported default lengths are infeasible;
- individual mode has no retrieval handles;
- terminal WT-equivalent material remains free but is represented in output;
- graph weights are integer eurocents;
- the registry and BOM are linked and internally consistent;
- both non-WT and non-common averages are reported;
- the notebook demonstrates custom overhang exclusion;
- the MthUPO workflow runs end to end;
- ordinary GGAssembler behavior remains unchanged;
- no Rust/search or synonymous-optimization work appears in the diff.

## Deferred work

- Oligopool/retrieval pricing and handle design.
- Exact or improved colorful-path search.
- Near-optimal and Pareto path reporting.
- Synonymous-codon optimization.
- Case-specific terminal-fragment costing.
- Echo transfer-file generation.
