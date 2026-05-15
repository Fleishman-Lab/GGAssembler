use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
mod colourful_dijkstra_impl;
mod scored;

use ordered_float::NotNan;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, VecDeque};

use crate::colourful_dijkstra_impl::{
    dijkstra, dijkstra_with_node_colors_and_bounds_options, Predecessor, SearchOptions,
};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::graphmap::DiGraphMap;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn node_count_from_edges(edges: &[(usize, usize, f32)], start: usize, goal: usize) -> usize {
    edges
        .iter()
        .fold(start.max(goal), |max_node, (src, dst, _)| {
            max_node.max(*src).max(*dst)
        })
        + 1
}

fn reconstruct_usize_path(
    predecessor: &Predecessor<usize, usize>,
    goal: usize,
    node: usize,
    node_color: usize,
) -> Vec<usize> {
    let mut path: Vec<usize> = vec![];
    if node == goal {
        let mut current = node;
        let mut current_color = node_color;
        path.push(current);
        while let Some((pred, pred_color)) = predecessor.get(&(current, current_color)) {
            path.push(*pred);
            current = *pred;
            current_color = *pred_color;
        }
    }

    path.reverse();
    path
}

fn reconstruct_node_index_path(
    predecessor: &Predecessor<NodeIndex, usize>,
    goal: NodeIndex,
    node: NodeIndex,
    node_color: usize,
) -> Vec<usize> {
    let mut path: Vec<usize> = vec![];
    if node == goal {
        let mut current = node;
        let mut current_color = node_color;
        path.push(current.index());
        while let Some((pred, pred_color)) = predecessor.get(&(current, current_color)) {
            path.push(pred.index());
            current = *pred;
            current_color = *pred_color;
        }
    }

    path.reverse();
    path
}

fn validate_gate_colors(
    node_count: usize,
    gate_color_ids_by_node: &[Vec<usize>],
    all_color_count: usize,
) -> PyResult<()> {
    if gate_color_ids_by_node.len() != node_count {
        return Err(PyValueError::new_err(
            "gate_color_ids_by_node length must match node count",
        ));
    }
    for color_ids in gate_color_ids_by_node {
        for color_id in color_ids {
            if *color_id >= all_color_count {
                return Err(PyValueError::new_err(
                    "gate color id is outside all_color_count",
                ));
            }
        }
    }
    Ok(())
}

fn recolor_node_masks(
    gate_color_ids_by_node: &[Vec<usize>],
    all_color_count: usize,
    no_colors: usize,
    rng: &mut impl Rng,
) -> Vec<usize> {
    let mut color_mapping = Vec::with_capacity(all_color_count);
    for _ in 0..all_color_count {
        let bit = rng.gen_range(0, no_colors);
        color_mapping.push(1usize << bit);
    }

    let mut node_masks = Vec::with_capacity(gate_color_ids_by_node.len());
    for color_ids in gate_color_ids_by_node {
        let mut mask = 0usize;
        for color_id in color_ids {
            mask |= color_mapping[*color_id];
        }
        node_masks.push(mask);
    }
    node_masks
}

fn validate_no_colors(no_colors: usize) -> PyResult<()> {
    if no_colors == 0 || no_colors > usize::BITS as usize {
        return Err(PyValueError::new_err(
            "no_colors must be between 1 and the number of bits in usize",
        ));
    }
    Ok(())
}

fn validate_find_many_args(
    min_gates: usize,
    max_gates: usize,
    no_colors: Option<usize>,
) -> PyResult<()> {
    if min_gates > max_gates {
        return Err(PyValueError::new_err("min_gates must be <= max_gates"));
    }
    let max_no_colors = max_gates
        .checked_add(1)
        .ok_or_else(|| PyValueError::new_err("max_gates is too large"))?;
    if max_no_colors > usize::BITS as usize {
        return Err(PyValueError::new_err(
            "max_gates + 1 must fit in the number of bits in usize",
        ));
    }
    if let Some(no_colors) = no_colors {
        validate_no_colors(no_colors)?;
    }
    Ok(())
}

fn search_options(use_a_star: bool, use_dominance: bool) -> SearchOptions {
    SearchOptions {
        use_a_star,
        use_dominance,
    }
}

fn reverse_bounds_from_edges(
    edges: &[(usize, usize, f32)],
    node_count: usize,
    goal: usize,
) -> (Vec<Option<f32>>, Vec<Option<usize>>) {
    let mut reverse_edges = vec![Vec::new(); node_count];
    for (src, dst, weight) in edges {
        if *src < node_count && *dst < node_count {
            reverse_edges[*dst].push((*src, *weight));
        }
    }

    let mut min_cost_to_goal = vec![None; node_count];
    let mut cost_heap = BinaryHeap::new();
    if goal < node_count {
        min_cost_to_goal[goal] = Some(0.0);
        cost_heap.push((Reverse(NotNan::new(0.0).unwrap()), goal));
    }
    while let Some((Reverse(cost), node)) = cost_heap.pop() {
        if min_cost_to_goal[node] != Some(cost.into_inner()) {
            continue;
        }
        for (previous, edge_cost) in &reverse_edges[node] {
            let next_cost = cost + *edge_cost;
            if min_cost_to_goal[*previous]
                .map(|known_cost| next_cost < known_cost)
                .unwrap_or(true)
            {
                min_cost_to_goal[*previous] = Some(next_cost);
                cost_heap.push((Reverse(NotNan::new(next_cost).unwrap()), *previous));
            }
        }
    }

    let mut min_hops_to_goal = vec![None; node_count];
    let mut hop_queue = VecDeque::new();
    if goal < node_count {
        min_hops_to_goal[goal] = Some(0);
        hop_queue.push_back(goal);
    }
    while let Some(node) = hop_queue.pop_front() {
        let next_hops = min_hops_to_goal[node].unwrap() + 1;
        for (previous, _edge_cost) in &reverse_edges[node] {
            if min_hops_to_goal[*previous].is_none() {
                min_hops_to_goal[*previous] = Some(next_hops);
                hop_queue.push_back(*previous);
            }
        }
    }

    (min_cost_to_goal, min_hops_to_goal)
}

fn dense_dag_adjacency_from_edges(
    edges: &[(usize, usize, f32)],
    node_count: usize,
    start: usize,
    goal: usize,
) -> PyResult<(Vec<Vec<(usize, f32)>>, bool)> {
    if start >= node_count || goal >= node_count {
        return Err(PyValueError::new_err(
            "start and goal must be valid node indices",
        ));
    }

    let mut adj = vec![Vec::new(); node_count];
    let mut is_dense_topological_dag = true;
    for (src, dst, weight) in edges {
        if *src >= node_count || *dst >= node_count {
            return Err(PyValueError::new_err("edge endpoint is outside node_count"));
        }
        if !weight.is_finite() || *weight < 0.0 {
            return Err(PyValueError::new_err(
                "edge weights must be finite and non-negative",
            ));
        }
        if src >= dst {
            is_dense_topological_dag = false;
        }
        adj[*src].push((*dst, *weight));
    }

    Ok((adj, is_dense_topological_dag))
}

/*
#[pyclass]
struct ColourfulPathFinder {
    graph: DiGraphMap<usize, f32>,
    start: usize,
    goal: usize,
    gate_color_ids_by_node: Option<Vec<Vec<usize>>>,
    all_color_count: usize,
    min_cost_to_goal: Vec<Option<f32>>,
    min_hops_to_goal: Vec<Option<usize>>,
}

impl ColourfulPathFinder {
    fn lower_bound_for_node(&self, node: usize) -> Option<(f32, usize)> {
        match (
            self.min_cost_to_goal.get(node),
            self.min_hops_to_goal.get(node),
        ) {
            (Some(Some(min_cost)), Some(Some(min_hops))) => Some((*min_cost, *min_hops)),
            _ => None,
        }
    }

    fn find_shortest_path_from_map(
        &self,
        node_color_map: HashMap<usize, usize>,
        limit: Option<usize>,
        options: SearchOptions,
    ) -> Vec<usize> {
        let (predecessor, node, node_color) = dijkstra_with_node_colors_and_bounds_options(
            &self.graph,
            self.start,
            self.goal,
            |e| *e.2,
            |node| node_color_map[node],
            |node| self.lower_bound_for_node(*node),
            limit,
            options,
        );
        reconstruct_usize_path(&predecessor, self.goal, node, node_color)
    }

    fn find_shortest_path_from_node_colors(
        &self,
        node_colors: &[usize],
        limit: Option<usize>,
        options: SearchOptions,
    ) -> Vec<usize> {
        let (predecessor, node, node_color) = dijkstra_with_node_colors_and_bounds_options(
            &self.graph,
            self.start,
            self.goal,
            |e| *e.2,
            |node| node_colors[*node],
            |node| self.lower_bound_for_node(*node),
            limit,
            options,
        );
        reconstruct_usize_path(&predecessor, self.goal, node, node_color)
    }

    fn find_many_inner(
        &self,
        gate_color_ids_by_node: &[Vec<usize>],
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        rng_seed: u64,
        no_colors: Option<usize>,
        options: SearchOptions,
    ) -> Vec<(usize, usize, Vec<usize>)> {
        let mut rng = StdRng::seed_from_u64(rng_seed);
        let mut results = Vec::new();

        for max_gates_for_run in min_gates..=max_gates {
            let no_colors = no_colors.unwrap_or(max_gates_for_run + 1);
            for retry_index in 0..retries {
                let node_colors = recolor_node_masks(
                    gate_color_ids_by_node,
                    self.all_color_count,
                    no_colors,
                    &mut rng,
                );
                let path = self.find_shortest_path_from_node_colors(
                    &node_colors,
                    Some(max_gates_for_run),
                    options,
                );
                if !path.is_empty() {
                    results.push((max_gates_for_run, retry_index, path));
                }
            }
        }
        results
    }
}

#[pymethods]
impl ColourfulPathFinder {
    #[new]
    fn new(edges: Vec<(usize, usize, f32)>, start: usize, goal: usize) -> Self {
        let node_count = node_count_from_edges(&edges, start, goal);
        let (min_cost_to_goal, min_hops_to_goal) =
            reverse_bounds_from_edges(&edges, node_count, goal);
        let graph = DiGraphMap::<usize, f32>::from_edges(&edges);
        Self {
            graph,
            start,
            goal,
            gate_color_ids_by_node: None,
            all_color_count: 0,
            min_cost_to_goal,
            min_hops_to_goal,
        }
    }

    #[staticmethod]
    fn with_gate_colors(
        edges: Vec<(usize, usize, f32)>,
        start: usize,
        goal: usize,
        node_count: usize,
        gate_color_ids_by_node: Vec<Vec<usize>>,
        all_color_count: usize,
    ) -> PyResult<Self> {
        validate_gate_colors(node_count, &gate_color_ids_by_node, all_color_count)?;
        let (min_cost_to_goal, min_hops_to_goal) =
            reverse_bounds_from_edges(&edges, node_count, goal);
        let graph = DiGraphMap::<usize, f32>::from_edges(&edges);
        Ok(Self {
            graph,
            start,
            goal,
            gate_color_ids_by_node: Some(gate_color_ids_by_node),
            all_color_count,
            min_cost_to_goal,
            min_hops_to_goal,
        })
    }

    #[pyo3(signature = (node_color_map, limit=None, use_a_star=true, use_dominance=true))]
    fn find_shortest_path(
        &self,
        py: Python<'_>,
        node_color_map: HashMap<usize, usize>,
        limit: Option<usize>,
        use_a_star: bool,
        use_dominance: bool,
    ) -> Vec<usize> {
        let options = search_options(use_a_star, use_dominance);
        py.allow_threads(move || self.find_shortest_path_from_map(node_color_map, limit, options))
    }

    #[pyo3(signature = (node_colors, limit=None, use_a_star=true, use_dominance=true))]
    fn find_shortest_path_with_node_colors(
        &self,
        py: Python<'_>,
        node_colors: Vec<usize>,
        limit: Option<usize>,
        use_a_star: bool,
        use_dominance: bool,
    ) -> Vec<usize> {
        let options = search_options(use_a_star, use_dominance);
        py.allow_threads(move || {
            self.find_shortest_path_from_node_colors(&node_colors, limit, options)
        })
    }

    #[pyo3(signature = (min_gates, max_gates, retries, seed=None, use_a_star=true, use_dominance=true, no_colors=None))]
    fn find_many(
        &self,
        py: Python<'_>,
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        seed: Option<u64>,
        use_a_star: bool,
        use_dominance: bool,
        no_colors: Option<usize>,
    ) -> PyResult<Vec<(usize, usize, Vec<usize>)>> {
        validate_find_many_args(min_gates, max_gates, no_colors)?;
        let gate_color_ids_by_node = self
            .gate_color_ids_by_node
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("finder was not constructed with gate colors"))?;
        let rng_seed = seed.unwrap_or_else(|| rand::thread_rng().gen());
        let options = search_options(use_a_star, use_dominance);
        Ok(py.allow_threads(move || {
            self.find_many_inner(
                gate_color_ids_by_node,
                min_gates,
                max_gates,
                retries,
                rng_seed,
                no_colors,
                options,
            )
        }))
    }
}
*/

#[pyclass]
struct ColourfulPathFinderDiGraph {
    graph: DiGraph<(), f32>,
    dag_adj: Vec<Vec<(usize, f32)>>,
    is_dense_topological_dag: bool,
    start: NodeIndex,
    goal: NodeIndex,
    gate_color_ids_by_node: Option<Vec<Vec<usize>>>,
    all_color_count: usize,
    min_cost_to_goal: Vec<Option<f32>>,
    min_hops_to_goal: Vec<Option<usize>>,
}

#[derive(Clone, Copy, Debug)]
struct DagLabel {
    mask: usize,
    cost: f32,
    depth: usize,
}

type DagPredecessor = HashMap<(usize, usize), (usize, usize)>;

fn dag_label_dominates(existing: DagLabel, candidate: DagLabel) -> bool {
    (existing.mask & candidate.mask) == existing.mask
        && existing.cost <= candidate.cost
        && existing.depth <= candidate.depth
}

fn insert_dag_nondominated_label(labels: &mut Vec<DagLabel>, candidate: DagLabel) -> bool {
    if labels
        .iter()
        .any(|existing| dag_label_dominates(*existing, candidate))
    {
        return false;
    }
    labels.retain(|existing| !dag_label_dominates(candidate, *existing));
    labels.push(candidate);
    true
}

impl ColourfulPathFinderDiGraph {
    fn lower_bound_for_node(&self, node: NodeIndex) -> Option<(f32, usize)> {
        let node = node.index();
        match (
            self.min_cost_to_goal.get(node),
            self.min_hops_to_goal.get(node),
        ) {
            (Some(Some(min_cost)), Some(Some(min_hops))) => Some((*min_cost, *min_hops)),
            _ => None,
        }
    }

    fn find_shortest_path_from_map(
        &self,
        node_color_map: HashMap<usize, usize>,
        limit: Option<usize>,
        options: SearchOptions,
    ) -> Vec<usize> {
        let indexed_node_color_map: HashMap<NodeIndex, usize> = node_color_map
            .into_iter()
            .map(|(node, color)| (NodeIndex::new(node), color))
            .collect();
        let (predecessor, node, node_color) = dijkstra_with_node_colors_and_bounds_options(
            &self.graph,
            self.start,
            self.goal,
            |e| *e.weight(),
            |node| indexed_node_color_map[node],
            |node| self.lower_bound_for_node(*node),
            limit,
            options,
        );
        reconstruct_node_index_path(&predecessor, self.goal, node, node_color)
    }

    fn find_shortest_path_from_node_colors(
        &self,
        node_colors: &[usize],
        limit: Option<usize>,
        options: SearchOptions,
    ) -> Vec<usize> {
        let (predecessor, node, node_color) = dijkstra_with_node_colors_and_bounds_options(
            &self.graph,
            self.start,
            self.goal,
            |e| *e.weight(),
            |node| node_colors[node.index()],
            |node| self.lower_bound_for_node(*node),
            limit,
            options,
        );
        reconstruct_node_index_path(&predecessor, self.goal, node, node_color)
    }

    fn find_shortest_path_dag_dp_from_node_colors(
        &self,
        node_colors: &[usize],
        limit: Option<usize>,
        use_dominance: bool,
        use_color_count_bound: bool,
    ) -> PyResult<Vec<usize>> {
        if !self.is_dense_topological_dag {
            return Err(PyValueError::new_err(
                "DAG DP requires dense node order to be topological: every edge must have src < dst",
            ));
        }
        let node_count = self.dag_adj.len();
        if node_colors.len() != node_count {
            return Err(PyValueError::new_err(
                "node_colors length must match node count",
            ));
        }

        let start = self.start.index();
        let goal = self.goal.index();
        let limit = limit.unwrap_or(0);
        let start_mask = node_colors[start];
        let mut scores: Vec<HashMap<usize, f32>> =
            (0..node_count).map(|_| HashMap::new()).collect();
        let mut labels: Vec<Vec<DagLabel>> = (0..node_count).map(|_| Vec::new()).collect();
        let mut predecessor: DagPredecessor = HashMap::new();

        scores[start].insert(start_mask, 0.0);
        labels[start].push(DagLabel {
            mask: start_mask,
            cost: 0.0,
            depth: 0,
        });

        for u in 0..node_count {
            if scores[u].is_empty() {
                continue;
            }
            let current_labels = labels[u].clone();
            for label in current_labels {
                if scores[u].get(&label.mask) != Some(&label.cost) {
                    continue;
                }
                if use_color_count_bound && limit > 0 {
                    if let Some(min_hops) = self.min_hops_to_goal[u] {
                        if label.mask.count_ones() as usize + min_hops > limit + 1 {
                            continue;
                        }
                    }
                }
                if limit > 0 && label.depth + 1 > limit {
                    continue;
                }
                for (v, weight) in &self.dag_adj[u] {
                    let v_color = node_colors[*v];
                    if (label.mask & v_color) != 0 {
                        continue;
                    }
                    let next_depth = label.depth + 1;
                    let next_mask = label.mask | v_color;
                    if use_color_count_bound && limit > 0 {
                        if let Some(min_hops) = self.min_hops_to_goal[*v] {
                            if next_mask.count_ones() as usize + min_hops > limit + 1 {
                                continue;
                            }
                        }
                    }
                    let next_cost = label.cost + *weight;
                    let should_update = scores[*v]
                        .get(&next_mask)
                        .map(|known| next_cost < *known)
                        .unwrap_or(true);
                    if !should_update {
                        continue;
                    }

                    let next_label = DagLabel {
                        mask: next_mask,
                        cost: next_cost,
                        depth: next_depth,
                    };
                    if use_dominance {
                        if !insert_dag_nondominated_label(&mut labels[*v], next_label) {
                            continue;
                        }
                    } else {
                        labels[*v].push(next_label);
                    }

                    scores[*v].insert(next_mask, next_cost);
                    predecessor.insert((*v, next_mask), (u, label.mask));
                }
            }
        }

        let Some((&best_mask, _best_cost)) =
            scores[goal].iter().min_by(|(_, a), (_, b)| a.total_cmp(b))
        else {
            return Ok(Vec::new());
        };

        let mut path = vec![goal];
        let mut current = goal;
        let mut current_mask = best_mask;
        while let Some((previous, previous_mask)) = predecessor.get(&(current, current_mask)) {
            path.push(*previous);
            current = *previous;
            current_mask = *previous_mask;
        }
        path.reverse();
        if path.first() == Some(&start) {
            Ok(path)
        } else {
            Ok(Vec::new())
        }
    }

    fn find_many_inner(
        &self,
        gate_color_ids_by_node: &[Vec<usize>],
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        rng_seed: u64,
        no_colors: Option<usize>,
        options: SearchOptions,
    ) -> Vec<(usize, usize, Vec<usize>)> {
        let mut rng = StdRng::seed_from_u64(rng_seed);
        let mut results = Vec::new();

        for max_gates_for_run in min_gates..=max_gates {
            let no_colors = no_colors.unwrap_or(max_gates_for_run + 1);
            for retry_index in 0..retries {
                let node_colors = recolor_node_masks(
                    gate_color_ids_by_node,
                    self.all_color_count,
                    no_colors,
                    &mut rng,
                );
                let path = self.find_shortest_path_from_node_colors(
                    &node_colors,
                    Some(max_gates_for_run),
                    options,
                );
                if !path.is_empty() {
                    results.push((max_gates_for_run, retry_index, path));
                }
            }
        }
        results
    }

    fn find_many_dag_dp_inner(
        &self,
        gate_color_ids_by_node: &[Vec<usize>],
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        rng_seed: u64,
        no_colors: Option<usize>,
        use_dominance: bool,
        use_color_count_bound: bool,
    ) -> PyResult<Vec<(usize, usize, Vec<usize>)>> {
        let mut rng = StdRng::seed_from_u64(rng_seed);
        let mut results = Vec::new();

        for max_gates_for_run in min_gates..=max_gates {
            let no_colors = no_colors.unwrap_or(max_gates_for_run + 1);
            for retry_index in 0..retries {
                let node_colors = recolor_node_masks(
                    gate_color_ids_by_node,
                    self.all_color_count,
                    no_colors,
                    &mut rng,
                );
                let path = self.find_shortest_path_dag_dp_from_node_colors(
                    &node_colors,
                    Some(max_gates_for_run),
                    use_dominance,
                    use_color_count_bound,
                )?;
                if !path.is_empty() {
                    results.push((max_gates_for_run, retry_index, path));
                }
            }
        }
        Ok(results)
    }
}

#[pymethods]
impl ColourfulPathFinderDiGraph {
    #[new]
    fn new(
        edges: Vec<(usize, usize, f32)>,
        start: usize,
        goal: usize,
        node_count: Option<usize>,
    ) -> PyResult<Self> {
        let final_node_count =
            node_count.unwrap_or_else(|| node_count_from_edges(&edges, start, goal));
        let (dag_adj, is_dense_topological_dag) =
            dense_dag_adjacency_from_edges(&edges, final_node_count, start, goal)?;
        let (min_cost_to_goal, min_hops_to_goal) =
            reverse_bounds_from_edges(&edges, final_node_count, goal);
        let indexed_edges = edges
            .iter()
            .map(|(src, dst, weight)| (NodeIndex::new(*src), NodeIndex::new(*dst), *weight));
        let mut graph = DiGraph::<(), f32>::from_edges(indexed_edges);
        while graph.node_count() < final_node_count {
            graph.add_node(());
        }
        Ok(Self {
            graph,
            dag_adj,
            is_dense_topological_dag,
            start: NodeIndex::new(start),
            goal: NodeIndex::new(goal),
            gate_color_ids_by_node: None,
            all_color_count: 0,
            min_cost_to_goal,
            min_hops_to_goal,
        })
    }

    #[staticmethod]
    fn with_gate_colors(
        edges: Vec<(usize, usize, f32)>,
        start: usize,
        goal: usize,
        node_count: usize,
        gate_color_ids_by_node: Vec<Vec<usize>>,
        all_color_count: usize,
    ) -> PyResult<Self> {
        validate_gate_colors(node_count, &gate_color_ids_by_node, all_color_count)?;
        let (dag_adj, is_dense_topological_dag) =
            dense_dag_adjacency_from_edges(&edges, node_count, start, goal)?;
        let (min_cost_to_goal, min_hops_to_goal) =
            reverse_bounds_from_edges(&edges, node_count, goal);
        let indexed_edges = edges
            .iter()
            .map(|(src, dst, weight)| (NodeIndex::new(*src), NodeIndex::new(*dst), *weight));
        let mut graph = DiGraph::<(), f32>::from_edges(indexed_edges);
        while graph.node_count() < node_count {
            graph.add_node(());
        }
        Ok(Self {
            graph,
            dag_adj,
            is_dense_topological_dag,
            start: NodeIndex::new(start),
            goal: NodeIndex::new(goal),
            gate_color_ids_by_node: Some(gate_color_ids_by_node),
            all_color_count,
            min_cost_to_goal,
            min_hops_to_goal,
        })
    }

    #[pyo3(signature = (node_color_map, limit=None, use_a_star=true, use_dominance=true))]
    fn find_shortest_path(
        &self,
        py: Python<'_>,
        node_color_map: HashMap<usize, usize>,
        limit: Option<usize>,
        use_a_star: bool,
        use_dominance: bool,
    ) -> Vec<usize> {
        let options = search_options(use_a_star, use_dominance);
        py.allow_threads(move || self.find_shortest_path_from_map(node_color_map, limit, options))
    }

    #[pyo3(signature = (node_colors, limit=None, use_a_star=true, use_dominance=true))]
    fn find_shortest_path_with_node_colors(
        &self,
        py: Python<'_>,
        node_colors: Vec<usize>,
        limit: Option<usize>,
        use_a_star: bool,
        use_dominance: bool,
    ) -> Vec<usize> {
        let options = search_options(use_a_star, use_dominance);
        py.allow_threads(move || {
            self.find_shortest_path_from_node_colors(&node_colors, limit, options)
        })
    }

    #[pyo3(signature = (node_colors, limit=None, use_dominance=true, use_color_count_bound=false))]
    fn find_shortest_path_dag_dp_with_node_colors(
        &self,
        py: Python<'_>,
        node_colors: Vec<usize>,
        limit: Option<usize>,
        use_dominance: bool,
        use_color_count_bound: bool,
    ) -> PyResult<Vec<usize>> {
        py.allow_threads(move || {
            self.find_shortest_path_dag_dp_from_node_colors(
                &node_colors,
                limit,
                use_dominance,
                use_color_count_bound,
            )
        })
    }

    #[pyo3(signature = (min_gates, max_gates, retries, seed=None, use_a_star=true, use_dominance=true, no_colors=None))]
    fn find_many(
        &self,
        py: Python<'_>,
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        seed: Option<u64>,
        use_a_star: bool,
        use_dominance: bool,
        no_colors: Option<usize>,
    ) -> PyResult<Vec<(usize, usize, Vec<usize>)>> {
        validate_find_many_args(min_gates, max_gates, no_colors)?;
        let gate_color_ids_by_node = self
            .gate_color_ids_by_node
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("finder was not constructed with gate colors"))?;
        let rng_seed = seed.unwrap_or_else(|| rand::thread_rng().gen());
        let options = search_options(use_a_star, use_dominance);
        Ok(py.allow_threads(move || {
            self.find_many_inner(
                gate_color_ids_by_node,
                min_gates,
                max_gates,
                retries,
                rng_seed,
                no_colors,
                options,
            )
        }))
    }

    #[pyo3(signature = (min_gates, max_gates, retries, seed=None, use_dominance=true, no_colors=None, use_color_count_bound=false))]
    fn find_many_dag_dp(
        &self,
        py: Python<'_>,
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        seed: Option<u64>,
        use_dominance: bool,
        no_colors: Option<usize>,
        use_color_count_bound: bool,
    ) -> PyResult<Vec<(usize, usize, Vec<usize>)>> {
        validate_find_many_args(min_gates, max_gates, no_colors)?;
        let gate_color_ids_by_node = self
            .gate_color_ids_by_node
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("finder was not constructed with gate colors"))?;
        let rng_seed = seed.unwrap_or_else(|| rand::thread_rng().gen());
        py.allow_threads(move || {
            self.find_many_dag_dp_inner(
                gate_color_ids_by_node,
                min_gates,
                max_gates,
                retries,
                rng_seed,
                no_colors,
                use_dominance,
                use_color_count_bound,
            )
        })
    }
}

/// Finds shortest path between source and target and an optional path length limit.
/// ----------
///
/// :param List[Tuple[(uint, uint uint)...] edges: A list of triplets [(a,b, w) ...] representing a directed edge from a to b with weight w.
///
/// :param uint start:  The source node
///
/// :param uint goal: The target node
///
/// :param Dict[uint, uint] node_color_map: A mapping of node to it's colour
///
/// :param uint limit: The path length limit
///
/// :rtype: List[uint]
/// :return: A list of nodes of the shortest path if one is found or an empty list if no path was found.

#[pyfunction]
fn colourful_shortest_path(
    py: Python<'_>,
    edges: Vec<(usize, usize, f32)>,
    start: usize,
    goal: usize,
    node_color_map: HashMap<usize, usize>,
    limit: Option<usize>,
) -> Vec<usize> {
    py.allow_threads(move || {
        let graph = DiGraphMap::<usize, f32>::from_edges(&edges);
        let (predecessor, node, node_color) = dijkstra(
            &graph,
            start.into(),
            goal.into(),
            |e| *e.2,
            node_color_map,
            limit,
        );
        reconstruct_usize_path(&predecessor, goal, node, node_color)
    })
}

/// A Python module implemented in Rust.
#[pymodule]
fn colourful_dijkstra(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(colourful_shortest_path, m)?)?;
    m.add_class::<ColourfulPathFinderDiGraph>()?;
    Ok(())
}
