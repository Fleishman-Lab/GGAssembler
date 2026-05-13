use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;
mod colourful_dijkstra_impl;
mod scored;

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, VecDeque};

use crate::colourful_dijkstra_impl::{dijkstra, dijkstra_with_node_colors_and_bounds};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::graphmap::DiGraphMap;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn node_count_from_edges(edges: &[(usize, usize, usize)], start: usize, goal: usize) -> usize {
    edges
        .iter()
        .fold(start.max(goal), |max_node, (src, dst, _)| {
            max_node.max(*src).max(*dst)
        })
        + 1
}

fn reconstruct_usize_path(
    predecessor: &HashMap<(usize, usize), (usize, usize)>,
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
    predecessor: &HashMap<(NodeIndex, usize), (NodeIndex, usize)>,
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

fn validate_find_many_args(min_gates: usize, max_gates: usize) -> PyResult<()> {
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
    Ok(())
}

fn reverse_bounds_from_edges(
    edges: &[(usize, usize, usize)],
    node_count: usize,
    goal: usize,
) -> (Vec<Option<usize>>, Vec<Option<usize>>) {
    let mut reverse_edges = vec![Vec::new(); node_count];
    for (src, dst, weight) in edges {
        if *src < node_count && *dst < node_count {
            reverse_edges[*dst].push((*src, *weight));
        }
    }

    let mut min_cost_to_goal = vec![None; node_count];
    let mut cost_heap = BinaryHeap::new();
    if goal < node_count {
        min_cost_to_goal[goal] = Some(0);
        cost_heap.push((Reverse(0usize), goal));
    }
    while let Some((Reverse(cost), node)) = cost_heap.pop() {
        if min_cost_to_goal[node] != Some(cost) {
            continue;
        }
        for (previous, edge_cost) in &reverse_edges[node] {
            let next_cost = cost + *edge_cost;
            if min_cost_to_goal[*previous]
                .map(|known_cost| next_cost < known_cost)
                .unwrap_or(true)
            {
                min_cost_to_goal[*previous] = Some(next_cost);
                cost_heap.push((Reverse(next_cost), *previous));
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

#[pyclass]
struct ColourfulPathFinder {
    graph: DiGraphMap<usize, usize>,
    start: usize,
    goal: usize,
    gate_color_ids_by_node: Option<Vec<Vec<usize>>>,
    all_color_count: usize,
    min_cost_to_goal: Vec<Option<usize>>,
    min_hops_to_goal: Vec<Option<usize>>,
}

impl ColourfulPathFinder {
    fn lower_bound_for_node(&self, node: usize) -> Option<(usize, usize)> {
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
    ) -> Vec<usize> {
        let (predecessor, node, node_color): (
            HashMap<(usize, usize), (usize, usize)>,
            usize,
            usize,
        ) = dijkstra_with_node_colors_and_bounds(
            &self.graph,
            self.start,
            self.goal,
            |e| *e.2,
            |node| node_color_map[node],
            |node| self.lower_bound_for_node(*node),
            limit,
        );
        reconstruct_usize_path(&predecessor, self.goal, node, node_color)
    }

    fn find_shortest_path_from_node_colors(
        &self,
        node_colors: &[usize],
        limit: Option<usize>,
    ) -> Vec<usize> {
        let (predecessor, node, node_color): (
            HashMap<(usize, usize), (usize, usize)>,
            usize,
            usize,
        ) = dijkstra_with_node_colors_and_bounds(
            &self.graph,
            self.start,
            self.goal,
            |e| *e.2,
            |node| node_colors[*node],
            |node| self.lower_bound_for_node(*node),
            limit,
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
    ) -> Vec<(usize, usize, Vec<usize>)> {
        let mut rng = StdRng::seed_from_u64(rng_seed);
        let mut results = Vec::new();

        for max_gates_for_run in min_gates..=max_gates {
            let no_colors = max_gates_for_run + 1;
            for retry_index in 0..retries {
                let node_colors = recolor_node_masks(
                    gate_color_ids_by_node,
                    self.all_color_count,
                    no_colors,
                    &mut rng,
                );
                let path =
                    self.find_shortest_path_from_node_colors(&node_colors, Some(max_gates_for_run));
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
    fn new(edges: Vec<(usize, usize, usize)>, start: usize, goal: usize) -> Self {
        let node_count = node_count_from_edges(&edges, start, goal);
        let (min_cost_to_goal, min_hops_to_goal) =
            reverse_bounds_from_edges(&edges, node_count, goal);
        let graph = DiGraphMap::<usize, usize>::from_edges(&edges);
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
        edges: Vec<(usize, usize, usize)>,
        start: usize,
        goal: usize,
        node_count: usize,
        gate_color_ids_by_node: Vec<Vec<usize>>,
        all_color_count: usize,
    ) -> PyResult<Self> {
        validate_gate_colors(node_count, &gate_color_ids_by_node, all_color_count)?;
        let (min_cost_to_goal, min_hops_to_goal) =
            reverse_bounds_from_edges(&edges, node_count, goal);
        let graph = DiGraphMap::<usize, usize>::from_edges(&edges);
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

    fn find_shortest_path(
        &self,
        py: Python<'_>,
        node_color_map: HashMap<usize, usize>,
        limit: Option<usize>,
    ) -> Vec<usize> {
        py.allow_threads(move || self.find_shortest_path_from_map(node_color_map, limit))
    }

    fn find_shortest_path_with_node_colors(
        &self,
        py: Python<'_>,
        node_colors: Vec<usize>,
        limit: Option<usize>,
    ) -> Vec<usize> {
        py.allow_threads(move || self.find_shortest_path_from_node_colors(&node_colors, limit))
    }

    fn find_many(
        &self,
        py: Python<'_>,
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        seed: Option<u64>,
    ) -> PyResult<Vec<(usize, usize, Vec<usize>)>> {
        validate_find_many_args(min_gates, max_gates)?;
        let gate_color_ids_by_node = self
            .gate_color_ids_by_node
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("finder was not constructed with gate colors"))?;
        let rng_seed = seed.unwrap_or_else(|| rand::thread_rng().gen());
        Ok(py.allow_threads(move || {
            self.find_many_inner(
                gate_color_ids_by_node,
                min_gates,
                max_gates,
                retries,
                rng_seed,
            )
        }))
    }
}

#[pyclass]
struct ColourfulPathFinderDiGraph {
    graph: DiGraph<(), usize>,
    start: NodeIndex,
    goal: NodeIndex,
    gate_color_ids_by_node: Option<Vec<Vec<usize>>>,
    all_color_count: usize,
    min_cost_to_goal: Vec<Option<usize>>,
    min_hops_to_goal: Vec<Option<usize>>,
}

impl ColourfulPathFinderDiGraph {
    fn lower_bound_for_node(&self, node: NodeIndex) -> Option<(usize, usize)> {
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
    ) -> Vec<usize> {
        let indexed_node_color_map: HashMap<NodeIndex, usize> = node_color_map
            .into_iter()
            .map(|(node, color)| (NodeIndex::new(node), color))
            .collect();
        let (predecessor, node, node_color): (
            HashMap<(NodeIndex, usize), (NodeIndex, usize)>,
            NodeIndex,
            usize,
        ) = dijkstra_with_node_colors_and_bounds(
            &self.graph,
            self.start,
            self.goal,
            |e| *e.weight(),
            |node| indexed_node_color_map[node],
            |node| self.lower_bound_for_node(*node),
            limit,
        );
        reconstruct_node_index_path(&predecessor, self.goal, node, node_color)
    }

    fn find_shortest_path_from_node_colors(
        &self,
        node_colors: &[usize],
        limit: Option<usize>,
    ) -> Vec<usize> {
        let (predecessor, node, node_color): (
            HashMap<(NodeIndex, usize), (NodeIndex, usize)>,
            NodeIndex,
            usize,
        ) = dijkstra_with_node_colors_and_bounds(
            &self.graph,
            self.start,
            self.goal,
            |e| *e.weight(),
            |node| node_colors[node.index()],
            |node| self.lower_bound_for_node(*node),
            limit,
        );
        reconstruct_node_index_path(&predecessor, self.goal, node, node_color)
    }

    fn find_many_inner(
        &self,
        gate_color_ids_by_node: &[Vec<usize>],
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        rng_seed: u64,
    ) -> Vec<(usize, usize, Vec<usize>)> {
        let mut rng = StdRng::seed_from_u64(rng_seed);
        let mut results = Vec::new();

        for max_gates_for_run in min_gates..=max_gates {
            let no_colors = max_gates_for_run + 1;
            for retry_index in 0..retries {
                let node_colors = recolor_node_masks(
                    gate_color_ids_by_node,
                    self.all_color_count,
                    no_colors,
                    &mut rng,
                );
                let path =
                    self.find_shortest_path_from_node_colors(&node_colors, Some(max_gates_for_run));
                if !path.is_empty() {
                    results.push((max_gates_for_run, retry_index, path));
                }
            }
        }
        results
    }
}

#[pymethods]
impl ColourfulPathFinderDiGraph {
    #[new]
    fn new(
        edges: Vec<(usize, usize, usize)>,
        start: usize,
        goal: usize,
        node_count: Option<usize>,
    ) -> Self {
        let final_node_count =
            node_count.unwrap_or_else(|| node_count_from_edges(&edges, start, goal));
        let (min_cost_to_goal, min_hops_to_goal) =
            reverse_bounds_from_edges(&edges, final_node_count, goal);
        let indexed_edges = edges
            .iter()
            .map(|(src, dst, weight)| (NodeIndex::new(*src), NodeIndex::new(*dst), *weight));
        let mut graph = DiGraph::<(), usize>::from_edges(indexed_edges);
        while graph.node_count() < final_node_count {
            graph.add_node(());
        }
        Self {
            graph,
            start: NodeIndex::new(start),
            goal: NodeIndex::new(goal),
            gate_color_ids_by_node: None,
            all_color_count: 0,
            min_cost_to_goal,
            min_hops_to_goal,
        }
    }

    #[staticmethod]
    fn with_gate_colors(
        edges: Vec<(usize, usize, usize)>,
        start: usize,
        goal: usize,
        node_count: usize,
        gate_color_ids_by_node: Vec<Vec<usize>>,
        all_color_count: usize,
    ) -> PyResult<Self> {
        validate_gate_colors(node_count, &gate_color_ids_by_node, all_color_count)?;
        let (min_cost_to_goal, min_hops_to_goal) =
            reverse_bounds_from_edges(&edges, node_count, goal);
        let indexed_edges = edges
            .iter()
            .map(|(src, dst, weight)| (NodeIndex::new(*src), NodeIndex::new(*dst), *weight));
        let mut graph = DiGraph::<(), usize>::from_edges(indexed_edges);
        while graph.node_count() < node_count {
            graph.add_node(());
        }
        Ok(Self {
            graph,
            start: NodeIndex::new(start),
            goal: NodeIndex::new(goal),
            gate_color_ids_by_node: Some(gate_color_ids_by_node),
            all_color_count,
            min_cost_to_goal,
            min_hops_to_goal,
        })
    }

    fn find_shortest_path(
        &self,
        py: Python<'_>,
        node_color_map: HashMap<usize, usize>,
        limit: Option<usize>,
    ) -> Vec<usize> {
        py.allow_threads(move || self.find_shortest_path_from_map(node_color_map, limit))
    }

    fn find_shortest_path_with_node_colors(
        &self,
        py: Python<'_>,
        node_colors: Vec<usize>,
        limit: Option<usize>,
    ) -> Vec<usize> {
        py.allow_threads(move || self.find_shortest_path_from_node_colors(&node_colors, limit))
    }

    fn find_many(
        &self,
        py: Python<'_>,
        min_gates: usize,
        max_gates: usize,
        retries: usize,
        seed: Option<u64>,
    ) -> PyResult<Vec<(usize, usize, Vec<usize>)>> {
        validate_find_many_args(min_gates, max_gates)?;
        let gate_color_ids_by_node = self
            .gate_color_ids_by_node
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("finder was not constructed with gate colors"))?;
        let rng_seed = seed.unwrap_or_else(|| rand::thread_rng().gen());
        Ok(py.allow_threads(move || {
            self.find_many_inner(
                gate_color_ids_by_node,
                min_gates,
                max_gates,
                retries,
                rng_seed,
            )
        }))
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
    edges: Vec<(usize, usize, usize)>,
    start: usize,
    goal: usize,
    node_color_map: HashMap<usize, usize>,
    limit: Option<usize>,
) -> Vec<usize> {
    py.allow_threads(move || {
        let graph = DiGraphMap::<usize, usize>::from_edges(&edges);
        let (predecessor, node, node_color): (
            HashMap<(usize, usize), (usize, usize)>,
            usize,
            usize,
        ) = dijkstra(
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
    //m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(colourful_shortest_path, m)?)?;
    m.add_class::<ColourfulPathFinder>()?;
    m.add_class::<ColourfulPathFinderDiGraph>()?;
    Ok(())
}
