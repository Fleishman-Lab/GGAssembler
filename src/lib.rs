use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;
mod colourful_dijkstra_impl;
mod scored;

use std::collections::HashMap;

use crate::colourful_dijkstra_impl::dijkstra;

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
    edges: Vec<(usize, usize, usize)>,
    start: usize,
    goal: usize,
    node_color_map: HashMap<usize, usize>,
    limit: Option<usize>,
) -> Vec<usize> {
    let node_count = edges
        .iter()
        .flat_map(|(source, target, _)| [*source, *target])
        .chain(node_color_map.keys().copied())
        .chain([start, goal])
        .max()
        .unwrap_or(0)
        + 1;
    let mut graph = vec![Vec::new(); node_count];
    for (source, target, weight) in edges {
        graph[source].push((target, weight));
    }
    let mut node_colors = vec![0; node_count];
    for (node, color) in node_color_map {
        node_colors[node] = color;
    }
    let (predecessor, (node, node_color, node_depth)) =
        dijkstra(&graph, start, goal, &node_colors, limit);
    // predecessor
    let mut path: Vec<usize> = vec![];
    if node == goal {
        let mut current = node;
        let mut current_color = node_color;
        let mut current_depth = node_depth;
        path.push(current);
        while let Some((pred, pred_color, pred_depth)) =
            predecessor.get(&(current, current_color, current_depth))
        {
            path.push(*pred);
            current = *pred;
            current_color = *pred_color;
            current_depth = *pred_depth;
        }
    }

    path.reverse();
    path
}

/// A Python module implemented in Rust.
#[pymodule]
fn colourful_dijkstra(m: &Bound<'_, PyModule>) -> PyResult<()> {
    //m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(colourful_shortest_path, m)?)?;
    Ok(())
}
