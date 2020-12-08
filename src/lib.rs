use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
mod colourful_dijkstra;
mod scored;

use std::collections::HashMap;

use petgraph::graphmap::DiGraphMap;
use crate::colourful_dijkstra::dijkstra;

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
    edges: Vec<(usize,usize,usize)>,
    start: usize,
    goal: usize,
    node_color_map: HashMap<usize, usize>,
    limit: Option<usize>
) -> Vec<usize>
{
    let graph = DiGraphMap::<usize, usize>::from_edges(&edges);
    let (predecessor, node, node_color): (HashMap<(usize, usize), (usize,usize)>, usize, usize) = dijkstra(&graph, start.into(), goal.into(), |e| *e.2, node_color_map, limit);
    // predecessor
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

/// A Python module implemented in Rust.
#[pymodule]
fn colourful_dijkstra(_py: Python, m: &PyModule) -> PyResult<()> {
    //m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(colourful_shortest_path, m)?)?;
    Ok(())
}