use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::collections::{BinaryHeap, HashMap};

use std::hash::Hash;

use petgraph::visit::{EdgeRef, IntoEdges, VisitMap, Visitable};
use petgraph::algo::Measure;
use crate::scored::MinScored;

/// \[Generic\] Dijkstra's shortest path algorithm.
///
/// Compute the length of the shortest path from `start` to every reachable
/// node.
///
/// The graph should be `Visitable` and implement `IntoEdges`. The function
/// `edge_cost` should return the cost for a particular edge, which is used
/// to compute path costs. Edge costs must be non-negative.
///
/// If `goal` is not `None`, then the algorithm terminates once the `goal` node's
/// cost is calculated.
///
/// Returns a `HashMap` that maps `NodeId` to path cost.
/// # Example
/// ```rust
/// use petgraph::Graph;
/// use petgraph::algo::dijkstra;
/// use petgraph::prelude::*;
/// use std::collections::HashMap;
///
/// let mut graph : Graph<(),(),Directed>= Graph::new();
/// let a = graph.add_node(()); // node with no weight
/// let b = graph.add_node(());
/// let c = graph.add_node(());
/// let d = graph.add_node(());
/// let e = graph.add_node(());
/// let f = graph.add_node(());
/// let g = graph.add_node(());
/// let h = graph.add_node(());
/// // z will be in another connected component
/// let z = graph.add_node(());
///
/// graph.extend_with_edges(&[
///     (a, b),
///     (b, c),
///     (c, d),
///     (d, a),
///     (e, f),
///     (b, e),
///     (f, g),
///     (g, h),
///     (h, e)
/// ]);
/// // a ----> b ----> e ----> f
/// // ^       |       ^       |
/// // |       v       |       v
/// // d <---- c       h <---- g
///
/// let expected_res: HashMap<NodeIndex, usize> = [
///      (a, 3),
///      (b, 0),
///      (c, 1),
///      (d, 2),
///      (e, 1),
///      (f, 2),
///      (g, 3),
///      (h, 4)
///     ].iter().cloned().collect();
/// let res = dijkstra(&graph,b,None, |_| 1);
/// assert_eq!(res, expected_res);
/// // z is not inside res because there is not path from b to z.
/// ```
pub fn dijkstra<G, F, K, C>(
    graph: G,
    start: G::NodeId,
    goal: G::NodeId,
    mut edge_cost: F,
    node_color_map: HashMap<G::NodeId, C>,
    limit: Option<usize>
) -> (HashMap<(G::NodeId,C), (G::NodeId,C)>, G::NodeId, C)
    where
        G: IntoEdges + Visitable,
        G::NodeId: Eq + Hash,
        F: FnMut(G::EdgeRef) -> K,
        K: Measure + Copy,
        C: num_traits::PrimInt + num_traits::Unsigned + Hash + std::fmt::Display,
{
    let limit = limit.unwrap_or(0);
    let mut visited = graph.visit_map();
    let mut scores: HashMap<(G::NodeId, C), K> = HashMap::new();
    let mut predecessor: HashMap<(G::NodeId,C), (G::NodeId,C)> = HashMap::new();
    let mut visit_next: BinaryHeap<MinScored<K, (C, usize, G::NodeId)>> = BinaryHeap::new();
    let zero_score = K::default();
    let start_color = node_color_map[&start];
    scores.insert((start, start_color), zero_score);
    visit_next.push(MinScored(zero_score, (start_color, 0 as usize, start)));
    while let Some(MinScored(node_score, (node_color, node_depth, node))) = visit_next.pop() {
        if visited.is_visited(&node) {
            continue;
        }
        if goal == node {
            return (predecessor, node.clone(), node_color)
        }
        if limit > 0 && node_depth + 1 > limit {
            continue;
        }
        for edge in graph.edges(node) {
            let next = edge.target();
            let next_color = node_color_map[&next];
            if visited.is_visited(&next) {
                continue;
            }
            if (next_color & node_color) != C::min_value() {
                continue;
            }
            let next_score = node_score + edge_cost(edge);
            let next_color = node_color | next_color;
            match scores.entry((next, next_color) ) {
                Occupied(ent) => {
                    if next_score < *ent.get() {
                        *ent.into_mut() = next_score;
                        visit_next.push(MinScored(next_score, (next_color, node_depth + 1, next)));
                        predecessor.insert((next.clone(), next_color), (node.clone(), node_color));
                    }
                }
                Vacant(ent) => {
                    ent.insert(next_score);
                    visit_next.push(MinScored(next_score, (next_color, node_depth + 1, next)));
                    predecessor.insert((next.clone(), next_color), (node.clone(), node_color));
                }
            }
        }
        visited.visit(node);
    }
    (predecessor, start.clone(), start_color)
}
