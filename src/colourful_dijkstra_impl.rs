use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::collections::{BinaryHeap, HashMap};

use std::hash::Hash;

use crate::scored::MinScored;
use petgraph::algo::Measure;
use petgraph::visit::{EdgeRef, IntoEdges, Visitable};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Label<K, C> {
    mask: C,
    cost: K,
    depth: usize,
}

fn label_dominates<K, C>(existing: &Label<K, C>, candidate: &Label<K, C>) -> bool
where
    K: PartialOrd,
    C: num_traits::PrimInt,
{
    (existing.mask & candidate.mask) == existing.mask
        && existing.cost <= candidate.cost
        && existing.depth <= candidate.depth
}

fn insert_nondominated_label<K, C>(labels: &mut Vec<Label<K, C>>, candidate: Label<K, C>) -> bool
where
    K: Copy + PartialOrd,
    C: Copy + num_traits::PrimInt,
{
    if labels
        .iter()
        .any(|existing| label_dominates(existing, &candidate))
    {
        return false;
    }

    labels.retain(|existing| !label_dominates(&candidate, existing));
    labels.push(candidate);
    true
}

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
    edge_cost: F,
    node_color_map: HashMap<G::NodeId, C>,
    limit: Option<usize>,
) -> (HashMap<(G::NodeId, C), (G::NodeId, C)>, G::NodeId, C)
where
    G: IntoEdges + Visitable,
    G::NodeId: Eq + Hash,
    F: FnMut(G::EdgeRef) -> K,
    K: Measure + Copy,
    C: num_traits::PrimInt + num_traits::Unsigned + Hash + std::fmt::Display,
{
    dijkstra_with_node_colors(
        graph,
        start,
        goal,
        edge_cost,
        |node| node_color_map[node],
        limit,
    )
}

pub fn dijkstra_with_node_colors<G, F, K, C, ColorFn>(
    graph: G,
    start: G::NodeId,
    goal: G::NodeId,
    edge_cost: F,
    node_color_fn: ColorFn,
    limit: Option<usize>,
) -> (HashMap<(G::NodeId, C), (G::NodeId, C)>, G::NodeId, C)
where
    G: IntoEdges + Visitable,
    G::NodeId: Eq + Hash,
    F: FnMut(G::EdgeRef) -> K,
    ColorFn: FnMut(&G::NodeId) -> C,
    K: Measure + Copy,
    C: num_traits::PrimInt + num_traits::Unsigned + Hash + std::fmt::Display,
{
    dijkstra_with_node_colors_and_bounds(
        graph,
        start,
        goal,
        edge_cost,
        node_color_fn,
        |_| Some((K::default(), 0)),
        limit,
    )
}

pub fn dijkstra_with_node_colors_and_bounds<G, F, K, C, ColorFn, BoundFn>(
    graph: G,
    start: G::NodeId,
    goal: G::NodeId,
    mut edge_cost: F,
    mut node_color_fn: ColorFn,
    mut lower_bound_fn: BoundFn,
    limit: Option<usize>,
) -> (HashMap<(G::NodeId, C), (G::NodeId, C)>, G::NodeId, C)
where
    G: IntoEdges + Visitable,
    G::NodeId: Eq + Hash,
    F: FnMut(G::EdgeRef) -> K,
    ColorFn: FnMut(&G::NodeId) -> C,
    BoundFn: FnMut(&G::NodeId) -> Option<(K, usize)>,
    K: Measure + Copy,
    C: num_traits::PrimInt + num_traits::Unsigned + Hash + std::fmt::Display,
{
    let limit = limit.unwrap_or(0);
    let mut scores: HashMap<(G::NodeId, C), K> = HashMap::new();
    let mut predecessor: HashMap<(G::NodeId, C), (G::NodeId, C)> = HashMap::new();
    let mut labels: HashMap<G::NodeId, Vec<Label<K, C>>> = HashMap::new();
    let mut visit_next: BinaryHeap<MinScored<K, (K, C, usize, G::NodeId)>> = BinaryHeap::new();
    let zero_score = K::default();
    let start_color = node_color_fn(&start);
    let start_bound = match lower_bound_fn(&start) {
        Some(bound) => bound,
        None => return (predecessor, start.clone(), start_color),
    };
    scores.insert((start, start_color), zero_score);
    labels.entry(start).or_default().push(Label {
        mask: start_color,
        cost: zero_score,
        depth: 0,
    });
    visit_next.push(MinScored(
        zero_score + start_bound.0,
        (zero_score, start_color, 0 as usize, start),
    ));
    while let Some(MinScored(_priority, (node_score, current_color, node_depth, node))) =
        visit_next.pop()
    {
        if scores.get(&(node, current_color)) != Some(&node_score) {
            continue;
        }
        let current_label = Label {
            mask: current_color,
            cost: node_score,
            depth: node_depth,
        };
        if !labels
            .get(&node)
            .map(|node_labels| node_labels.contains(&current_label))
            .unwrap_or(false)
        {
            continue;
        }
        if goal == node {
            return (predecessor, node.clone(), current_color);
        }
        match lower_bound_fn(&node) {
            Some((_min_cost, min_hops)) => {
                if limit > 0 && node_depth + min_hops > limit {
                    continue;
                }
            }
            None => continue,
        }
        if limit > 0 && node_depth + 1 > limit {
            continue;
        }
        for edge in graph.edges(node) {
            let next = edge.target();
            let next_color = node_color_fn(&next);
            if (next_color & current_color) != C::min_value() {
                continue;
            }
            let next_depth = node_depth + 1;
            let next_bound = match lower_bound_fn(&next) {
                Some(bound) => bound,
                None => continue,
            };
            if limit > 0 && next_depth + next_bound.1 > limit {
                continue;
            }
            let next_score = node_score + edge_cost(edge);
            let next_color = current_color | next_color;
            let next_priority = next_score + next_bound.0;
            let next_label = Label {
                mask: next_color,
                cost: next_score,
                depth: next_depth,
            };
            match scores.entry((next, next_color)) {
                Occupied(ent) => {
                    if next_score < *ent.get() {
                        if !insert_nondominated_label(labels.entry(next).or_default(), next_label) {
                            continue;
                        }
                        *ent.into_mut() = next_score;
                        visit_next.push(MinScored(
                            next_priority,
                            (next_score, next_color, next_depth, next),
                        ));
                        predecessor
                            .insert((next.clone(), next_color), (node.clone(), current_color));
                    }
                }
                Vacant(ent) => {
                    if !insert_nondominated_label(labels.entry(next).or_default(), next_label) {
                        continue;
                    }
                    ent.insert(next_score);
                    visit_next.push(MinScored(
                        next_priority,
                        (next_score, next_color, next_depth, next),
                    ));
                    predecessor.insert((next.clone(), next_color), (node.clone(), current_color));
                }
            }
        }
    }
    (predecessor, start.clone(), start_color)
}

#[cfg(test)]
mod tests {
    use super::{
        dijkstra_with_node_colors, dijkstra_with_node_colors_and_bounds, insert_nondominated_label,
        Label,
    };
    use petgraph::graphmap::DiGraphMap;
    use std::collections::HashMap;

    #[test]
    fn insert_nondominated_label_rejects_dominated_candidate() {
        let mut labels = vec![Label {
            mask: 0b001u8,
            cost: 5usize,
            depth: 2,
        }];

        assert!(!insert_nondominated_label(
            &mut labels,
            Label {
                mask: 0b011,
                cost: 5,
                depth: 3,
            }
        ));
        assert_eq!(
            labels,
            vec![Label {
                mask: 0b001,
                cost: 5,
                depth: 2,
            }]
        );
    }

    #[test]
    fn insert_nondominated_label_removes_dominated_existing_label() {
        let mut labels = vec![Label {
            mask: 0b011u8,
            cost: 5usize,
            depth: 3,
        }];

        assert!(insert_nondominated_label(
            &mut labels,
            Label {
                mask: 0b001,
                cost: 4,
                depth: 2,
            }
        ));
        assert_eq!(
            labels,
            vec![Label {
                mask: 0b001,
                cost: 4,
                depth: 2,
            }]
        );
    }

    #[test]
    fn insert_nondominated_label_keeps_non_comparable_labels() {
        let mut labels = vec![Label {
            mask: 0b001u8,
            cost: 6usize,
            depth: 1,
        }];

        assert!(insert_nondominated_label(
            &mut labels,
            Label {
                mask: 0b010,
                cost: 4,
                depth: 3,
            }
        ));
        assert_eq!(
            labels,
            vec![
                Label {
                    mask: 0b001,
                    cost: 6,
                    depth: 1,
                },
                Label {
                    mask: 0b010,
                    cost: 4,
                    depth: 3,
                },
            ]
        );
    }

    fn reconstruct_path(
        predecessor: &HashMap<(usize, usize), (usize, usize)>,
        goal: usize,
        node: usize,
        mask: usize,
    ) -> Vec<usize> {
        let mut path = Vec::new();
        if node == goal {
            let mut current = node;
            let mut current_mask = mask;
            path.push(current);
            while let Some((previous, previous_mask)) = predecessor.get(&(current, current_mask)) {
                path.push(*previous);
                current = *previous;
                current_mask = *previous_mask;
            }
        }
        path.reverse();
        path
    }

    #[test]
    fn dijkstra_with_bounds_matches_unbounded_result() {
        let edges = [
            (0usize, 1usize, 10usize),
            (1, 4, 10),
            (0, 2, 1),
            (2, 3, 1),
            (3, 4, 1),
        ];
        let graph = DiGraphMap::<usize, usize>::from_edges(&edges);
        let node_colors = [1usize, 2, 4, 8, 16];

        let (unbounded_predecessor, unbounded_node, unbounded_mask) = dijkstra_with_node_colors(
            &graph,
            0,
            4,
            |edge| *edge.2,
            |node| node_colors[*node],
            Some(4),
        );
        let (bounded_predecessor, bounded_node, bounded_mask) =
            dijkstra_with_node_colors_and_bounds(
                &graph,
                0,
                4,
                |edge| *edge.2,
                |node| node_colors[*node],
                |node| match *node {
                    0 => Some((3, 3)),
                    1 => Some((10, 1)),
                    2 => Some((2, 2)),
                    3 => Some((1, 1)),
                    4 => Some((0, 0)),
                    _ => None,
                },
                Some(4),
            );

        assert_eq!(
            reconstruct_path(&bounded_predecessor, 4, bounded_node, bounded_mask),
            reconstruct_path(&unbounded_predecessor, 4, unbounded_node, unbounded_mask)
        );
    }

    #[test]
    fn dijkstra_with_bounds_skips_nodes_disconnected_from_goal() {
        let edges = [(0usize, 1usize, 1usize), (1, 2, 1), (0, 3, 1)];
        let graph = DiGraphMap::<usize, usize>::from_edges(&edges);
        let node_colors = [1usize, 2, 4, 8];

        let (predecessor, node, mask) = dijkstra_with_node_colors_and_bounds(
            &graph,
            0,
            3,
            |edge| *edge.2,
            |node| node_colors[*node],
            |node| match *node {
                0 => Some((1, 1)),
                3 => Some((0, 0)),
                _ => None,
            },
            Some(3),
        );

        assert_eq!(reconstruct_path(&predecessor, 3, node, mask), vec![0, 3]);
    }
}
