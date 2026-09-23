use std::collections::{BinaryHeap, HashMap};

use crate::scored::MinScored;

type State = (usize, usize, usize);
type Predecessors = HashMap<State, State>;
type Frontiers = HashMap<(usize, usize), Vec<(usize, usize)>>;
type Queue = BinaryHeap<MinScored<usize, State>>;

/// Find the cheapest colorful path through an integer-indexed adjacency list.
pub fn dijkstra(
    graph: &[Vec<(usize, usize)>],
    start: usize,
    goal: usize,
    node_colors: &[usize],
    limit: Option<usize>,
) -> (Predecessors, State) {
    let limit = limit.unwrap_or(0);
    let mut scores: HashMap<State, usize> = HashMap::new();
    let mut frontiers: Frontiers = HashMap::new();
    let mut predecessor: Predecessors = HashMap::new();
    let mut visit_next: Queue = BinaryHeap::new();
    let zero_score = 0;
    let start_color = node_colors[start];
    scores.insert((start, start_color, 0), zero_score);
    frontiers.insert((start, start_color), vec![(0, zero_score)]);
    visit_next.push(MinScored(zero_score, (start_color, 0usize, start)));
    while let Some(MinScored(node_score, (node_color, node_depth, node))) = visit_next.pop() {
        if scores.get(&(node, node_color, node_depth)) != Some(&node_score) {
            continue;
        }
        if goal == node {
            return (predecessor, (node, node_color, node_depth));
        }
        if limit > 0 && node_depth >= limit {
            continue;
        }
        for &(next, edge_cost) in &graph[node] {
            let next_color = node_colors[next];
            if (next_color & node_color) != 0 {
                continue;
            }
            let next_score = node_score + edge_cost;
            let next_color = node_color | next_color;
            let next_depth = node_depth + 1;
            let frontier_key = (next, next_color);
            if let Some(labels) = frontiers.get(&frontier_key) {
                if labels
                    .iter()
                    .any(|(depth, score)| *depth <= next_depth && *score <= next_score)
                {
                    continue;
                }
            }

            let dominated_depths: Vec<usize> = frontiers
                .get(&frontier_key)
                .into_iter()
                .flatten()
                .filter(|(depth, score)| next_depth <= *depth && next_score <= *score)
                .map(|(depth, _)| *depth)
                .collect();
            if let Some(labels) = frontiers.get_mut(&frontier_key) {
                labels.retain(|(depth, score)| !(next_depth <= *depth && next_score <= *score));
            }
            for depth in dominated_depths {
                scores.remove(&(next, next_color, depth));
            }
            frontiers
                .entry(frontier_key)
                .or_default()
                .push((next_depth, next_score));

            scores.insert((next, next_color, next_depth), next_score);
            visit_next.push(MinScored(next_score, (next_color, next_depth, next)));
            predecessor.insert(
                (next, next_color, next_depth),
                (node, node_color, node_depth),
            );
        }
    }
    (predecessor, (start, start_color, 0))
}
