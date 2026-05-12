import networkx as nx

from dawdlib.dijkstra import colorful, colourful_dijkstra


def test_random_recolor_ors_colliding_random_colors(monkeypatch):
    finder = object.__new__(colorful.ShortestPathFinder)
    finder.graph = nx.Graph()
    finder.graph.add_node("n0")
    finder.node_map = {"n0": 0}
    finder.gate_colors = {"n0": ["a", "b"]}
    finder.all_gate_colors = {"a", "b"}
    finder.num_of_colors = [2]

    monkeypatch.setattr(colorful.np.random, "choice", lambda selected: selected[1])

    assert finder.random_recolor(len_cutoff=None, no_colors=3) == {0: 2}


def test_colourful_shortest_path_revisits_node_with_different_color_mask():
    edges = [
        (0, 1, 1),
        (1, 3, 1),
        (0, 2, 3),
        (2, 3, 1),
        (3, 4, 1),
    ]
    node_colors = {
        0: 1,
        1: 2,
        2: 4,
        3: 0,
        4: 2,
    }

    path = colourful_dijkstra.colourful_shortest_path(
        edges, 0, 4, node_colors, 4
    )

    assert path == [0, 2, 3, 4]


def test_colourful_shortest_path_keeps_existing_public_api(monkeypatch):
    graph = nx.DiGraph()
    graph.add_edge("source", "a", weight=1)
    graph.add_edge("a", "target", weight=1)

    finder = object.__new__(colorful.ShortestPathFinder)
    finder.graph = graph
    finder.source = "source"
    finder.target = "target"
    finder.node_map = {node: index for index, node in enumerate(graph.nodes)}
    finder.idx_node_map = {index: node for node, index in finder.node_map.items()}
    finder.graph_edges = [
        (finder.node_map[u], finder.node_map[v], data["weight"])
        for u, v, data in graph.edges(data=True)
    ]

    monkeypatch.setattr(
        finder,
        "random_recolor",
        lambda len_cutoff, no_colors: {
            finder.node_map["source"]: 1,
            finder.node_map["a"]: 2,
            finder.node_map["target"]: 4,
        },
    )

    assert finder.find_shortest_path(len_cutoff=2, no_colors=3) == [
        "source",
        "a",
        "target",
    ]


def test_colourful_shortest_path_respects_len_cutoff():
    edges = [
        (0, 1, 1),
        (1, 2, 1),
        (2, 3, 1),
    ]
    node_colors = {
        0: 1,
        1: 2,
        2: 4,
        3: 8,
    }

    assert colourful_dijkstra.colourful_shortest_path(
        edges, 0, 3, node_colors, 2
    ) == []
    assert colourful_dijkstra.colourful_shortest_path(
        edges, 0, 3, node_colors, 3
    ) == [0, 1, 2, 3]


def test_colourful_shortest_path_prefers_lower_weight_over_fewer_edges():
    edges = [
        (0, 1, 10),
        (1, 4, 10),
        (0, 2, 1),
        (2, 3, 1),
        (3, 4, 1),
    ]
    node_colors = {
        0: 1,
        1: 2,
        2: 4,
        3: 8,
        4: 16,
    }

    assert colourful_dijkstra.colourful_shortest_path(
        edges, 0, 4, node_colors, 4
    ) == [0, 2, 3, 4]
