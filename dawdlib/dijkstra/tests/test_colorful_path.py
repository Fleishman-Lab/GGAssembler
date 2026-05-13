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


def test_random_recolor_dense_matches_dict_with_fixed_seed():
    graph = nx.Graph()
    graph.add_nodes_from(["n0", "n1", "n2"])

    finder = object.__new__(colorful.ShortestPathFinder)
    finder.graph = graph
    finder.node_map = {"n0": 0, "n1": 1, "n2": 2}
    finder.gate_colors = {
        "n0": ["a", "b"],
        "n1": ["b", "c"],
        "n2": [],
    }
    finder.all_gate_colors = {"a", "b", "c"}
    finder._all_gate_colors_order = ["a", "b", "c"]
    finder._gate_color_id_map = {"a": 0, "b": 1, "c": 2}
    finder._gate_color_ids_by_node = [(0, 1), (1, 2), ()]
    finder.num_of_colors = [3]

    colorful.np.random.seed(123)
    dict_colors = finder.random_recolor(len_cutoff=None, no_colors=4)
    colorful.np.random.seed(123)
    dense_colors = finder._random_recolor_dense(len_cutoff=None, no_colors=4)

    assert dense_colors == [
        dict_colors[0],
        dict_colors[1],
        dict_colors[2],
    ]


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


def test_rust_resident_finder_matches_module_function():
    edges = [
        (0, 1, 1),
        (1, 4, 10),
        (0, 2, 2),
        (2, 3, 2),
        (3, 4, 2),
    ]
    node_colors = {
        0: 1,
        1: 2,
        2: 4,
        3: 8,
        4: 16,
    }

    expected = colourful_dijkstra.colourful_shortest_path(
        edges, 0, 4, node_colors, 4
    )
    finder = colourful_dijkstra.ColourfulPathFinder(edges, 0, 4)

    assert finder.find_shortest_path(node_colors, 4) == expected


def test_rust_resident_digraph_finder_matches_module_function():
    edges = [
        (0, 1, 1),
        (1, 4, 10),
        (0, 2, 2),
        (2, 3, 2),
        (3, 4, 2),
    ]
    node_colors = {
        0: 1,
        1: 2,
        2: 4,
        3: 8,
        4: 16,
    }

    expected = colourful_dijkstra.colourful_shortest_path(
        edges, 0, 4, node_colors, 4
    )
    finder = colourful_dijkstra.ColourfulPathFinderDiGraph(edges, 0, 4, 5)

    assert finder.find_shortest_path(node_colors, 4) == expected


def test_rust_resident_finders_accept_dense_node_colors():
    edges = [
        (0, 1, 1),
        (1, 4, 10),
        (0, 2, 2),
        (2, 3, 2),
        (3, 4, 2),
    ]
    node_colors = [1, 2, 4, 8, 16]
    node_color_map = dict(enumerate(node_colors))
    expected = colourful_dijkstra.colourful_shortest_path(
        edges, 0, 4, node_color_map, 4
    )

    graphmap_finder = colourful_dijkstra.ColourfulPathFinder(edges, 0, 4)
    digraph_finder = colourful_dijkstra.ColourfulPathFinderDiGraph(edges, 0, 4, 5)

    assert (
        graphmap_finder.find_shortest_path_with_node_colors(node_colors, 4)
        == expected
    )
    assert (
        digraph_finder.find_shortest_path_with_node_colors(node_colors, 4)
        == expected
    )


def test_rust_resident_search_modes_match_on_simple_path():
    edges = [
        (0, 1, 10),
        (1, 4, 10),
        (0, 2, 1),
        (2, 3, 1),
        (3, 4, 1),
    ]
    node_colors = [1, 2, 4, 8, 16]
    expected = [0, 2, 3, 4]
    finder = colourful_dijkstra.ColourfulPathFinderDiGraph(edges, 0, 4, 5)

    for use_a_star in (False, True):
        for use_dominance in (False, True):
            assert finder.find_shortest_path_with_node_colors(
                node_colors,
                4,
                use_a_star,
                use_dominance,
            ) == expected


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
    finder._rust_finder = colourful_dijkstra.ColourfulPathFinder(
        finder.graph_edges,
        finder.node_map[finder.source],
        finder.node_map[finder.target],
    )

    monkeypatch.setattr(
        finder,
        "_random_recolor_dense",
        lambda len_cutoff, no_colors: [1, 2, 4],
    )

    assert finder.find_shortest_path(len_cutoff=2, no_colors=3) == [
        "source",
        "a",
        "target",
    ]


def test_shortest_path_finder_uses_resident_rust_finder(monkeypatch):
    graph = nx.DiGraph()
    graph.add_edge("source", "a", weight=1)
    graph.add_edge("a", "target", weight=1)

    gate_colors = {
        "source": [],
        "a": ["a"],
        "target": ["target"],
    }
    monkeypatch.setattr(
        colorful,
        "color_gates",
        lambda graph, ggdata: (gate_colors, [len(gate_colors)]),
    )

    finder = colorful.ShortestPathFinder(
        graph,
        ggdata=object(),
        source="source",
        target="target",
    )

    assert hasattr(finder, "_rust_finder")

    monkeypatch.setattr(
        finder,
        "_random_recolor_dense",
        lambda len_cutoff, no_colors: [1, 2, 4],
    )

    assert finder.find_shortest_path(len_cutoff=2, no_colors=3) == [
        "source",
        "a",
        "target",
    ]


def test_shortest_path_finder_find_many_returns_original_nodes(monkeypatch):
    graph = nx.DiGraph()
    graph.add_edge("source", "a", weight=1)
    graph.add_edge("a", "target", weight=1)

    gate_colors = {
        "source": [],
        "a": [],
        "target": [],
    }
    monkeypatch.setattr(
        colorful,
        "color_gates",
        lambda graph, ggdata: (gate_colors, [0]),
    )

    finder = colorful.ShortestPathFinder(
        graph,
        ggdata=object(),
        source="source",
        target="target",
    )

    results = finder.find_many(
        min_gates=2,
        max_gates=2,
        retries=3,
        seed=123,
    )

    assert results == [
        (2, 0, ["source", "a", "target"]),
        (2, 1, ["source", "a", "target"]),
        (2, 2, ["source", "a", "target"]),
    ]
    assert results == finder.find_many(
        min_gates=2,
        max_gates=2,
        retries=3,
        seed=123,
    )


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
