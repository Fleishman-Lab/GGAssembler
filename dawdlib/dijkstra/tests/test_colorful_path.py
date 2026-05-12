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

