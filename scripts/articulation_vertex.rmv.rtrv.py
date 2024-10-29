import networkx as nx
from itertools import combinations
import random
from util import hash_as_str

random.seed(0)

def edge_distance(vec1, vec2):
    assert len(vec1) == len(vec2)
    assert hash_as_str(vec1) != hash_as_str(vec2)
    penalty = 0
    for i in range(len(vec1)):
        if vec1[i] != vec2[i]:
            if i >= 1 and  vec1[i] == vec1[i - 1] and vec2[i] == vec2[i - 1]:
                penalty += 0
            else:
                penalty += 1
    return penalty


def analyze_isoform_graph(isoforms):
    """
    Create and analyze an isoform graph, identifying connected components
    and articulation points.
    """
    # Create the graph
    G = nx.Graph()
    for i, isoform in enumerate(isoforms):
        G.add_node(i, binary=isoform)
    
    # Add edges between isoforms with edge_distance = 1
    for (i1, iso1), (i2, iso2) in combinations(enumerate(isoforms), 2):
        if edge_distance(iso1, iso2) == 1:
            G.add_edge(i1, i2)
    
    # Find connected components
    components = list(nx.connected_components(G))
    
    # Find articulation points
    articulation_points = list(nx.articulation_points(G))
        
    return G, components, articulation_points

# Example usage with sample isoforms
isoforms = [
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, 1, 1],
    [0, 1, 1, 0],
]

# Analyze the graph
G, components, articulation_points = analyze_isoform_graph(isoforms)

# Print analysis results
print(f"\nGraph Analysis Results:")
print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")
print(f"\nConnected Components ({len(components)}):")
for i, component in enumerate(components, 1):
    print(f"Component {i}: {sorted(component)} -> {[G.nodes[n]['binary'] for n in sorted(component)]}")

print(f"\nArticulation Points: {sorted(articulation_points)}")
if articulation_points:
    print("Articulation points with their binary vectors:")
    for point in sorted(articulation_points):
        print(f"Node {point}: {G.nodes[point]['binary']}")
        # Show which components would be created if this point were removed
        H = G.copy()
        H.remove_node(point)
        subcomponents = list(nx.connected_components(H))
        print(f"  If removed, creates {len(subcomponents)} subcomponents:")
        for i, subcomp in enumerate(subcomponents, 1):
            print(f"  Subcomponent {i}: {sorted(subcomp)}")
else:
    print("No articulation points found in the graph.")