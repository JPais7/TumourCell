"""Conservative assignment of cells to externally supported clone profiles."""
def assign_clone(distances, maximum_distance):
    if not distances: return "UNASSIGNED"
    clone, distance=min(distances.items(), key=lambda x:x[1])
    return clone if distance <= maximum_distance else "UNASSIGNED"
