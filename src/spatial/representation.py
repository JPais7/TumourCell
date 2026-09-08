"""Coordinates, neighbourhood graphs, composition, niches and uncertainty."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree
@dataclass(frozen=True)
class SpatialContext:
    neighbours: tuple[int,...]; composition: dict[str,float]; niche: str; uncertainty: float
def radius_graph(coordinates, radius):
    x=np.asarray(coordinates,float)
    if x.ndim != 2 or x.shape[1] not in (2,3): raise ValueError("coordinates must be n x 2 or n x 3")
    return [tuple(j for j in js if j != i) for i,js in enumerate(cKDTree(x).query_ball_point(x,radius))]
def neighbourhood_context(graph, labels, index, minimum_neighbours=5):
    neighbours=graph[index]; counts={k:sum(labels[j]==k for j in neighbours) for k in set(labels)}
    n=len(neighbours); comp={k:v/n for k,v in counts.items()} if n else {}
    niche=max(comp,key=comp.get) if comp else "UNRESOLVED"
    return SpatialContext(neighbours,comp,niche,float("inf") if n==0 else max(1,n/minimum_neighbours)**-0.5)
