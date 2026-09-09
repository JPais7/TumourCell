"""Coordinates, neighbourhood graphs, composition, niches and uncertainty."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree
from dataclasses import field
@dataclass(frozen=True)
class SpatialObservation:
    cell_id: str; sample_id: str; x: float; y: float; latent_state: tuple[float,...]
    clone_id: str | None = None; neighbourhood: tuple[int,...] = (); modality: str = "unknown"; qc: str = "UNASSESSED"
@dataclass(frozen=True)
class SpatialContext:
    neighbours: tuple[int,...]; composition: dict[str,float]; niche: str; uncertainty: float
def radius_graph(coordinates, radius):
    if coordinates is None: raise ValueError("cell coordinates are required for spatial analysis")
    x=np.asarray(coordinates,float)
    if x.ndim != 2 or x.shape[1] not in (2,3): raise ValueError("coordinates must be n x 2 or n x 3")
    return [tuple(j for j in js if j != i) for i,js in enumerate(cKDTree(x).query_ball_point(x,radius))]
def validate_coordinates(coordinates):
    x=np.asarray(coordinates,float)
    if not np.isfinite(x).all(): raise ValueError("coordinates contain missing or non-finite values")
    if x.ndim != 2 or x.shape[1] not in (2,3): raise ValueError("coordinates must be n x 2 or n x 3")
    return x

def build_neighbourhood_graph(coordinates, method="knn", k=6, radius=None):
    x=validate_coordinates(coordinates)
    if method=="radius":
        if radius is None or radius<=0: raise ValueError("positive radius required")
        return radius_graph(x,radius)
    if method!="knn" or k<1 or k>=len(x): raise ValueError("method must be knn with 1 <= k < n")
    d=cKDTree(x).query(x,k=k+1)[1][:,1:]
    return [tuple(row.tolist()) for row in d]

def neighbourhood_features(latent_state, graph, labels=None, distances=None):
    z=np.asarray(latent_state,float); out=[]
    for i,neigh in enumerate(graph):
        if not neigh: out.append({"local_mean_state":None,"local_variance_state":None,"composition":{},"n_neighbours":0,"isolated":True}); continue
        local=z[list(neigh)]; weights=np.ones(len(neigh))
        if distances is not None: weights=1/np.maximum(np.asarray(distances[i]),1e-12)
        mean=np.average(local,axis=0,weights=weights)
        comp={}
        if labels is not None:
            vals,counts=np.unique(np.asarray(labels)[list(neigh)],return_counts=True); comp={str(v):float(c/len(neigh)) for v,c in zip(vals,counts)}
        out.append({"local_mean_state":mean.tolist(),"local_variance_state":local.var(0).tolist(),"composition":comp,"n_neighbours":len(neigh),"isolated":False})
    return out

def spatial_confounding(global_before, global_after, regions_before, regions_after):
    """Descriptive decomposition; never assigns causality."""
    a=np.asarray(global_after,float); b=np.asarray(global_before,float)
    rb={k:np.mean(np.asarray(v),axis=0) for k,v in regions_before.items()}; ra={k:np.mean(np.asarray(v),axis=0) for k,v in regions_after.items()}
    common=set(rb)&set(ra); within=np.mean([ra[k]-rb[k] for k in common],axis=0) if common else None
    result={"global_state_shift":(a-b).tolist(),"within_region_shift":None if within is None else within.tolist(),"between_region_composition_shift":None,"identifiability_status":"NOT_IDENTIFIABLE" if not common else "INDETERMINATE"}
    return result
def neighbourhood_context(graph, labels, index, minimum_neighbours=5):
    neighbours=graph[index]; counts={k:sum(labels[j]==k for j in neighbours) for k in set(labels)}
    n=len(neighbours); comp={k:v/n for k,v in counts.items()} if n else {}
    niche=max(comp,key=comp.get) if comp else "UNRESOLVED"
    return SpatialContext(neighbours,comp,niche,float("inf") if n==0 else max(1,n/minimum_neighbours)**-0.5)
