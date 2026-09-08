"""Discrete transition counting with smoothing."""
import numpy as np
def transition_matrix(current,next_state,n_states,alpha=1.0):
    counts=np.full((n_states,n_states),alpha,float)
    for a,b in zip(current,next_state): counts[int(a),int(b)]+=1
    return counts/counts.sum(1,keepdims=True)
