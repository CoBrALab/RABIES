import sys
import numpy as np

""" countCommonEdge
Function that calculates the number pf pairs of points clustered together (or,
the number of common edges) (Ben-Hur, 2002) in a memory-efficient manner.

It takes as input two identically-long vectors (matrix with a single column).
    Each index across the two vectors are assumed to represent the same point /
    subject / thing.

Input
    - vec1: numpy vector where each element is the group number of that particular point
    - vec2: numpy vector where each element is the group number of that particular point
    - pointQuant: the total number of points (should be identical to len(vec1))
Output
    an integer representing the number of common edges / pairs of points clustered
        together across the two cluster outcomes
"""
def countCommonEdge (vec1, vec2, pointQuant, verbose=True):
    #Initialize variable to count the stable edges across the two sets
    stableEdges = 0

    #Iterate over each element of each vectors
    for idx in range(0,pointQuant):
        #User indication
        if verbose:
            if idx%100 == 0 or idx == pointQuant-1:
                sys.stdout.write("Counting edges: %d out of %d vertices\r"%(idx+1,pointQuant))
                sys.stdout.flush()
        #Make a intra-vector correlation vector
        corVec1 = (vec1 == vec1[idx])
        corVec2 = (vec2 == vec2[idx])
        #Set the element's correlation itself to False (since it is always True)
        corVec1[idx] = False
        corVec2[idx] = False
        #Calculate the common edges for this element
        overlapVec = np.logical_and(corVec1, corVec2)
        #Aggregate the the common edge count
        stableEdges += np.sum(overlapVec)
        #Clear up the vectors for memory conservation
        del corVec1
        del corVec2
        del overlapVec

    #Return the common edge count
    return stableEdges
