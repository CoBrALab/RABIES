def compute_edge_mask(mask_array, num_edge_voxels=1):
    import numpy as np
    #custom function for computing edge mask from an input brain mask
    shape = mask_array.shape

    #iterate through all voxels from the three dimensions and look if it contains surrounding voxels
    edge_mask = np.zeros(shape, dtype=bool)
    num_voxel = 0
    while num_voxel < num_edge_voxels:
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    #only look if the voxel is part of the mask
                    if mask_array[x, y, z]:
                        if (mask_array[x-1:x+2, y-1:y+2, z-1:z+2] == 0).sum() > 0:
                            edge_mask[x, y, z] = 1
        mask_array = mask_array-edge_mask
        num_voxel += 1

    return edge_mask