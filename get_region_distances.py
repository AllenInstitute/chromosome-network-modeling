import os

import numpy as np
import pandas as pd
from scipy.spatial import distance
from mcmodels.core import VoxelModelCache

MANIFEST_FILE = os.path.join(os.path.expanduser('~'), 'connectivity',
                             'voxel_model_manifest.json')
MAJOR_BRAIN_SET_ID = 687527670
OUTPUT_DIR = 'model'


def main():
    # initialize cache object
    cache = VoxelModelCache(manifest_file=MANIFEST_FILE)
    rs = cache.get_reference_space()
    rs.remove_unassigned(update_self=True)

    # major id children only
    structures = []
    for s in rs.structure_tree.get_structures_by_set_id([MAJOR_BRAIN_SET_ID]):
        structures.extend(rs.structure_tree.descendants([s['id']])[0])

    print('getting masks')
    centroids = []
    for s in structures:
        mask = cache.get_structure_mask(s['id'])[0]

        # NOTE: ipsi
        mask[..., :mask.shape[-1]//2] = 0

        centroids.append(np.argwhere(mask).mean(axis=0))

    print('computing distances')
    distances = distance.squareform(distance.pdist(np.asarray(centroids)))

    # convert from 100um to mm
    distances *= 0.1

    structure_acronyms = [s['acronym'] for s in structures]
    distances = pd.DataFrame(
        distances, index=structure_acronyms, columns=structure_acronyms)

    # save
    distances.to_csv(os.path.join(OUTPUT_DIR, 'distances.csv'))


if __name__ == '__main__':
    main()
