import os

import numpy as np
from mcmodels.core import VoxelModelCache, Mask
from mcmodels.models.voxel import RegionalizedModel


MANIFEST_FILE = os.path.join(os.path.expanduser('~'), 'connectivity',
                             'voxel_model_manifest.json')
ROI_PATH = 'behavior_rois.json'
OUTPUT_DIR = 'model'

TARGET_SEARCH_SET_ID = 184527634
SUMMARY_STRUCTURES_SET_ID = 687527945
METRICS = ('connection_density', 'connection_strength',
           'normalized_connection_density', 'normalized_connection_strength')

def get_summary_structure_map(tree, structure_ids):
    all_targets = tree.get_structures_by_set_id([SUMMARY_STRUCTURES_SET_ID])
    return {s['id'] : s['acronym'] for s in all_targets if s['id'] in structure_ids}

def get_leaves_map(tree, structure_ids):
    all_targets = tree.get_structures_by_set_id([TARGET_SEARCH_SET_ID])
    n_children = list(map(len, tree.child_ids([s['id'] for s in all_targets])))
    return {s['id'] : s['acronym'] for s, n in zip(all_targets, n_children)
            if s['id'] in structure_ids and n == 0}

def main():

    # initialize cache object
    cache = VoxelModelCache(manifest_file=MANIFEST_FILE)

    # get mapping from acronym to id
    tree = cache.get_structure_tree()
    structure_ids = np.unique(cache.get_reference_space().annotation)

    # load in voxel model
    voxel_array, source_mask, target_mask = cache.get_voxel_connectivity_array()

    for structure_set in ('ss', 'leaves'):
        if structure_set == 'ss':
            id_acronym_map = get_summary_structure_map(tree, structure_ids)
        else:
            id_acronym_map = get_leaves_map(tree, structure_ids)

        # regionalize (NOTE: bilateral sum)
        source_key = source_mask.get_key(structure_ids=id_acronym_map.keys())
        target_key = target_mask.get_key(structure_ids=id_acronym_map.keys())
        regional_model = RegionalizedModel.from_voxel_array(
            voxel_array, source_key, target_key, dataframe=True)

        tables = [getattr(regional_model, metric) for metric in METRICS]
        for table, metric in zip(tables, METRICS):
            # get model metric and rename index/columns
            table.index = [id_acronym_map[v] for v in table.index.values]
            table.columns = [id_acronym_map[v] for v in table.columns.values]

            table.to_csv(os.path.join(OUTPUT_DIR, '%s_%s.csv' % (structure_set, metric)))


if __name__ == '__main__':
    main()
