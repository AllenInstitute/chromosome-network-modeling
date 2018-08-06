import os
import json
from mcmodels.core import VoxelModelCache

MANIFEST_FILE = os.path.join(os.path.expanduser('~'), 'connectivity',
                             'voxel_model_manifest.json')

ROI_PATH = 'behavior_rois.json'
LEAVES_PATH = 'deepest_rois.json'


def main():
    # read in ROIs
    with open(ROI_PATH, 'r') as f:
        task_rois = json.load(f)

    # initialize cache object
    cache = VoxelModelCache(manifest_file=MANIFEST_FILE)

    # get mapping from acronym to id
    tree = cache.get_structure_tree()
    acronym_id_map = tree.get_id_acronym_map()
    id_acronym_map = {v : k for k, v in acronym_id_map.items()}

    result = dict()
    for task, rois in task_rois.items():
        roi_ids = [acronym_id_map[k] for k in rois]
        overlapping_ancestors = tree.has_overlaps(roi_ids)

        leaves = set(roi_ids).difference(overlapping_ancestors)
        result[task] = [id_acronym_map[k] for k in leaves]

    with open(LEAVES_PATH, 'w') as f:
        json.dump(result, f, indent=2)


if __name__ == '__main__':
    main()
