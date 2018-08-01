import os
import json
import itertools
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

ROI_PATH = 'behavior_rois.json'
MODEL_DIR = 'model'
OUTPUT_DIR = 'output'

ALPHA = 0.01
LEVELS = ('leaves', 'ss')
METRICS = ('connection_density', 'connection_strength',
           'normalized_connection_density', 'normalized_connection_strength')


def main():
    # read in ROIs
    with open(ROI_PATH, 'r') as f:
        task_rois = json.load(f)

    for task, rois in task_rois.items():
        for level in LEVELS:
            result = defaultdict(dict)
            for metric in METRICS:

                # read in table
                path = os.path.join(MODEL_DIR, '%s_%s.csv' % (level, metric))
                model = pd.read_csv(path, index_col=0)

                # need small constant for log
                eps = np.percentile(model.values[model.values.nonzero()], 2)

                # get rois found in this table
                common_rois = set(model.index).intersection(rois)
                if len(common_rois) < 2:
                    continue

                # get network connectivity
                network = model.loc[common_rois, common_rois].values.ravel()
                network = np.log10(network + eps)

                popmean = np.log10(model.values[model.values > eps]).mean()

                # perform t test
                prob = stats.ttest_1samp(network, popmean)[1]

                result[metric]['p_value'] = prob
                result[metric]['significant'] = str(prob < ALPHA).lower()
                result[metric]['alpha'] = ALPHA
                result[metric]['network_rois'] = len(common_rois)
                result[metric]['total_rois'] = len(rois)

            if not result:
                continue

            # meta data
            result['task'] = task
            result['level'] = level

            path = os.path.join(OUTPUT_DIR, '%s_%s_t-test.json' % (task, level))
            with open(path, 'w') as f:
                json.dump(result, f, indent=2)


if __name__ == '__main__':
    main()
