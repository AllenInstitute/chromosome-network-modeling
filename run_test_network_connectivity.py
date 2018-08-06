import os
import json
import itertools
from collections import defaultdict

from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


ROI_PATH = 'deepest_rois.json'
MODEL_DIR = 'model'
OUTPUT_DIR = 'output'

ALPHA = 0.01
LEVELS = ('leaves', 'ss', 'deepest')
METRICS = ('connection_density', 'connection_strength',
           'normalized_connection_density', 'normalized_connection_strength')

def plot_dists(population, network, title):
    fig, ax1 = plt.subplots()
    sns.distplot(population, ax=ax1, color=sns.xkcd_rgb['dusky blue'])
    ax1.set_xlabel('logarithmic connectivity weight')
    ax1.set_ylabel('population count')

    ax2 = ax1.twinx()
    sns.distplot(network, ax=ax2, color=sns.xkcd_rgb['faded orange'])
    ax2.set_ylabel('ROI network count')

    plt.title(title)
    return fig


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
                eps = np.percentile(model.values[model.values.nonzero()], 3)

                # get rois found in this table
                common_rois = set(model.index).intersection(rois)

                # get network connectivity
                network = model.loc[common_rois, common_rois].values.ravel()
                network = np.log10(network[network > eps] + eps)
                population = np.log10(model.values[model.values > eps])

                fig = plot_dists(population, network, '%s ROI %s (%s)' % (task, metric, level))
                fig.savefig(os.path.join(
                    OUTPUT_DIR, 'figures', '%s_%s_%s.png' % (task, level, metric)))
                plt.close(fig)

                # perform t test
                #prob = stats.ttest_1samp(network, population.mean())[1]

                # perform nonparametric test
                prob = stats.wilcoxon(
                    network, np.full(network.shape, np.median(population)))[1]

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

            path = os.path.join(OUTPUT_DIR, '%s_%s_test_results.json' % (task, level))
            with open(path, 'w') as f:
                json.dump(result, f, indent=2)



if __name__ == '__main__':
    main()
