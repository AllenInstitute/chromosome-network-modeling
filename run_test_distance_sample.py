from multiprocessing import Pool, cpu_count
from collections import Counter
from functools import partial
import os
import json
import binascii

from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from mcmodels.core import VoxelModelCache

MANIFEST_FILE = os.path.join(os.path.expanduser('~'), 'connectivity',
                             'voxel_model_manifest.json')
ROI_PATH = 'behavior_rois.json'
MODEL_DIR = 'model'
OUTPUT_DIR = 'output'
FIGURE_DIR = os.path.join(OUTPUT_DIR, 'figures')
DISTANCES_PATH = os.path.join(MODEL_DIR, 'distances.csv')
WEIGHTS_PATH = os.path.join(MODEL_DIR, 'normalized_connection_density.csv')

SUMMARY_STRUCTURES_SET_ID = 687527945

N_SAMPLES = 1000
TASKS = ('oft', '3-chamber')


def sample_network(weights, population, network, max_iter=10000, thresh=1, alpha=0.01):
    def get_random_seed():
        return int(binascii.b2a_hex(os.urandom(4)), 16)

    rng = np.random.RandomState(get_random_seed())

    n = population.shape[0] # NOTE: should be square
    k = network.shape[0]

    tri_idx = np.tril_indices_from(network)
    network_ = network[np.tril_indices_from(network)]

    for _ in range(max_iter):
        idx = rng.choice(n, size=k, replace=False)
        sample = population[np.ix_(idx, idx)]
        D, pvalue = stats.ks_2samp(sample[tri_idx], network_)

        if D < thresh and pvalue < alpha:
            return weights[np.ix_(idx, idx)]

    raise RuntimeError('a suitable sample was not found in %d iterations' % max_iter)


def get_ss(tree, acs):
    ss = [s['acronym'] for s in tree.get_structures_by_set_id([SUMMARY_STRUCTURES_SET_ID])]
    return list(set(ss).intersection(acs))


def get_p_value(samples, test):
    samples = np.log10(samples)
    test = np.log10(test)

    return stats.ttest_1samp(samples, test)[1]

def plot(samples, test, task):
    samples = np.log10(samples)
    test = np.log10(test)

    fig, ax = plt.subplots()
    sns.set_style('white')
    sns.distplot(samples, ax=ax, color=sns.xkcd_rgb['dusky blue'])
    plt.plot([test, test], [0, 0.9 * ax.get_ylim()[1]], 'r-')

    # remove y
    sns.despine(left=True)
    ax.set_xlim(-6.5, -3.5)
    ax.set_yticks([])

    # label
    ax.set_xlabel('Log normalized connection density')
    ax.set_title(task, loc='left')

    return fig


def main():
    # read in ROIs
    with open(ROI_PATH, 'r') as f:
        task_rois = json.load(f)

    # initialize cache object
    cache = VoxelModelCache(manifest_file=MANIFEST_FILE)
    rs = cache.get_reference_space()
    rs.remove_unassigned(update_self=True)
    tree = rs.structure_tree

    # load distances
    distances = pd.read_csv(DISTANCES_PATH, index_col=0)
    model = pd.read_csv(WEIGHTS_PATH, index_col=0)

    # get level
    level = get_ss(tree, model.index.values)

    for task in TASKS:
        net = set(level).intersection(task_rois[task])

        weights = model.loc[level, level].values
        population = distances.loc[level, level].values
        network = distances.loc[net, net].values

        valid = np.ix_(~np.isnan(weights).all(axis=1), ~np.isnan(weights).all(axis=0))
        weights = weights[valid]
        population = population[valid]

        # run
        pool = Pool(processes=cpu_count())
        func = partial(sample_network, weights, population, network,
                       max_iter=int(1e6), thresh=0.25)
        mapper = [pool.apply_async(func) for _ in range(N_SAMPLES)]
        samples = np.array([np.median(f.get()) for f in mapper])

        test = np.median(model.loc[net, net].values)
        print("unique graphs: ", Counter(np.unique(samples, return_counts=True)[1]))

        # compute stats
        p_value = get_p_value(samples, test)
        pctl = 100 * np.sum(samples < test) / samples.size

        # plot
        fig = plot(samples, test, task)
        fig.savefig(os.path.join(FIGURE_DIR, '%s_distance_analysis.svg' % task))
        plt.close(fig)

        # save results
        result = dict(task=task,
                      rois=task_rois[task],
                      p_value=p_value,
                      percentile=pctl)
        with open(os.path.join(OUTPUT_DIR, '%s_results' % task), 'w') as f:
            json.dump(result, f, indent=2)

        # save ROI graph
        model.loc[level, level].to_csv(
            os.path.join(MODEL_DIR, '%s_normalized_connection_density.csv') % task)


if __name__ == '__main__':
    main()
