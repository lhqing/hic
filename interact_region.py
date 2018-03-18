import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")


def merge_interaction(df, binsize, verbose=False, eps=0.001, min_samples=5):
    x = df['PeakID(1)'].apply(lambda i: int(i.split('-')[1]) // binsize)
    y = df['PeakID(2)'].apply(lambda i: int(i.split('-')[1]) // binsize)
    coord_df = pd.DataFrame([x, y]).T
    coord_df['chr(1)'] = df['chr(1)'][0]
    coord_df['chr(2)'] = df['chr(2)'][0]

    X = np.array(list(zip(x, y)))
    X = StandardScaler().fit_transform(X)
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
    labels = db.labels_

    if verbose:
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        print('Estimated number of clusters in %s vs %s: %d' %
              (df['chr(1)'][0], df['chr(2)'][0], n_clusters_))

    coord_df['cluster_dbscan'] = labels
    return coord_df


def merge_homer_inter_table(fp, binsize=20000, eps=0.001, min_samples=5):
    df = pd.read_table(fp, index_col=0)
    cluster = df.groupby(['chr(1)', 'chr(2)']).apply(merge_interaction, binsize=binsize,
                                                     eps=eps, min_samples=min_samples)
    return cluster


