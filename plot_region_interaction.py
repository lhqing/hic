import pandas as pd
import seaborn as sns
import subprocess
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


def transform_bedgraph(fp):
    # return format: dataframe, col1: chrom-start, col2: value
    df = pd.read_table(fp, skiprows=1, header=None, names=['chr', 'start', 'end', 'value'],
                       dtype={'chr': str, 'start': str, 'value': float})
    new_df = pd.DataFrame([df['chr'] + '-' + df['start'], df['value']], index=['region', 'value']).T
    return new_df


c_pc1_df = transform_bedgraph('/home/hliu/data/proB-C-matrix/c.PCA.25000.PC1.bedGraph')
g_pc1_df = transform_bedgraph('/home/hliu/data/proB-G-matrix/g.PCA.25000.PC1.bedGraph')
c_pc1_large_df = transform_bedgraph('/home/hliu/data/proB-C-matrix/c.PCA.500000.PC1.bedGraph')
g_pc1_large_df = transform_bedgraph('/home/hliu/data/proB-G-matrix/g.PCA.500000.PC1.bedGraph')
cg_corrdiff_df = transform_bedgraph('/home/hliu/data/compare/c_vs_g.25000.corrDiff.bedGraph')


def condense_regions(df, ratio=20, overlap=5, method='sum'):
    if overlap >= ratio:
        raise ValueError('overlap %d >= ratio %d' % (overlap, ratio))
    row_count = df.shape[0]
    if ratio > row_count:
        return df
    starts = [i for i in list(range(0, row_count, ratio - overlap)) if (i + ratio < row_count)]
    ends = [i + ratio for i in starts]
    if ends[-1] != row_count:
        starts.append(starts[-1] + ratio - overlap)
        ends.append(row_count)

    bins = list(zip(starts, ends))
    regions = []
    region_index = []
    for (start, end) in bins:
        region_df = df.iloc[start:end]
        if method == 'sum':
            regions.append(region_df.sum())
        elif method == 'mean':
            regions.append(region_df.mean())
        elif method == 'first':
            regions.append(region_df.iloc[0])
        else:
            raise ValueError('Unknown method:', method)
        region_index.append(region_df.index[0].split('-')[1] + '-' + region_df.index[-1].split('-')[1])
    condense_df = pd.DataFrame(regions, index=region_index)
    return condense_df


def read_region_matrix(fp, condense_inter=True, condense_intra=False, ratio=20, overlap=5, region_expend=None):
    df = pd.read_table(fp)
    tmp = df.iloc[:, 2:].T
    tmp.rename(columns=df.iloc[:, 1], inplace=True)
    col_chr_l = list(set(map(lambda i: i.split('-')[0], tmp.columns.tolist())))
    if len(col_chr_l) != 1:
        raise ValueError('Region not located in same chromosome:', col_chr_l)
    else:
        col_chr = col_chr_l[0]

    intra_region = [i for i in tmp.index.tolist() if col_chr == i.split('-')[0]]
    inter_region = [i for i in tmp.index.tolist() if col_chr != i.split('-')[0]]

    intra_df = tmp.loc[intra_region]
    inter_df = tmp.loc[inter_region]
    if region_expend:
        pos = os.path.split(fp)[-1].split('.')[1:4]
        pos_list = intra_df.index.map(lambda i: int(i.split('-')[1]))
        intra_df = intra_df[(pos_list > (int(pos[1]) - region_expend)) & (pos_list < (int(pos[2]) + region_expend))]

    if condense_inter:
        condense_df = inter_df.groupby(lambda i: i.split('-')[0]).apply(condense_regions,
                                                                        ratio=ratio, overlap=overlap)
        condense_df.reset_index(inplace=True)
        condense_df.rename(columns={'level_0': 'chr', 'level_1': 'region'}, inplace=True)
        condense_df.set_index(condense_df['chr'] + '-' + condense_df['region'], inplace=True)
        inter_df = condense_df.iloc[:, 2:]
    if condense_intra:
        condense_df = intra_df.groupby(lambda i: i.split('-')[0]).apply(condense_regions,
                                                                        ratio=ratio, overlap=overlap)
        condense_df.reset_index(inplace=True)
        condense_df.rename(columns={'level_0': 'chr', 'level_1': 'region'}, inplace=True)
        condense_df.set_index(condense_df['chr'] + '-' + condense_df['region'], inplace=True)
        intra_df = condense_df.iloc[:, 2:]
    return intra_df, inter_df


def call_region_matrix(chrom, start, end, tag_dirs, result_dir, prefixes,
                       resolution=25000, expend=0, region_name='', ratio=200, overlap=5,
                       region_expend=None, condense_intra=True):
    # change start end to bins
    chrom_length = {'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
                    'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
                    'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
                    'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
                    'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299,
                    'chrY': 91744698, 'chrM': 16299}
    ture_bin_start = ((start // resolution) * resolution)
    ture_bin_end = ((end // resolution + 1) * resolution)
    use_bin_start = max(ture_bin_start - expend * resolution, 0)
    use_bin_end = min(ture_bin_end + expend * resolution, chrom_length[chrom])

    matrix_fps = []
    for i in range(len(tag_dirs)):
        tag_dir = tag_dirs[i]
        prefix = prefixes[i]
        output = os.path.join(result_dir,
                              '.'.join([prefix, chrom, str(use_bin_start),
                                        str(use_bin_end), 'vsGenome',
                                        str(resolution), region_name, 'tsv']))
        matrix_fps.append(output)
        if os.path.exists(output):
            continue
        subprocess.check_call(['analyzeHiC', tag_dir, '-res', str(resolution),
                               '-chr', chrom, '-start', str(use_bin_start), '-end', str(use_bin_end),
                               '-vsGenome', '-cpu', '12'],  # , '-raw'
                              stdout=open(output, 'w'))
    output = os.path.join(result_dir,
                          '.'.join(['matrix', chrom, str(use_bin_start), str(use_bin_end),
                                    'vsGenome', region_name, 'png']))
    print(output)
    make_region_plot(matrix_fps[0], matrix_fps[1], savep=output, ratio=ratio, overlap=overlap,
                     region_expend=region_expend, condense_intra=condense_intra)
    return


def get_bedgraph_value(chrom, start, end, bedgraph_df):
    # return region dataframe
    chrom_list = bedgraph_df['region'].apply(lambda i: i.split('-')[0])
    start_list = bedgraph_df['region'].apply(lambda i: int(i.split('-')[1]))
    return bedgraph_df[(chrom_list == chrom) &
                       (start_list >= int(start)) &
                       (start_list < int(end))]


def make_region_plot(fpc, fpg, savep, c_pc1_df=c_pc1_df, g_pc1_df=g_pc1_df, cg_corrdiff_df=cg_corrdiff_df,
                     condense_inter=True, condense_intra=True, ratio=200, overlap=5,
                     region_expend=None):
    pos = os.path.split(fpc)[-1].split('.')[1:4]
    c_pc1_region = get_bedgraph_value(*pos, c_pc1_df)
    g_pc1_region = get_bedgraph_value(*pos, g_pc1_df)
    cg_corrdiff_region = get_bedgraph_value(*pos, cg_corrdiff_df)

    intra_c, _ = read_region_matrix(fpc, condense_inter, condense_intra, ratio, overlap, region_expend)
    intra_g, _ = read_region_matrix(fpg, condense_inter, condense_intra, ratio, overlap, region_expend)
    delta_intra = intra_c - intra_g

    if region_expend:
        chrom_c_pc1_df = c_pc1_df[c_pc1_df['region'].apply(lambda i: i.split('-')[0]) == pos[0]]
        chrom_cg_corrdiff_df = cg_corrdiff_df[cg_corrdiff_df['region'].apply(lambda i: i.split('-')[0]) == pos[0]]
        pos_list = chrom_c_pc1_df['region'].apply(lambda i: int(i.split('-')[1]))
        chrom_c_pc1_df = chrom_c_pc1_df[(pos_list > (int(pos[1]) - region_expend)) &
                                        (pos_list < (int(pos[2]) + region_expend))]
        pos_list = chrom_cg_corrdiff_df['region'].apply(lambda i: int(i.split('-')[1]))
        chrom_cg_corrdiff_df = chrom_cg_corrdiff_df[(pos_list > (int(pos[1]) - region_expend)) &
                                                    (pos_list < (int(pos[2]) + region_expend))]
    else:
        chrom_c_pc1_df = c_pc1_large_df[c_pc1_large_df['region'].apply(lambda i: i.split('-')[0]) == pos[0]]
        chrom_cg_corrdiff_df = cg_corrdiff_df[cg_corrdiff_df['region'].apply(lambda i: i.split('-')[0]) == pos[0]]

    fig = plt.figure()
    fig.set_size_inches(17, 10)
    gs = gridspec.GridSpec(10, 17)
    ax1 = plt.subplot(gs[:9, :5])
    ax2 = plt.subplot(gs[:9, 5:10])
    ax3 = plt.subplot(gs[:9, 10:15])
    ax4 = plt.subplot(gs[:9, 15])
    ax5 = plt.subplot(gs[:9, 16])
    ax6 = plt.subplot(gs[9, :5])
    ax7 = plt.subplot(gs[9, 5:10])
    ax8 = plt.subplot(gs[9, 10:15])

    heatmap_range = [-10, 10]
    sns.heatmap(intra_c, xticklabels=False, yticklabels=False, cbar=False,
                ax=ax1, cmap="RdBu_r", vmin=heatmap_range[0], vmax=heatmap_range[1])
    sns.heatmap(intra_g, xticklabels=False, yticklabels=False, cbar=False,
                ax=ax2, cmap="RdBu_r", vmin=heatmap_range[0], vmax=heatmap_range[1])
    sns.heatmap(delta_intra, xticklabels=False, yticklabels=False, cbar=False,
                ax=ax3, cmap="RdBu_r", vmin=heatmap_range[0], vmax=heatmap_range[1])
    bar_df = pd.DataFrame([chrom_c_pc1_df['region'].apply(lambda i: int(i.split('-')[1])).tolist(),
                           list(map(float, chrom_c_pc1_df['value'].tolist()))],
                          index=['region', 'value']).T
    sns.barplot(x='value', y='region', data=bar_df, ci=None, ax=ax4, color="DarkCyan", orient='h')
    bar_df = pd.DataFrame([chrom_cg_corrdiff_df['region'].apply(lambda i: int(i.split('-')[1])).tolist(),
                           list(map(lambda i: 1 - float(i), chrom_cg_corrdiff_df['value'].tolist()))],
                          index=['region', 'value']).T
    sns.barplot(x='value', y='region', data=bar_df, ci=None, ax=ax5, color="Salmon", orient='h')
    bar_df = pd.DataFrame([c_pc1_region['region'].apply(lambda i: int(i.split('-')[1])).tolist(),
                           list(map(float, c_pc1_region['value'].tolist()))], index=['region', 'value']).T
    sns.barplot(x='region', y='value', data=bar_df, ci=None, ax=ax6, color="DarkCyan")
    bar_df = pd.DataFrame([g_pc1_region['region'].apply(lambda i: int(i.split('-')[1])).tolist(),
                           list(map(float, g_pc1_region['value'].tolist()))], index=['region', 'value']).T
    sns.barplot(x='region', y='value', data=bar_df, ci=None, ax=ax7, color="Salmon")
    bar_df = pd.DataFrame([cg_corrdiff_region['region'].apply(lambda i: int(i.split('-')[1])).tolist(),
                           list(map(lambda i: 1 - float(i), cg_corrdiff_region['value'].tolist()))],
                          index=['region', 'value']).T
    sns.barplot(x='region', y='value', data=bar_df, ci=None, ax=ax8, color="DodgerBlue")

    ax1.set_title('Interaction of proB-C')
    ax2.set_title('Interaction of proB-G')
    ax3.set_title('Delta Interaction of proB-C vs proB-G')
    sns.despine(ax=ax4, top=True, right=True, left=False, bottom=True)
    ax4.get_yaxis().set_visible(False)
    ax4.xaxis.set_ticklabels([])
    ax4.set_xlabel('proB-C\nGenome\nPC1')
    sns.despine(ax=ax5, top=True, right=True, left=False, bottom=True)
    ax5.get_yaxis().set_visible(False)
    ax5.xaxis.set_ticklabels([])
    ax5.set_xlabel('C-G\nCorrDiff')
    sns.despine(ax=ax6, top=True, right=True, left=True, bottom=False)
    ax6.get_yaxis().set_visible(False)
    ax6.xaxis.set_ticklabels([])
    ax6.set_xlabel('proB-C Region PC1')
    sns.despine(ax=ax7, top=True, right=True, left=True, bottom=False)
    ax7.get_yaxis().set_visible(False)
    ax7.xaxis.set_ticklabels([])
    ax7.set_xlabel('proB-G Region PC1')
    sns.despine(ax=ax8, top=True, right=True, left=True, bottom=False)
    ax8.get_yaxis().set_visible(False)
    ax8.xaxis.set_ticklabels([])
    ax8.set_xlabel('Region Correlation Difference')
    plt.savefig(savep)
    plt.clf()
    return


if __name__ == '__main__':
    # chrom='chr1', start=7410679, end=7411993
    # chrom='chr2', start=101625000, end=101975000
    call_region_matrix(chrom='chr1', start=7410679, end=7411993, expend=15,
                       tag_dirs=['/home/hliu/data/proB-C-tag', '/home/hliu/data/proB-G-tag'],
                       result_dir='/home/hliu/data/result/region_matrix',
                       prefixes=['c', 'g'], region_name='testgene', ratio=3, overlap=0, region_expend=3000000,
                       condense_intra=True)



