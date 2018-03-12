import pandas as pd


def filter_corrDiff(corr_diff, global_threshold=0.15, local_threshold=0.08, gap_tolerance=1, length_threshold=3,
                    corr_transfer=lambda i: 1 - abs(i)):
    tmp_corr_diff = corr_diff.apply(corr_transfer)
    remain_item = []

    temp_list = []
    gap = 0
    cur_length = 0
    pass_global = False
    for i, corr in tmp_corr_diff.iteritems():
        if corr > global_threshold:
            temp_list.append(i)
            pass_global = True
        elif corr > local_threshold:
            temp_list.append(i)
        else:
            if gap < gap_tolerance:
                temp_list.append(i)
                gap += 1
            else:
                if pass_global and (len(temp_list) >= length_threshold):
                    remain_item += temp_list

                temp_list = []
                gap = 0
                cur_length = 0
                pass_global = False
    return corr_diff.loc[remain_item]


def condense_regions(df, ratio=20, overlap=5):
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
        regions.append(region_df.sum())
        region_index.append(region_df.index[0].split('-')[1] + '-' + region_df.index[-1].split('-')[1])
    condense_df = pd.DataFrame(regions, index=region_index)
    return condense_df


def read_region_matrix(fp, condense_inter=True, ratio=20, overlap=5):
    df = pd.read_table(fp)
    tmp = df.iloc[:, 2:].T
    tmp.rename(columns=df.iloc[:, 1], inplace=True)
    col_chr_l = list(set(map(lambda i: i.split('-')[0], tmp.columns.tolist())))
    if len(col_chr_l) != 1:
        raise ValueError('Region not located in same chromosome:', col_chr_l)
    else:
        col_chr = col_chr_l[0]
    
    intra_region = [i for i in tmp.index.tolist() if col_chr in i]
    inter_region = [i for i in tmp.index.tolist() if col_chr not in i]
    
    intra_df = tmp.loc[intra_region]
    inter_df = tmp.loc[inter_region]
    
    if condense_inter:
        condense_df = inter_df.groupby(lambda i: i.split('-')[0]).apply(condense_regions,
                                                                        ratio=ratio, overlap=overlap)
        condense_df.reset_index(inplace=True)
        condense_df.rename(columns={'level_0': 'chr', 'level_1': 'region'}, inplace=True)
        condense_df.set_index(condense_df['chr'] + '-' + condense_df['region'], inplace=True)
        condensed_inter_df = condense_df.iloc[:, 2:]
        return intra_df, condensed_inter_df
    else:
        return intra_df, inter_df