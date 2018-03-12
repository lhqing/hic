import subprocess
import os
import fire


def homer_matrix(i, o, name, resolution, cpu, bychr=False, *args):
    tag_dir = i
    out_dir = o
    if resolution == '':
        resolution = 1000000
    if cpu == '':
        cpu = 8
    if bychr == '':
        bychr = False

    result_path = '.'.join(list(map(str, [name, *(str(i).strip('-') for i in args), resolution, 'tsv'])))
    save_path = os.path.join(out_dir, result_path)
    command = ['analyzeHiC', tag_dir, '-res', str(resolution), *args, '-cpu', str(cpu)]
    if bychr:
        for chr in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                    'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                    'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                    'chr19', 'chrX', 'chrY', 'chrM']:

            chr_command = command + ['-chr', chr]
            chr_save_path = save_path[:-3] + chr + '.tsv'
            subprocess.check_call(chr_command,
                                  stdout=open(chr_save_path, 'w'),
                                  stderr=open(chr_save_path[:-3]+'err', 'w'))
    else:
        print(' '.join(map(str, command)))

        print(save_path)

        #subprocess.check_call(command,
        #                      stdout=open(save_path, 'w'),
        #                      stderr=open(save_path[:-3] + 'err', 'w'))
    return


if __name__ == '__main__':
    fire.Fire()