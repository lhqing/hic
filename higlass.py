import cooler
import pandas as pd
import subprocess
import os
import shlex
import glob


def homer2cooler(fp, save_dir, bin_size, bin_file, chrom_size_file):
    df = pd.read_table(fp, index_col=[0, 1])
    file_name = os.path.splitext(os.path.split(fp)[-1])[0]
    save_path = os.path.join(save_dir, file_name + '.cool_in')
    fs = open(save_path, 'w')
    for i, row in df.iterrows():
        chrom1, start1 = i[0].split('-')
        end1 = int(start1) + bin_size

        for j, v in row.iteritems():
            chrom2, start2 = j.split('-')
            end2 = int(start2) + bin_size
            if v > 0:
                line_list = [chrom1, start1, str(end1), chrom2, start2, str(end2), '1', str(v)]
                # 1 for count, value actually store in the last column weight... That's stupid
                fs.write('\t'.join(line_list) + '\n')
    fs.close()

    # cooler csort
    subprocess.check_call(['cooler', 'csort', '-c1', '1', '-p1', '2',
                           '-c2', '4', '-p2', '5', '-o', save_path + '.sorted',
                           save_path, chrom_size_file])
    # cooler cload
    subprocess.check_call(['cooler', 'cload', 'pairix', bin_file, save_path + '.sorted',
                           os.path.join(save_dir, file_name + '.cool')])
    # cooler zoomify
    subprocess.check_call(['cooler', 'zoomify', '--no-balance',
                           os.path.join(save_dir, file_name + '.cool')])
    subprocess.check_call(['rm', '-f', os.path.join(save_dir, file_name + '.cool_in*')])
    return


def ingest_cool(fp, container, map_dir):
    fl = glob.glob(fp)
    for f in fl:
        print(f)
        subprocess.check_call(['docker', 'exec', container, 'python', 'higlass-server/manage.py', 'ingest_tileset',
                               '--filename', os.path.join(map_dir, os.path.split(f)[1]),
                               '--datatype', 'matrix', '--filetype', 'cooler'])
    return


if __name__ == '__main__':
    #for f in glob.glob('/home/hliu/data/proB-*-matrix/*.25000.chr*.tsv'):
    #    homer2cooler(f,
    #                 save_dir='/home/hliu/hg-data', bin_size=25000,
    #                 bin_file='/home/hliu/ref/bins.25kb.bed',
    #                 chrom_size_file='/home/hliu/ref/mm10_chromsizes.tsv')

    ingest_cool('/home/hliu/hg-data/*.mcool', container='higlass-container-hq', map_dir='/data')

