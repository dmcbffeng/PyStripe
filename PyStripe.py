from utils import load_chrom_sizes, hic2txt, txt2horizontal, txt2vertical, _stripe_caller
import argparse
import numpy as np


def stripe_caller(hic_file, output_file, reference_genome='hg38', chroms='all',
                  resolution=25000, max_range=3000000,
                  interval=500000, min_length=500000, closeness=1000000,
                  stripe_width=1, merge=1,
                  window_size=5, threshold=1e-1):
    """
    :param hic_file: (str) .hic file path
    :param output_file: (str) output bedpe path
    :param reference_genome: (str) reference genome
    :param chroms: (str or list) chromosomes to calculate, default: 'all'
    :param resolution: (int) resolution
    :param max_range: (int) only consider this range off the diagonal
    :param interval: (int) minimum interval between two stripe anchors (for accelerating calculation)
    :param min_length: (int) minimum length of stripes
    :param closeness: (int) maximum distance off the diagonal
    :param stripe_width: (int) stripe width (# of bins)
    :param merge: (int) merge stripes within this range (# of bins)
    :param window_size: (int) size of the window for calculating enrichment
    :param threshold: (float) threshold of p values
    """
    chrom_lengths = load_chrom_sizes(reference_genome)
    if isinstance(chroms, str) and chroms.lower() == 'all':
        chroms = chrom_lengths.keys()
    else:
        chroms = list(chroms)

    f = open(output_file, 'w')
    f.write('chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tenrichment\t-log(pvalue)\n')
    for chrom in chroms:
        # Step 0: dump .hic file
        print(f'Dumping .hic file - {chrom}')
        hic2txt(hic_file, chrom, resolution=resolution, output='temp.txt')

        # Step 1: obtain all horizontal and vertical lines
        print(f'Loading contact maps - {chrom}')
        mat_h = txt2horizontal('temp.txt', chrom_lengths[chrom], max_range, resolution)
        mat_v = txt2vertical('temp.txt', chrom_lengths[chrom], max_range, resolution)
        # mat_h = np.load('../HSPC_10kb_KR.npy')

        print(f'Calculating stripes - {chrom} - horizontal')
        stripes_h = _stripe_caller(mat_h, max_range=max_range, resolution=resolution,
                                   interval=interval, min_length=min_length, closeness=closeness,
                                   stripe_width=stripe_width, merge=merge, window_size=window_size)
        for elm in stripes_h:
            st, ed, head, tail, enr, pval = elm
            if pval < threshold:
                x1, x2 = st * resolution, (ed + 1) * resolution
                y1, y2 = max((st + head) * resolution, x2), (ed + tail) * resolution  # avoid overlap
                f.write(f'{chrom}\t{x1}\t{x2}\t{chrom}\t{y1}\t{y2}\t0,255,0\t{enr}\t{-np.log(pval)}\n')  # green

        print(f'Calculating stripes - {chrom} - vertical')
        stripes_v = _stripe_caller(mat_v, max_range=max_range, resolution=resolution,
                                   interval=interval, min_length=min_length, closeness=closeness,
                                   stripe_width=stripe_width, merge=merge, window_size=window_size)
        for elm in stripes_v:
            st, ed, head, tail, enr, pval = elm
            if pval < threshold:
                x1, x2 = (st - tail) * resolution, (st - head) * resolution
                y1, y2 = st * resolution, (ed + 1) * resolution
                f.write(f'{chrom}\t{x1}\t{x2}\t{chrom}\t{y1}\t{y2}\t0,255,0\t{enr}\t{-np.log(pval)}\n')  # blue
    f.close()


def PyStripe(args):
    chromosomes = ['chr' + elm for elm in args.chromosomes.split(',')] if args.chromosomes.lower() != 'all' else 'all'
    stripe_caller(
        args.input,
        args.output,
        reference_genome=args.rg,
        chroms=chromosomes,
        resolution=args.resolution,
        max_range=args.max_distance,
        interval=300000,
        min_length=args.min_length,
        closeness=args.min_distance,
        stripe_width=args.width,
        merge=args.merge,
        window_size=args.window_size,
        threshold=args.threshold
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input',
        type=str,
        default=None,
        help='.hic file path'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='chr1_10kb_stripes_v2.txt',
        help='output bedpe path'
    )
    parser.add_argument(
        '-r', '--resolution',
        type=int,
        default=10000,
        help='resolution'
    )
    parser.add_argument(
        '--rg',
        type=str,
        default='hg38',
        help='reference genome'
    )
    parser.add_argument(
        '--chromosomes',
        type=str,
        default='all',
        help='chromosomes, separated by comma, e.g. 1,2,3. Can also be "all".'
    )
    parser.add_argument(
        '--max_distance',
        type=int,
        default=3000000,
        help='max distance off the diagonal to be calculated'
    )
    parser.add_argument(
        '--min_length',
        type=int,
        default=400000,
        help='minimum length of stripes'
    )
    parser.add_argument(
        '--min_distance',
        type=int,
        default=2000000,
        help='threshold for removing stripes too far away from the diagonal'
    )
    parser.add_argument(
        '--width',
        type=int,
        default=3,
        help='stripe width (# of bins)'
    )
    parser.add_argument(
        '--merge',
        type=int,
        default=3,
        help='merge stripes which are close to each other (# of bins)'
    )
    parser.add_argument(
        '--window_size',
        type=int,
        default=8,
        help='size of the window for calculating enrichment score'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=1e-3,
        help='threshold of p values'
    )
    args, _ = parser.parse_known_args()
    PyStripe(args)




