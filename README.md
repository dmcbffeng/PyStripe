# PyStripe
Python script for calling stripes from HiC contact maps.
- Input: .hic file
- Output: .bedpe 3D annotation

## Usage

```config
 >>> python PyStripe.py 4943.hic 4943_stripes.bedpe \
 ... -r 25000 --rg hg38 --chromosomes all \
 ... --max_distance 4000000 --min_length 500000 --min_distance 2000000 \
 ... --width 3 --merge 3  --window_size 8 \
 ... --threshold 1e-4
 ```


- '-i', '--input': .hic file path
- '-o', '--output': output bedpe path
- '-r', '--resolution': default=25000, resolution
- '--rg': default='hg38', reference genome
- '--chromosomes': default='all', chromosomes, separated by comma, e.g. 1,2,3. Can also be "all".
- '--max_distance': default=4000000, max distance off the diagonal to be calculated
- '--min_length': default=500000, minimum length of stripes
- '--min_distance': default=2000000, threshold for removing stripes too far away from the diagonal
- '--width': default=3, stripe width (# of bins)
- '--merge': default=3, merge stripes which are close to each other (# of bins)
- '--window_size': default=8, size of the window for calculating enrichment score (# of bins)
- '--threshold': default=1e-3, threshold of p values'
