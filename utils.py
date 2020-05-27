import os
import numpy as np
from scipy.stats import ttest_ind


def hic2txt(hic_file, ch, resolution=25000, output='temp.txt'):
    """
    Dump .hic file into contact lists
    :param hic_file: (str) .hic file path
    :param ch: (str) chromosome
    :param resolution: (int) resolution to use
    :param output: (str) temporary output path
    """
    juicer = 'juicer_tools_1.11.04_jcuda.0.8.jar'
    cmd = f'java -jar {juicer} dump oe KR {hic_file} {ch} {ch} BP {resolution} {output}'
    os.system(cmd)


def load_chrom_sizes(reference_genome):
    """
    Load chromosome sizes for a reference genome
    """
    my_path = os.path.abspath(os.path.dirname(__file__))
    f = open(os.path.join(my_path, reference_genome + '.chrom.sizes'))
    lengths = {}
    for line in f:
        [ch, l] = line.strip().split()
        lengths[ch] = int(l)
    return lengths


def txt2matrix(txt, length, resolution=25000):
    """
        Convert contact list into a numpy matrix, size: N * N

        :param txt: str, path of input .txt file
        :param length: chromosome length
        :param resolution: int, default: 25000
        """
    f = open(txt)
    n_bins = length // resolution + 1
    mat = np.zeros((n_bins, n_bins))

    for line in f:
        p1, p2, v = line.strip().split()
        if v == 'NaN':
            continue
        p1, p2, v = int(p1), int(p2), float(v)

        if max(p1, p2) >= n_bins * resolution:
            continue

        mat[p1 // resolution, p2 // resolution] += v
        if p1 // resolution != p2 // resolution:
            mat[p2 // resolution, p1 // resolution] += v

    return mat


def txt2horizontal(txt, length, max_range, resolution=25000):
    """
        :param txt: str, path of input .txt file
        :param length: chromosome length
        :param max_range: int, max distance
        :param resolution: int, default: 25000
        """
    assert max_range % resolution == 0
    f = open(txt)
    n_bins = length // resolution + 1
    rg = max_range // resolution
    mat = np.zeros((n_bins, rg))
    for line in f:
        p1, p2, v = line.strip().split()
        if v == 'NaN':
            continue
        p1, p2, v = int(p1), int(p2), float(v)
        if max(p1, p2) >= n_bins * resolution:
            continue
        if p1 > p2:
            p1, p2 = p2, p1
        p1, p2 = p1 // resolution, p2 // resolution
        if p2 - p1 >= rg:
            continue
        mat[p1, p2 - p1] += v
    return mat


def txt2vertical(txt, length, max_range, resolution=25000):
    """
        :param txt: str, path of input .txt file
        :param length: chromosome length
        :param max_range: int, max distance
        :param resolution: int, default: 25000
        """
    assert max_range % resolution == 0
    f = open(txt)
    n_bins = length // resolution + 1
    rg = max_range // resolution
    mat = np.zeros((n_bins, rg))
    for line in f:
        p1, p2, v = line.strip().split()
        if v == 'NaN':
            continue
        p1, p2, v = int(p1), int(p2), float(v)
        if max(p1, p2) >= n_bins * resolution:
            continue
        if p1 > p2:
            p1, p2 = p2, p1
        p1, p2 = p1 // resolution, p2 // resolution
        if p2 - p1 >= rg:
            continue
        mat[p2, p2 - p1] += v
    return mat


def pick_max_positions(mat, interval=500000, distance_range=(500000, 1000000), resolution=25000,
                       line_width=1, window_size=4):
    assert interval % resolution == 0
    assert distance_range[0] % resolution == 0
    assert distance_range[1] % resolution == 0

    st, ed = distance_range[0] // resolution, distance_range[1] // resolution
    size = interval // resolution
    length = mat.shape[0]
    stats = np.sum(mat[:, st:ed], axis=1)
    all_pos = []
    for i in range(0, length, size):
        region = stats[i: min(i + size, length)]
        idx = int(np.argmax(region) + i)
        # print(idx, window_size, mat.shape[0] - window_size)

        if idx < window_size or idx >= mat.shape[0] - window_size:
            continue

        previous = stats[max(0, idx - size): idx-1]
        later = stats[idx + 2: min(idx + size + 1, length)]
        # print(stats[idx], np.max(previous), np.max(later))

        if stats[idx] > np.max(previous) and stats[idx] > np.max(later):
            # print(idx)
            check = enrichment_score(mat, idx, line_width,
                                     (st, ed), window_size)
            if np.sum(check) > 0:
                all_pos.append(idx)
    return all_pos


def enrichment_score(mat, idx, line_width=1, distance_range=(20, 40), window_size=4):
    st, ed = max(distance_range[0], window_size), min(distance_range[1], mat.shape[1] - window_size)
    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width

    new_mat = np.zeros((ed - st,))
    for j in range(st, ed):
        y = j - st
        line_min = min(np.mean(mat[x1:x2, j-window_size:j-1]), np.mean(mat[x1:x2, j+2:j+window_size+1]))
        neighbor_mean = max(np.mean(mat[idx-window_size:x1-1, j-window_size:j+window_size+1]),
                            np.mean(mat[x2+1:idx+window_size+1, j-window_size:j+window_size+1]))
        new_mat[y] = line_min - neighbor_mean
    return new_mat


def find_max_slice(arr):
    _max, head, tail = 0, 0, 0
    _max_ending, h, t = 0, 0, 0
    i = 0
    while i < len(arr):
        _max_ending = _max_ending + arr[i]
        if _max_ending < 0:
            h, t = i + 1, i + 1
            _max_ending = 0
        else:
            t = i + 1
        if _max_ending > _max:
            head, tail, _max = h, t, _max_ending
        i += 1
    return head, tail, _max


def stat_test(mat, idx, line_width, head, tail, window_size):
    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width
    r1 = mat[idx-window_size:x1-1, head:tail].flatten()
    r2 = mat[x2+1:idx+window_size+line_width+1, head:tail].flatten()
    r = mat[x1:x2, head:tail].flatten()

    t1, p1 = ttest_ind(r, r1)
    t2, p2 = ttest_ind(r, r2)
    return max(p1 / 2, p2 / 2)


def _stripe_caller(mat, max_range=3000000, resolution=25000,
                   interval=200000, min_length=500000, closeness=1000000,
                   stripe_width=1, merge=1, window_size=5):
    # Step 2: for different distance ranges pick the "local maximum" positions
    print(' Finding local maximum for different contact distances...')
    positions = {}
    for dis in range(0, max_range - min_length + 1, min_length // 2):
        print(f'  {dis}-{dis + interval}')
        distance_range = (dis, dis + interval)
        pos_h = pick_max_positions(mat, interval=interval, distance_range=distance_range,
                                   resolution=resolution, line_width=stripe_width, window_size=window_size)
        for p in pos_h:
            if p not in positions:
                positions[p] = []
            positions[p].append(distance_range)
    print(len(positions))

    # Step 3: merge nearby positions and find the accurate range of stripe
    print(' Finding the spanning range for each stripe...')
    all_positions = {}
    lst = sorted(positions.keys())
    ppp, temp, N = 0, None, 0
    for i, idx in enumerate(lst):
        if idx <= window_size or idx >= mat.shape[0] - window_size:
            continue
        arr = enrichment_score(mat, idx, line_width=stripe_width,
                               distance_range=(0, max_range),
                               window_size=window_size)
        if i < len(lst) - 1 and lst[i + 1] - idx <= merge:
            if N != 0:
                temp += arr
                N += 1
                ppp += idx
            else:
                temp = arr
                N = 1
                ppp = idx
        else:
            if N != 0:
                arr = (arr + temp) / (N + 1)
                pos = int((ppp + idx) // (N + 1))
                ppp, temp, N = 0, None, 0
            else:
                pos = idx
            head, tail, _max = find_max_slice(arr)
            if (tail - head) * resolution >= min_length and head * resolution < closeness:
                all_positions[pos] = (head, tail, _max)

    # Step 4: Statistical test
    print(' Statistical Tests...')
    for idx in all_positions:
        st, ed, mx = all_positions[idx]
        p = stat_test(mat, idx, stripe_width, st, ed, window_size)
        all_positions[idx] = (st, ed, mx, p)
    return all_positions



