import os
import numpy as np
from scipy.stats import kruskal


def hic2txt(hic_file, ch, resolution=25000, output='temp.txt'):
    """
    Dump .hic file into contact lists
    :param hic_file: (str) .hic file path
    :param ch: (str) chromosome
    :param resolution: (int) resolution to use
    :param output: (str) temporary output path
    """
    juicer = 'juicer_tools_1.11.04_jcuda.0.8.jar'
    cmd = f'java -jar {juicer} dump observed KR {hic_file} {ch} {ch} BP {resolution} {output}'
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
    # st, ed = max(distance_range[0], window_size), min(distance_range[1], mat.shape[1] - window_size)
    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width

    new_mat = np.zeros((distance_range[1] - distance_range[0],))
    for j in range(distance_range[0], distance_range[1]):
        if j < window_size + half or j >= mat.shape[1] - window_size - half:
            continue
        y = j - distance_range[0]
        line_min = min(np.mean(mat[x1:x2, j-window_size-half:j-half]),
                       np.mean(mat[x1:x2, j+1+half:j+window_size+half+1]))
        neighbor_mean = max(np.mean(mat[idx-window_size:x1, j-window_size-half:j+window_size+half+1]),
                            np.mean(mat[x2+1:idx+window_size+1, j-window_size-half:j+window_size+half+1]))
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


def merge_positions(lst, merge_range):
    def _merge(small_lst):
        # print(small_lst)
        # it is wrong!
        merged = []  # [st, ed, head, tail, score]
        for elm in small_lst:
            found = False
            for j, m in enumerate(merged):
                if not (elm[2] > m[3] or elm[3] < m[2]):
                    merged[j][1] = elm[1]
                    merged[j][2] = min(merged[j][2], elm[2])
                    merged[j][3] = max(merged[j][3], elm[3])
                    merged[j][4] = (merged[j][4] + elm[4]) / 2
                    found = True
                    break
            if not found:
                merged.append(elm)
        if len(merged) > 1:
            merged.sort(key=lambda x: x[4], reverse=True)
        # print(merged[0])
        return merged[0]

    new_lst = []
    temp = []
    for i, (idx, head, tail, score) in enumerate(lst):
        if i != len(lst) - 1 and lst[i + 1][0] - idx < merge_range:
            temp.append([idx, idx, head, tail, score])
        else:
            if len(temp) != 0:
                temp.append([idx, idx, head, tail, score])
                new_lst.append(_merge(temp))
                temp = []
            else:
                new_lst.append([idx, idx, head, tail, score])
    return new_lst


def stat_test(mat, st, ed, line_width, head, tail, window_size):
    half = int(line_width // 2)
    x1, x2 = st - half, ed + half + 1
    r1 = mat[x1-window_size:x1, head:tail].flatten()
    r2 = mat[x2:x2+window_size, head:tail].flatten()
    r = mat[x1:x2, head:tail].flatten()

    t1, p1 = kruskal(r, r1)
    t2, p2 = kruskal(r, r2)
    return max(p1, p2)


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

    # Step 3: find the accurate range of stripe
    print(' Finding the spanning range for each stripe...')
    all_positions = []
    lst = sorted(positions.keys())
    for i, idx in enumerate(lst):
        # print(i, idx)
        if idx <= window_size or idx >= mat.shape[0] - window_size:
            continue
        arr = enrichment_score(mat, idx, line_width=stripe_width,
                               distance_range=(0, max_range // resolution),
                               window_size=window_size)
        head, tail, _max = find_max_slice(arr)
        all_positions.append((idx, head, tail, _max))

    # Step 4: Merging
    print(' Merging...')
    all_positions = merge_positions(all_positions, merge)
    print(len(all_positions))

    new_positions = []
    for elm in all_positions:
        # print(elm, end=' ')
        if (elm[3] - elm[2]) * resolution >= min_length and elm[2] * resolution <= closeness:
            # print(True)
            new_positions.append(elm)
        else:
            # print(False)
            pass
    print(len(new_positions))

    # Step 5: Statistical test
    results = []
    print(' Statistical Tests...')
    for elm in new_positions:
        [st, ed, head, tail, score] = elm
        p = stat_test(mat, st, ed, stripe_width, head, tail, window_size)
        # print(idx * resolution, p)
        results.append([st, ed, head, tail, score, p])
    return results



