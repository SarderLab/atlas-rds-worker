#!/usr/bin/env python3

import argparse
import itertools
import random
import sys
import time

import h5py
import numpy as np

# Inspired from https://stackoverflow.com/a/43374773


def scan_node(src, dest=None, analyze=False, showattrs=False, convert=None, sample=None, indent=0):  # noqa
    print('%s%s' % ('  ' * indent, src.name))
    for ak in src.attrs:
        if showattrs:
            print('%s:%s: %r' % ('  ' * indent, ak, src.attrs[ak]))
        if dest:
            dest.attrs[ak] = src.attrs[ak]
    for k, v in src.items():
        if isinstance(v, h5py.Dataset):
            print('%s - %s %s %r %r %s' % (
                '  ' * (indent + 1), v.name, v.dtype, v.shape,
                v.chunks, v.compression))
            if showattrs:
                for ak in v.attrs:
                    print('%s   :%s: %r' % ('  ' * (indent + 1), ak, v.attrs[ak]))
            lasttime = time.time()
            if v.dtype.kind in {'f', 'i'} and analyze:
                minv = maxv = None
                sumv = 0
                for coor in itertools.product(*(
                        range(0, v.shape[idx], v.chunks[idx]) for idx in range(len(v.shape)))):
                    field = tuple(
                        slice(coor[idx], min(coor[idx] + v.chunks[idx], v.shape[idx]))
                        for idx in range(len(v.shape)))
                    part = v[field]
                    if minv is None:
                        minv = np.amin(part)
                        maxv = np.amax(part)
                    else:
                        minv = min(minv, np.amin(part))
                        maxv = max(maxv, np.amax(part))
                    if part.dtype == np.float16:
                        part = part.astype(np.float32)
                    sumv += part.sum()
                avgv = sumv / v.size
                print('%s   [%g,%g] %g' % (
                    '  ' * (indent + 1), minv, maxv, avgv))
            if dest:
                conv = convert and (v.dtype == np.float64 or (
                    v.dtype == np.float32 and convert == 'float16'))
                if conv:
                    conv = np.float32 if convert == 'float32' or max(
                        abs(minv), maxv) >= 65504 else np.float16
                    if conv == v.dtype:
                        conv = False
                vshape = v.shape
                if sample is not None:
                    if v.shape[-1] == sample.shape[-1]:
                        vshape = tuple(list(vshape)[:-1] + [sample.sum()])
                    elif v.shape[0] == sample.shape[-1]:
                        vshape = tuple([sample.sum()] + list(vshape)[1:])
                if conv:
                    destv = dest.create_dataset(
                        k, shape=vshape,
                        dtype=conv,
                        chunks=True, fillvalue=0,
                        compression='gzip', compression_opts=9, shuffle=True)
                else:
                    destv = dest.create_dataset(
                        k, shape=vshape,
                        dtype=v.dtype,
                        chunks=True, fillvalue=v.fillvalue,
                        compression='gzip', compression_opts=9, shuffle=v.shuffle)
                for ak in v.attrs:
                    destv.attrs[ak] = v.attrs[ak]
                steps = len(list(itertools.product(*(
                    range(0, v.shape[idx], destv.chunks[idx])
                    for idx in range(len(v.shape))))))
                skip = 0
                if sample is not None:
                    if v.shape[-1] == sample.shape[-1]:
                        if len(v.shape) == 1:
                            v = v[sample]
                        else:
                            v = v[..., sample]
                    elif v.shape[0] == sample.shape[-1]:
                        v = v[sample, ...]
                for cidx, coor in enumerate(itertools.product(*(
                        range(0, v.shape[idx], destv.chunks[idx])
                        for idx in range(len(v.shape))))):
                    if time.time() - lasttime > 10:
                        sys.stdout.write('  %5.2f%% %r %r %r\r' % (
                            100.0 * cidx / steps, coor, v.shape, destv.chunks))
                        sys.stdout.flush()
                        lasttime = time.time()
                    field = tuple(
                        slice(coor[idx], min(coor[idx] + destv.chunks[idx], v.shape[idx]))
                        for idx in range(len(v.shape)))
                    part = v[field]
                    if conv:
                        if not part.any():
                            skip += 1
                            continue
                        part = part.astype(conv)
                    destv[field] = part
                print('%s > %s %s %r %r %s%s' % (
                    '  ' * (indent + 1), destv.name, destv.dtype, destv.shape,
                    destv.chunks, destv.compression,
                    ' %d' % skip if skip else ''))

        elif isinstance(v, h5py.Group):
            destv = None
            if dest:
                destv = dest.create_group(k)
            scan_node(v, destv, analyze, showattrs, convert, sample, indent=indent + 1)
    if (dest and sample is not None and 'dims' in src.attrs and
            src.attrs['dims'][-1] == sample.shape[-1]):
        dims = src.attrs['dims']
        dims[-1] = sample.sum()
        dest.attrs['dims'] = dims
        sdata = src[src.name + '/data']
        sind = src[src.name + '/indices']
        sindptr = src[src.name + '/indptr']
        print('%s>... %s (%d, %d)' % ('  ' * indent, src.name, sdata.shape[-1], sindptr.shape[-1]))
        ddata = dest[dest.name + '/data']
        dind = dest[dest.name + '/indices']
        dindptr = dest[dest.name + '/indptr']
        dpos = 0
        didx = 0
        for sidx, samp in enumerate(sample):
            if not samp:
                continue
            spos = sindptr[sidx]
            slen = sindptr[sidx + 1] - spos
            dindptr[didx] = dpos
            ddata[dpos:dpos + slen] = sdata[spos:spos + slen]
            dind[dpos:dpos + slen] = sind[spos:spos + slen]
            didx += 1
            dpos += slen
        dindptr[didx] = dpos
        ddata.resize((dpos, ))
        dind.resize((dpos, ))
        dindptr.resize((didx + 1, ))
        print('%s...> %s (%d, %d)' % ('  ' * indent, dest.name, ddata.shape[-1], dindptr.shape[-1]))


def scan_hdf5(path, analyze=False, showattrs=False, outpath=None, convert=None,
              sample=None):
    if convert:
        analyze = True
    with h5py.File(path, 'r') as fptr:
        fptr2 = None
        if outpath:
            fptr2 = h5py.File(outpath, 'w')
        scan_node(fptr, fptr2, analyze, showattrs, convert, sample)


def command():
    parser = argparse.ArgumentParser(
        description='Scan and optionally subsample an hdf5 file and report on '
        'its groups, datasets, and attributes.  Optionally report minimum, '
        'maximum, and average values for datasets with integer or float '
        'datatypes.  Optionally rewrite the file with lower precision float '
        'datasets or subsampled values.')
    parser.add_argument(
        'source', type=str, help='Source file to read and analyze.')
    parser.add_argument(
        '--analyze', '-s', action='store_true',
        help='Analyze the min/max/average of datasets.')
    parser.add_argument(
        '--attrs', '-k', action='store_true',
        help='Show attributes on groups and datasets.')
    parser.add_argument(
        '--dest', help='Write a new output file')
    parser.add_argument(
        '--convert', choices=('float16', 'float32'),
        help='Reduce the precision of the output file.')
    parser.add_argument(
        '--category', '--cat',
        help='This is the path of a 1-d array that contains the unique values '
        'to use as sampling categories (e.g., /meta.data/subclass.l1).  For '
        'random sampling, this is only used to determine the length of the '
        'data that gets sampled.')
    parser.add_argument(
        '--sample', type=float,
        help='Reduce the main array to the specified count.  This must be a '
        'value in the range [0-1].')
    parser.add_argument(
        '--method', default='random', choices=('random', 'equalize', 'capped'),
        help='Sampling method to use when reducing the main array.  Random '
        'picks each entry with equal probability.  Equalize and capped keep '
        'an approximately equal number for each category.  If some categories '
        'have few enough values such that all of their values are used, '
        'equalize will increase other categories so that the sample '
        'proportion of the whole atlas is the specified value; capped will '
        'not do this increase and therefore end up with a smaller atlas.')
    opts = parser.parse_args()
    if opts.sample and opts.category:
        with h5py.File(opts.source, 'r') as fptr:
            v = fptr[opts.category]
            if len(v.shape) != 1:
                raise Exception('Surprising shape')
            sample = np.random.rand(v.shape[0]) < opts.sample
            if opts.method in {'equalize', 'capped'}:
                sampleset = dict(zip(*np.unique(v, return_counts=True)))
                avgcount = v.shape[0] / len(sampleset) * opts.sample
                # Adjust the proportions by a factor to account for some
                # categories being fully sampled
                factor = 1
                for _ in range(5 if opts.method == 'equalize' else 1):
                    proportions = {k: min(1, avgcount / val * factor)
                                   for k, val in sampleset.items()}
                    factor *= opts.sample * v.shape[0] / sum([
                        sampleset[k] * proportions[k] for k in sampleset])
                for i in range(v.shape[0]):
                    sample[i] = random.random() < proportions[v[i]]
        print('Sampling %d down to %d (%5.3f)' % (
            sample.shape[0], sample.sum(), sample.sum() / sample.shape[0]))
    scan_hdf5(opts.source, opts.analyze, opts.attrs, opts.dest, opts.convert,
              sample)


if __name__ == '__main__':
    command()
