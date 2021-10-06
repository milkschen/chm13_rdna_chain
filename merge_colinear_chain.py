'''
Merge chain records when possible
'''

import argparse
import copy
import intervaltree
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--chain', required=True,
        help='Path to the chain file'
    )
    parser.add_argument(
        '-o', '--output', default='',
        help='Path to the output merged chain file.'
    )
    args = parser.parse_args()
    return args


class ChainConst():
    CHDR_SCORE = 1
    CHDR_SOURCE = 2
    CHDR_SLEN = 3
    CHDR_SSTRAND = 4
    CHDR_SSTART = 5
    CHDR_SEND = 6
    CHDR_DEST = 7
    CHDR_DLEN = 8
    CHDR_DSTRAND = 9
    CHDR_DSTART = 10
    CHDR_DEND = 11
    CHDR_ID = 12
    C_SIZE = 0
    C_DS = 1
    C_DD = 2


class Chain(ChainConst):
    def __init__(self, fields=None):
        self.tr = intervaltree.IntervalTree()
        if not fields:
            return
        self.score = int(fields[self.CHDR_SCORE])
        self.source = fields[self.CHDR_SOURCE]
        self.slen = int(fields[self.CHDR_SLEN])
        self.sstart = int(fields[self.CHDR_SSTART])
        self.soffset = self.sstart
        self.send = int(fields[self.CHDR_SEND])
        self.sstrand = fields[self.CHDR_SSTRAND]
        assert fields[self.CHDR_SSTRAND] == '+'
        self.dest = fields[self.CHDR_DEST]
        self.dlen = int(fields[self.CHDR_DLEN])
        self.dstrand = fields[self.CHDR_DSTRAND]
        assert fields[self.CHDR_DSTRAND] == '+'
        self.dstart = int(fields[self.CHDR_DSTART])
        self.doffset = self.dstart
        self.dend = int(fields[self.CHDR_DEND])
        self.id = fields[self.CHDR_ID]
        self.ds = 0
        self.dd = 0

    def add_record_three(self, fields):
        segment_size = int(fields[self.C_SIZE])
        if segment_size > 0:
            self.tr[self.soffset : self.soffset + segment_size] = (self.doffset-self.soffset, self.ds, self.dd)
        self.ds = int(fields[self.C_DS])
        self.dd = int(fields[self.C_DD])
        self.soffset += (segment_size + self.ds)
        self.doffset += (segment_size + self.dd)

    def add_record_one(self, fields):
        segment_size = int(fields[self.C_SIZE])
        if segment_size > 0:
            self.tr[self.soffset : self.soffset + segment_size] = (self.doffset-self.soffset, self.ds, self.dd)

    def try_merge(self, c):
        # print(len(list(self.tr)))
        # print(self.tr.items())
        # print(c.tr.items())

        self_sorted_intervals = sorted(self.tr.all_intervals)
        c_sorted_intervals = sorted(c.tr.all_intervals)
        
        clone_tr = self.tr.copy()
        for c_intvl in c_sorted_intervals:
            for i_s, s_intvl in enumerate(self_sorted_intervals):
                if c_intvl.begin > s_intvl.end or c_intvl.end < s_intvl.begin:
                    continue
                else:
                    # If offsets are the same, good to merge
                    if c_intvl.data[0] == s_intvl.data[0]:
                        if c_intvl.begin < s_intvl.begin and c_intvl.end < s_intvl.end:
                            # Extend from beginning
                            # print(c_intvl)
                            # print(s_intvl)
                            start = min(s_intvl.begin, c_intvl.begin)
                            diff = s_intvl.begin - c_intvl.begin
                            data = (s_intvl.data[0], s_intvl.data[1]-diff, s_intvl.data[2]-diff)
                            updated_s_intvl = intervaltree.Interval(start, s_intvl.end, data)
                            clone_tr.remove(s_intvl)
                            c.tr.remove(c_intvl)
                            clone_tr.add(updated_s_intvl)
                        elif c_intvl.begin > s_intvl.begin and c_intvl.end > s_intvl.end:
                            # Extend from end
                            # This operation is trickier b/c the chain diffs are stored in the
                            # next interval. We need to pop the next interval and update it.
                            # print(c_intvl)
                            # print(s_intvl)
                            clone_tr.remove(s_intvl)
                            c.tr.remove(c_intvl)
                            end = max(s_intvl.end, c_intvl.end)
                            diff = c_intvl.end - s_intvl.end
                            clone_tr.add(intervaltree.Interval(s_intvl.begin, end, s_intvl[2]))
                            next_s_intvl = self_sorted_intervals[i_s+1]
                            next_data = (
                                next_s_intvl.data[0], next_s_intvl.data[1]-diff,
                                next_s_intvl.data[2]-diff)
                            updated_next_s_intvl = intervaltree.Interval(
                                next_s_intvl.begin, next_s_intvl.end, next_data)
                            clone_tr.remove(next_s_intvl)
                            clone_tr.add(updated_next_s_intvl)
                        elif c_intvl.begin < s_intvl.begina and c_intvl.end > s_intvl.end:
                            print('During merging, the new interval consumes the existing one',
                                  file=sys.stderr)
                            print('This operation is not yet supported', file=sys.stderr)
                            return False
                    else:
                        return False
        self.tr = clone_tr
        # print('\nmerged:')
        # print(self.tr.items())
        sorted_intervals = sorted(self.tr.all_intervals)

        # Update chain info
        self.score += c.score
        if sorted_intervals[0][0] < self.sstart:
            self.sstart = sorted_intervals[0][0]
            self.dstart = self.sstart + sorted_intervals[0][2]
        if sorted_intervals[-1][1] > self.send:
            self.send = sorted_intervals[-1][1]
            self.dend = self.send + sorted_intervals[-1][2]

        return True

    def print_chain(self, fout=sys.stdout):
        print((f'chain {self.score} '
               f'{self.source} {self.slen} {self.sstrand} {self.sstart} {self.send} '
               f'{self.dest} {self.dlen} {self.dstrand} {self.dstart} {self.dend} '
               f'{self.id}'),
              file=fout)
        intervals = sorted(self.tr.all_intervals)
        for i, intvl in enumerate(intervals):
            if i == 0:
                if intvl[2] != (0, 0, 0):
                    print(intervals[0].begin - self.sstart - intervals[0].data[1], intervals[0].data[1], intervals[0].data[2], file=fout)
            else:
                pintvl = intervals[i-1]
                print(pintvl.end - pintvl.begin, intvl[2][1], intvl[2][2], file=fout)
            # print(intvl[1]-intvl[0], intvl[0], intvl[1], intvl[2])
        if intvl.end == self.send and intvl.end + intvl.data[0] == self.dend:
            print(intvl.end - intvl.begin, file=fout)
        else:
            print(intvl.end - intvl.begin, self.send - intvl.end, self.dend - (intvl.end+intvl.data[0]), file=fout)
            print('0\n', file=fout)
        # print(pintvl.end - pintvl.begin, intvl[2][1], intvl[2][2], file=fout)




def merge_colinear_chain(args):
    tr_dict = {}
    f = open(args.chain, 'r')
    for line in f:
        fields = line.split()
        if len(fields) == 0:
            continue
        elif line.startswith('chain'):
            c = Chain(fields)
        elif len(fields) == 3:
            c.add_record_three(fields)
            pass
        elif len(fields) == 1:
            c.add_record_one(fields)
            if c.source in tr_dict:
                for tr in tr_dict[c.source]:
                    if not tr.try_merge(c):
                        tr_dict[c.source].append(c)
                        break
            else:
                tr_dict[c.source] = [c]
            c = None

    if args.output:
        fo = open(args.output, 'w')
    else:
        fo = sys.stdout
    for contig in tr_dict.keys():
        for tr in tr_dict[contig]:
            tr.print_chain(fo)


if __name__ == '__main__':
    args = parse_args()
    merge_colinear_chain(args)
