import argparse
import intervaltree
import pysam
import re
import sys

# PAF fields
QNAME = 0
QLEN = 1
QSTART = 2
QEND = 3
STRAND = 4
TNAME = 5
TLEN = 6
TSTART = 7
TEND = 8
NRES = 9
ALN = 10
MAPQ = 11


def parse_args():
    # START_CHAIN_IDX = 25
    # RDNA_V1_0_PAF = 'rDNA_to_v1.0.paf'
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p', '--input_paf', default='rDNA_to_v1.0.paf',
        help='Path to the PAF file where rDNA patches are aligned to an unpolished reference.'
    )
    parser.add_argument(
        '-ic', '--start_chain_idx', default=25, type=int,
        help='Starting index of the output chain records.'
    )
    parser.add_argument(
        '-a', '--rdna_sam', default='rDNA_to_v1.1.sam',
        help='Path to the SAM file where rDNA patches are aligned to the polished reference.'
    )
    parser.add_argument(
        '-o', '--output_chain', default='',
        help='Path to the output rDNA chain file.'
    )
    args = parser.parse_args()
    return args


''' Retrieve CIGAR info from a PAF record '''
def paf_retrieve_cigar(p: list) -> list:
    cg_tag = 'cg:Z:'
    for i in range(MAPQ, len(p)):
        if p[i].startswith(cg_tag):
            return p[i][len(cg_tag):]


''' Read a PAF file '''
def read_paf(input_paf) -> list:
    paf = []
    with open(input_paf, 'r') as f:
        for line in f:
            paf.append(line.split('\t'))

    query_trees = {}
    target_trees = {}

    paf_fields = [QNAME, QSTART, QEND, TNAME, TSTART, TEND, NRES, MAPQ]
    paf.sort(key=lambda x: int(x[NRES]), reverse=True)
    paf_final = []
    for p in paf:
        if p[QNAME].split(':')[0] != p[TNAME]:
            continue
        s = ''
        for i in paf_fields:
            s += p[i]
            s += '\t'
        cg = paf_retrieve_cigar(p)
        s += f'{cg}\t'
        qname = p[QNAME]
        qstart = int(p[QSTART])
        qend = int(p[QEND])
        num_res = int(p[NRES])
        if qname in query_trees:
            # We require PAF records to be sorted by NRES, so if there's an overlap
            # we exclude the new record.
            if query_trees[qname].overlaps(qstart, qend):
                continue
            else:
                query_trees[qname][qstart: qend] = 1
        else:
            query_trees[qname] = intervaltree.IntervalTree()
            query_trees[qname][qstart: qend] = 1

        tname = p[TNAME]
        tstart = int(p[TSTART])
        tend = int(p[TEND])
        if tname in target_trees:
            if target_trees[tname].overlaps(tstart, tend):
                continue
            else:
                target_trees[tname][tstart: tend] = 1
        else:
            target_trees[tname] = intervaltree.IntervalTree()
            target_trees[tname][tstart: tend] = 1
            # target_trees[tname].add((tstart, tend))

        # print(s.rstrip(), file=sys.stderr)
        if p[STRAND] != '+':
            print(p, file=sys.stderr)
            print('reverse alignments are not yet supported', file=sys.stderr)
            exit(1)
        paf_final.append(p)
    return paf_final


def rdna_paf_to_chain(args):
    paf_final = read_paf(args.input_paf)

    # Load offsets wrt the source reference
    # Positions are 0-based
    rdna_offsets = {}
    rdna_chrlen = {}
    with pysam.AlignmentFile(args.rdna_sam) as f:
        for r in f:
            if not r.is_supplementary and not r.is_secondary:
                if r.query_name not in rdna_offsets:
                    if r.query_name.split(':')[0] == r.reference_name:
                        rdna_offsets[r.query_name] = r.pos
                        rdna_chrlen[r.reference_name] = f.get_reference_length(r.reference_name)

    # print('v1.1 rdna offsets', rdna_offsets, file=sys.stderr)
    # print('v1.1 chr lengths', rdna_chrlen, file=sys.stderr)
    # print('', file=sys.stderr)

    # Write to chain
    if args.output_chain == '':
        f_out = sys.stdout
    else:
        f_out = open(args.output_chain, 'w')

    re_cigar = re.compile('[MIDS]+')
    for idx_p, p in enumerate(paf_final):
        qname = p[QNAME]
        chr_qname = qname.split(':')[0]
        qstart = int(p[QSTART]) + rdna_offsets[qname]
        qend = int(p[QEND]) + rdna_offsets[qname]
        print((f'chain {p[NRES]} '
               f'{chr_qname} {rdna_chrlen[chr_qname]} {p[STRAND]} {qstart} {qend} '
               f'{p[TNAME]} {p[TLEN]} + {p[TSTART]} {p[TEND]} '
               f'{args.start_chain_idx+idx_p}'), file=f_out)

        cg = paf_retrieve_cigar(p)
        cg_ops = re_cigar.findall(cg)
        cg_lens = re_cigar.split(cg)
        match = None
        for i_cg, op in enumerate(cg_ops):
            if op == 'M':
                match = cg_lens[i_cg]
            elif op == 'D':
                assert match is not None
                print(f'{match} 0 {cg_lens[i_cg]}', file=f_out)
                match = None
            elif op == 'I':
                assert match is not None
                print(f'{match} {cg_lens[i_cg]} 0', file=f_out)
                match = None
            else:
                print('error: unexpected cigar op', op, file=sys.stderr)
                exit(1)
        print(f'{match}\n', file=f_out)


if __name__ == '__main__':
    args = parse_args()
    rdna_paf_to_chain(args)
