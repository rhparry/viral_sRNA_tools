# signatures_updated.py
# Rhys Parry - Updated June 2025
# This script computes overlap signatures from BAM files using pair and probability
# approaches, writing both raw and z-score values. In iterative mode, it outputs
# four matrices:
#   - Number of pairs (raw)
#   - Number of pairs (z-score)
#   - Overlap probability (raw)
#   - Overlap probability (z-score)

import argparse
from collections import defaultdict
import numpy
import pysam

def Parser():
    the_parser = argparse.ArgumentParser()
    the_parser.add_argument('--input', type=str, help="bam alignment file")
    the_parser.add_argument('--minquery', type=int)
    the_parser.add_argument('--maxquery', type=int)
    the_parser.add_argument('--mintarget', type=int)
    the_parser.add_argument('--maxtarget', type=int)
    the_parser.add_argument('--minscope', type=int)
    the_parser.add_argument('--maxscope', type=int)
    the_parser.add_argument('--output_h', type=str, help="h-signature file")
    the_parser.add_argument('--output_z', type=str, help="z-signature file")
    the_parser.add_argument('--mode', type=str, default="single", choices=["single", "iterative"])
    the_parser.add_argument('--minsize', type=int, help="Minimum read size for iteration")
    the_parser.add_argument('--maxsize', type=int, help="Maximum read size for iteration")
    the_parser.add_argument('--iter_output', type=str, help="Output prefix for iterative mode")
    return the_parser.parse_args()

class Map:
    def __init__(self, bam_file, minquery=23, maxquery=29, mintarget=23,
                 maxtarget=29, minscope=1, maxscope=19, output_h='',
                 output_z=''):
        self.bam_object = pysam.AlignmentFile(bam_file, 'rb')
        self.query_range = range(minquery, maxquery + 1)
        self.target_range = range(mintarget, maxtarget + 1)
        self.scope = range(minscope, maxscope + 1)
        self.H = open(output_h, 'w') if output_h != '/dev/null' else None
        self.Z = open(output_z, 'w') if output_z != '/dev/null' else None
        self.chromosomes = dict(zip(self.bam_object.references, self.bam_object.lengths))
        self.map_dict = self.create_map(self.bam_object)
        self.query_positions = self.compute_query_positions()
        if self.Z: self.Z.write(self.compute_signature_pairs())
        if self.H: self.H.write(self.compute_signature_h())
        if self.H: self.H.close()
        if self.Z: self.Z.close()

    def create_map(self, bam_object):
        d = defaultdict(list)
        for chrom in self.chromosomes:
            d[(chrom, 1, 'F')] = []
            d[(chrom, self.chromosomes[chrom], 'F')] = []
        for chrom in self.chromosomes:
            for read in bam_object.fetch(chrom):
                if read.is_reverse:
                    d[(chrom, read.reference_end, 'R')].append(read.query_alignment_length)
                else:
                    d[(chrom, read.reference_start+1, 'F')].append(read.query_alignment_length)
        return d

    def compute_query_positions(self):
        pos = defaultdict(list)
        for key in list(self.map_dict.keys()):
            chrom, coord, pol = key
            for i in self.scope:
                if pol == 'F' and (chrom, coord+i-1, 'R') in self.map_dict and len(self.map_dict[(chrom, coord+i-1, 'R')]) > 0:
                    pos[chrom].append(coord)
                    break
        for chrom in pos:
            pos[chrom] = sorted(set(pos[chrom]))
        return pos

    def countpairs(self, uppers, lowers):
        u = [s for s in uppers if s in self.query_range or s in self.target_range]
        l = [s for s in lowers if s in self.query_range or s in self.target_range]
        paired = 0
        used = set()
        for i, up in enumerate(u):
            for j, down in enumerate(l):
                if j in used:
                    continue
                if (up in self.query_range and down in self.target_range) or \
                   (up in self.target_range and down in self.query_range):
                    paired += 1
                    used.add(j)
                    break
        return paired

    def compute_signature_pairs(self):
        ft = defaultdict(dict)
        for chrom in self.chromosomes:
            for o in self.scope:
                ft[chrom][o] = 0
        for chrom in self.query_positions:
            for coord in self.query_positions[chrom]:
                for o in self.scope:
                    uppers = self.map_dict[(chrom, coord, 'F')]
                    lowers = self.map_dict[(chrom, coord+o-1, 'R')]
                    ft[chrom][o] += self.countpairs(uppers, lowers)
        for o in self.scope:
            ft['all_chromosomes'][o] = sum(ft[chrom][o] for chrom in self.chromosomes)
        return self.stringify_table(ft)

    def signature_tables(self):
        qr, tr = self.query_range, self.target_range
        Q, T = defaultdict(dict), defaultdict(dict)
        for (chrom, coord, pol), sizes in self.map_dict.items():
            for size in sizes:
                if size in qr or size in tr:
                    c = coord if pol == 'F' else -coord
                    if size in qr:
                        Q[chrom][c] = Q[chrom].get(c, 0) + 1
                    if size in tr:
                        T[chrom][c] = T[chrom].get(c, 0) + 1
        return Q, T

    def compute_signature_h(self):
        Q, T = self.signature_tables()
        ft = defaultdict(dict)
        for chrom in self.chromosomes:
            for o in self.scope:
                ft[chrom][o] = 0
        for chrom in Q:
            tqn = sum(Q[chrom].values())
            for coord in Q[chrom]:
                local = {o: 0 for o in self.scope}
                nt = 0
                for o in self.scope:
                    local[o] += Q[chrom][coord] * T[chrom].get(-coord - o + 1, 0)
                    nt += T[chrom].get(-coord - o + 1, 0)
                for o in self.scope:
                    try:
                        ft[chrom][o] += local[o] / nt / float(tqn)
                    except ZeroDivisionError:
                        continue
        gt = {o: 0 for o in self.scope}
        total_reads = sum(self.bam_object.count(chrom) for chrom in ft)
        for chrom in ft:
            for o in ft[chrom]:
                try:
                    gt[o] += ft[chrom][o] / total_reads * self.bam_object.count(chrom)
                except ZeroDivisionError:
                    continue
        ft['all_chromosomes'] = gt
        return self.stringify_table(ft)

    def stringify_table(self, ft):
        out = []
        for chrom in sorted(ft):
            vals = [ft[chrom][o] for o in sorted(ft[chrom])]
            m, s = numpy.mean(vals), numpy.std(vals)
            for o in sorted(ft[chrom]):
                z = 0 if s == 0 else (ft[chrom][o] - m)/s
                out.append(f"{chrom}\t{o}\t{ft[chrom][o]}\t{z}\n")
        return ''.join(out)

def run_iterative_analysis_full(bam, minsize, maxsize, minscope, maxscope, outprefix):
    sizes = list(range(minsize, maxsize + 1))
    labels = [str(i) for i in sizes] + [f"{minsize}-{maxsize}"]
    overlaps = list(range(minscope, maxscope + 1))
    pair_mat, pair_zmat, prob_mat, prob_zmat = [], [], [], []

    for size in sizes + [None]:
        minq = minsize if size is None else size
        maxq = maxsize if size is None else size
        m = Map(bam, minq, maxq, minq, maxq, minscope, maxscope, '/dev/null', '/dev/null')
        raw_pairs, z_pairs = {}, {}
        for line in m.compute_signature_pairs().splitlines():
            chrom, o, val, z = line.strip().split('\t')
            if chrom == 'all_chromosomes':
                raw_pairs[int(o)] = float(val)
                z_pairs[int(o)] = float(z)
        raw_probs, z_probs = {}, {}
        for line in m.compute_signature_h().splitlines():
            chrom, o, val, z = line.strip().split('\t')
            if chrom == 'all_chromosomes':
                raw_probs[int(o)] = float(val)
                z_probs[int(o)] = float(z)
        pair_mat.append([raw_pairs.get(k, 0.0) for k in overlaps])
        pair_zmat.append([z_pairs.get(k, 0.0) for k in overlaps])
        prob_mat.append([raw_probs.get(k, 0.0) for k in overlaps])
        prob_zmat.append([z_probs.get(k, 0.0) for k in overlaps])

    def write_matrix(fname, matrix):
        with open(fname, 'w') as f:
            f.write("ReadSize\t" + "\t".join(map(str, overlaps)) + "\n")
            for label, row in zip(labels, matrix):
                f.write(label + "\t" + "\t".join(f"{v:.6f}" for v in row) + "\n")

    write_matrix(f"{outprefix}_pairs.tsv", pair_mat)
    write_matrix(f"{outprefix}_pairs_z.tsv", pair_zmat)
    write_matrix(f"{outprefix}_prob.tsv", prob_mat)
    write_matrix(f"{outprefix}_prob_z.tsv", prob_zmat)

if __name__ == "__main__":
    args = Parser()
    if args.mode == "single":
        Map(args.input, args.minquery, args.maxquery, args.mintarget,
            args.maxtarget, args.minscope, args.maxscope, args.output_h,
            args.output_z)
    elif args.mode == "iterative":
        assert args.minsize and args.maxsize and args.iter_output, "--minsize, --maxsize, and --iter_output required"
        run_iterative_analysis_full(args.input, args.minsize, args.maxsize,
                                    args.minscope, args.maxscope,
                                    args.iter_output)
