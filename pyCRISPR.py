#!/usr/bin/env python3
import csv
import collections
import itertools
import multiprocessing
import re
import subprocess
import functools
from quicksect import IntervalNode

Row = collections.namedtuple('Row', ['chrom', 'start', 'end'])
Guide = collections.namedtuple('Guide', ['chrom', 'start', 'end', 'sequence', 'strand'])

def split_target(target):
    # Split BED records longer than 3kb into 2kb tiles.
    chrom = target.chrom
    start = int(target.start)
    end = int(target.end)
    size = 3000
    window = 2000
    length = end - start
    if length <= size:
        yield target
    else:
        sub_regions = []
        for i in range(0, length, size):
            sub_chrom = chrom
            sub_start = start + i
            sub_end = start + i + window
            if sub_end > end:
                sub_end = end
                yield Row(chrom=sub_chrom, start=sub_start, end=sub_end)
            else:
                yield Row(chrom=sub_chrom, start=sub_start, end=sub_end)
            
def get_revcomp(dna_seq):
    # Return reverse complement of DNA sequence.
    return dna_seq[::-1].translate(str.maketrans('AaCcGgTtNn', 'TtGgCcAaNn'))

def get_dna_seq(chrom, start, end, revcomp):
    # Get DNA sequence.
    genomic_region = chrom + ':' + str(start) + '-' + str(end)
    samtools_command = ['samtools', 'faidx', 'mm10.fa', genomic_region]
    output_bytes = subprocess.check_output(samtools_command)
    output_text = output_bytes.decode('utf-8')
    fasta_header = output_text.split()[0]
    fasta_sequence = ''.join(output_text.split()[1:])
    if revcomp:
        return get_revcomp(fasta_sequence)
    return fasta_sequence

def create_oligo(guide_seq1, guide_seq2):
    # Create oligo structure for library.
    # Return DNA sequence of oligo structure.
    guide_seq1_no_pam = guide_seq1[:-3]
    guide_seq2_no_pam = guide_seq2[:-3]
    oligo_seq = 'GCAGATGGCTCTTTGTCCtaGAAGACAAcaccg{0}gtttGAAGAGCCATGACTGCTCTTCcctcg{1}GTTTTAGTCTTCTCGTCGC'.format(guide_seq1_no_pam, guide_seq2_no_pam)
    return oligo_seq

def find_guides(dna_seq):
    # Find gRNA site in DNA sequence.
    # Return match object if found.
    guide_regex = re.compile('(?=([ACGT]{20}GG))')
    matches = re.finditer(guide_regex, dna_seq)
    return matches

def find_poly_nucl(dna_seq):
    # Find poly A,C,G or T in DNA sequence.
    # Return TRUE if found, or FALSE if not.
    poly_regex = re.compile('(A{5,}|C{5,}|G{5,}|T{4,})')
    if re.search(poly_regex, dna_seq):
        return True
    return False

def find_repeat_kmers(dna_seq):
    pass

def find_gene_overlap(chrom, start, end):
    # Return TRUE if gRNA coordinates overlap with a gene.
    genomic_region = chrom + ':' + str(start) + '-' + str(end)
    tabix_command = ['tabix', 'genes.mm10.sorted.bed.gz', genomic_region]
    output_bytes = subprocess.check_output(tabix_command)
    output_text = output_bytes.decode('utf-8')
    if output_text:
        return True
    else:
        return False

def find_multiple_sapi_site(dna_seq):
    # Find more than two SapI site in DNA sequence.
    # Return TRUE if found, or FALSE if not.
    sap1_regex = re.compile('GAAGAGC|CTTCTCG')
    matches = re.findall(sap1_regex, dna_seq)
    if len(matches) > 2:
        return True
    return False

def find_sapi_site(dna_seq):
    # Find SapI site in DNA sequence.
    # Return TRUE if found, or FALSE if not.
    sap1_regex = re.compile('GAAGAGC|CTTCTCG')
    if re.search(sap1_regex, dna_seq):
        return True
    return False

def find_multiple_bbsi_site(dna_seq):
    # Find more than two BbsI site in DNA sequence.
    # Return TRUE if found, or FALSE if not.
    bbs1_regex = re.compile('GAAGAC|CTTCTG')
    matches = re.findall(bbs1_regex, dna_seq)
    if len(matches) > 2:
        return True
    return False

def find_bbsi_site(dna_seq):
    # Find BbsI site in DNA sequence.
    # Return TRUE if found, or FALSE if not.
    bbs1_regex = re.compile('GAAGAC|CTTCTG')
    if re.search(bbs1_regex, dna_seq):
        return True
    return False

def parse_bed_file(bed_file):
    # Parse BED file.
    # Return tuple of BED record.
    with open(bed_file) as f:
        f_bed = csv.reader(f, delimiter='\t')
        return [Row(*r) for r in f_bed]

def find_offtarget_site(guide_start, guide_sequence):
    'Find off-targets'
    num_offtargets = 0
    bowtie_command = ('bowtie -n 3 -l 10 -a -S --quiet mm10 -c {0} | samtools calmd -e - mm10.fa | samtools view - '.format(guide_sequence))
    output_bytes = subprocess.check_output(bowtie_command, shell=True)
    output_text = output_bytes.decode('utf-8')
    sam_alignments = output_text.rstrip().split('\n')
    for sam_alignment in sam_alignments:
        sam_fields = sam_alignment.split('\t')
        sam_start = sam_fields[3]
        if guide_start == int(sam_start):
            pass
        else:
            sam_alignment = sam_fields[9]
            perfect_region_mismatches = sam_alignment[10:].replace('=','')
            variable_region_mismatches = sam_alignment[0:10].replace('=','')
            if perfect_region_mismatches == 0 and variable_region_mismatches <= 3:
                num_offtargets += 1

    if num_offtargets > 0:
        return True
    else:
        return False                

def get_plus_guides(flank_chrom, flank_start, flank_end, flank_seq):
    # Collect valid guides.
    valid_plus_guides = []
    # Search on plus strand.
    for match in find_guides(flank_seq):
        # Create gRNA tuple.
        guide_chrom = flank_chrom
        guide_start = flank_start + match.start() - 1
        guide_end = flank_start + match.end() - 1 + 22
        guide_sequence = match.group(1)
        guide_strand = '+'
        guide_tuple = Guide(guide_chrom, guide_start, guide_end, guide_sequence, guide_strand)
        # Regex searches.
        if find_poly_nucl(guide_sequence):
            pass
        elif find_bbsi_site(guide_sequence):
            pass
        elif find_sapi_site(guide_sequence):
            pass
        elif find_gene_overlap(guide_chrom, guide_start, guide_end):
            pass
        elif find_offtarget_site(guide_start, guide_sequence):
            pass
        else:
            valid_plus_guides.append(guide_tuple)
    return valid_plus_guides

def get_minus_guides(flank_chrom, flank_start, flank_end, flank_seq):
    # Collect valid guides.
    valid_minus_guides = []
    # Search on minus strand.
    for match in find_guides(flank_seq):
        # Create gRNA tuple.
        flank_size = flank_end - flank_start
        guide_chrom = flank_chrom
        guide_start = flank_start + (flank_size - match.end()) - 22
        guide_end = flank_end - match.start()
        guide_sequence = match.group(1)
        guide_strand = '-'
        guide_tuple = Guide(guide_chrom, guide_start, guide_end, guide_sequence, guide_strand)
        # Regex searches.
        if find_poly_nucl(guide_sequence):
            pass
        elif find_bbsi_site(guide_sequence):
            pass
        elif find_sapi_site(guide_sequence):
            pass
        elif find_gene_overlap(guide_chrom, guide_start, guide_end):
            pass
        elif find_offtarget_site(guide_start, guide_sequence):
            pass
        else:
            valid_minus_guides.append(guide_tuple)
    return valid_minus_guides

def build_interval_tree(intervals):
    root = IntervalNode(intervals[0].start, intervals[0].end, other=intervals[0])
    return functools.reduce(lambda tree, x: tree.insert(x.start, x.end, other=x), intervals[1:], root)

class Interval(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.removed = False

def remove_overlapping_guides(sorted_guides, flank):
    
    # Define intervals
    intervals = [Interval(guide.start, guide.end) for guide in sorted_guides]
    
    # Sort by start coordinate
    if flank == 'upstream':
        intervals.sort(key=lambda x: (x.end, x.end), reverse=True)
    
    # Sort by end coordinate
    elif flank == 'downstream':
        intervals.sort(key=lambda x: (x.end, x.start), reverse=False)
    tree = build_interval_tree(intervals)
    results = []
    for smallest in intervals:
        if smallest.removed:
            continue
        smallest.removed = True
        results.append([smallest.start, smallest.end])
        
        tree.intersect(smallest.start, smallest.end, lambda x: setattr(x.other, 'removed', True))

    valid_guides = []
    for result in results:
        for guide in sorted_guides:
            if result[0] == guide.start and result[1] == guide.end:
                valid_guides.append(guide)

    return valid_guides

def main(row):
    
    # flanking region
    flank_size = 500
    
    # enhancer region
    enhancer_chrom = row.chrom
    enhancer_start = int(row.start)
    enhancer_end = int(row.end)

    #print("Enhancer region:", enhancer_chrom, enhancer_start, enhancer_end)

    # offset size
    offset_size = round((10/100)*(enhancer_end-enhancer_start))
    
    # upstream region
    upstream_chrom = enhancer_chrom
    upstream_start = enhancer_start - flank_size
    upstream_end = enhancer_start + offset_size
    upstream_plus_seq = get_dna_seq(upstream_chrom, upstream_start, upstream_end, revcomp=False)
    upstream_minus_seq = get_dna_seq(upstream_chrom, upstream_start, upstream_end, revcomp=True)
    
    #print("Upstream region:", upstream_chrom, upstream_start, upstream_end)

    # downstream region
    downstream_chrom = enhancer_chrom
    downstream_start = enhancer_end - offset_size
    downstream_end = enhancer_end + flank_size
    downstream_plus_seq = get_dna_seq(downstream_chrom, downstream_start, downstream_end, revcomp=False)
    downstream_minus_seq = get_dna_seq(downstream_chrom, downstream_start, downstream_end, revcomp=True)
    
    #print("Downstream region:", downstream_chrom, downstream_start, downstream_end)

    # upstream target and guides
    upstream_plus_guides = get_plus_guides(upstream_chrom, upstream_start, upstream_end, upstream_plus_seq)
    upstream_minus_guides = get_minus_guides(upstream_chrom, upstream_start, upstream_end, upstream_minus_seq)
    
    #print("Number of upstream plus strand guides:", len(upstream_plus_guides))
    #print("Number of upstream minus strand guides:", len(upstream_minus_guides))
    
    # downstream target and guides
    downstream_plus_guides = get_plus_guides(downstream_chrom, downstream_start, downstream_end, downstream_plus_seq)
    downstream_minus_guides = get_minus_guides(downstream_chrom, downstream_start, downstream_end, downstream_minus_seq)

    #print("Number of downstream plus strand guides:", len(downstream_plus_guides))
    #print("Number of downstream minus strand guides:", len(downstream_minus_guides))

    # combine plus and minus guides
    upstream_guides = upstream_plus_guides + upstream_minus_guides
    downstream_guides = downstream_plus_guides + downstream_minus_guides
    
    #print("Number of upstream guides:", len(upstream_guides))
    #print("Number of downstream guides:", len(downstream_guides))

    results = []
    
    if len(upstream_guides) != 0 and len(downstream_guides) != 0:
        
        #print("Removing overlapping guides...")
        
        # remove overlaps
        upstream_guides_sorted_overlap = remove_overlapping_guides(upstream_guides, flank='upstream')
        downstream_guides_sorted_overlap = remove_overlapping_guides(downstream_guides, flank='downstream')
        
        #print("Number of upstream guides with overlapping guides removed:", len(upstream_guides_sorted_overlap))
        #print("Number of downstream guides with overlapping guides removed:", len(downstream_guides_sorted_overlap))
        
        #print("Oligo library:")
        
        # Create combinations
        for up_guide, down_guide in zip(upstream_guides_sorted_overlap, downstream_guides_sorted_overlap):
            
            # Create library
            oligo_seq = create_oligo(up_guide.sequence, down_guide.sequence)
            
            # last minute check for enzyme sites
            if find_multiple_sapi_site(oligo_seq) or find_multiple_bbsi_site(oligo_seq):
                pass
            else:
                # write results
                target_region = '\t'.join([str(row.chrom), str(row.start), str(row.end)])
                up_region = '\t'.join([str(up_guide.chrom), str(up_guide.start), str(up_guide.end), str(up_guide.strand)])
                down_region = '\t'.join([str(down_guide.chrom), str(down_guide.start), str(down_guide.end), str(down_guide.strand)])
                result_row = '\t'.join([target_region, up_region, down_region, oligo_seq + '\n'])
                results.append(result_row)

    return results
    
def mp_handler():
    # Handle multi-processing output
    p = multiprocessing.Pool(20)
    targets = [target for target in parse_bed_file('Sample_ESC_Enhancers.sorted.bed')]
    sub_targets = []
    for target in targets:
        for sub_target in split_target(target):
            sub_targets.append(sub_target)
    with open('results.txt', 'w') as f:
        for result in p.imap(main, sub_targets):
            print(result)
            for guide in result:
                f.write(guide)

if __name__ == '__main__':
    mp_handler()
