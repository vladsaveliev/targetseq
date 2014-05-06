#!/usr/bin/env python
import numpy

import sys
import os
from os.path import join, splitext, basename, realpath
import subprocess
from datetime import datetime
#downlad hg19.genome
#https://github.com/arq5x/bedtools/tree/master/genomes

#TODO
# check on the input file format
# format result and calculation to the 2 decimal places on the header report    .00
# multi - sample report                                                         header only
#       sample1 sample2
#number 2       3
#bases  10      20

# check if samtools and bedtools exist
# log file
# yaml
# take folder name as a sample name ( first column on the report)
# give user an option to select type of the report to run ????

_header = ['Sample', 'Chr', 'Start', 'End', 'Exome Num', 'Direction', 'Length', 'MeanDepth', '1', '5', '10', '25', '50',
           '100', '500', '1000', '5000', '10000', '50000']
_groups = [[1, 0], [5, 0], [10, 0], [25, 0], [50, 0], [100, 0], [500, 0], [1000, 0], [5000, 0], [10000, 0], [50000, 0]]


def _call(cmdline, stdout=None):
    print ''
    print '*' * 70
    print datetime.now()
    print cmdline
    subprocess.call(cmdline.split(), stdout=stdout)


def _call_stdout(cmdline, stdout=subprocess.PIPE):
    print ''
    print '*' * 70
    print datetime.now()
    print cmdline
    return subprocess.Popen(cmdline.split(), stdout=stdout)


# TODO to check if input files a re not empty
def _call_and_write(cmdline, fpath, new_ext):
    base_name, ext = os.path.splitext(fpath)
    output_fpath = base_name + '.' + new_ext
    if not os.path.isfile(output_fpath):
        _call(cmdline, open(output_fpath, 'w'))
        #TODO check if we have file
    return output_fpath


def gnu_sort(bed_path):
    cmdline = 'sort -k1,1V  -k2,2n -k3,3n %s' % (bed_path)
    return _call_and_write(cmdline, bed_path, 'sorted.bed')


def samtool_depth_range(bam_path, range):
    cmdline = 'samtools depth %s -r %s' % (bam_path, range)
    return _call_stdout(cmdline, stdout=subprocess.PIPE)


def intersect_bed(bed_path_gene, bed_path_capture):
    cmdline = 'intersectBed -a %s -b %s -u' % (  bed_path_gene, bed_path_capture)
    return _call_and_write(cmdline, bed_path_gene, 'intersect.bed')


def number_of_mapped_reads(bam):
    cmdline = 'samtools view -c -F 4 %s' % ( bam)
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    return proc.stdout.readline()


def number_of_unmapped_reads(bam):
    cmdline = 'samtools view -c -f 4 %s' % ( bam)
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    return proc.stdout.readline()


def number_of_reads(bam):
    cmdline = 'samtools view -c %s' % bam
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    return proc.stdout.readline()


def reads_on_target(bed, bam):
    cmdline = 'samtools view -c -F 4 %s -L %s' % (  bam, bed)
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    return proc.stdout.readline()

# TODO very slow
def total_aligned_base_reads(bam):
    cmdline = 'samtools depth %s ' % (bam)
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    count = 0
    while True:
        coverage_line = proc.stdout.readline()
        if coverage_line != '':
            values = coverage_line.strip().split("\t")
            count += int(values[2])
        else:
            break
    return count


def get_analytics_target_depth(bed, bam):
    cmdline = 'samtools depth %s -b %s ' % ( bam, bed)
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    covered_bases = 0
    total_reads_on_target = 0
    max_depth = 0
    while True:
        coverage_line = proc.stdout.readline()
        if coverage_line != '':
            values = coverage_line.strip().split("\t")
            depth_value = int(values[2])
            total_reads_on_target += depth_value
            if max_depth < depth_value: max_depth = depth_value
            covered_bases += 1

            for i in _groups:
                if depth_value >= i[0]:
                    i[1] += 1

        else:
            break

    return _groups, covered_bases, total_reads_on_target, max_depth

# TODO how to pass the data stream to samtools vs. creating file
def get_pedded_bed_file(bed, genome, padding):
    cmdline = 'bedtools slop -i %s -g %s -b %s' % (bed, genome, padding)
    return _call_and_write(cmdline, bed, "padded.bed")


def count_by_group(depths):
    for depth_value in depths:
        for i in _groups:
            if depth_value >= i[0]:
                i[1] += 1
    return [100 * k[1] / len(depths) for k in _groups]


def get_depth_by_bed_range(bam, range):
    proc = samtool_depth_range(bam, range)
    depths = []
    while True:
        coverage_line = proc.stdout.readline()
        values = coverage_line.strip().split("\t")
        if coverage_line != '' and len(values) > 2:
            depths.append(values[2])
        else:
            break

    return depths

###REPORTS
def print_line(sample_name, line, length, mean, start_end_count):
    return "{0} \t {1} \t {2} \t {3} \t ".format(sample_name, line, length, mean) + "\t".join(start_end_count) + '\n'


def print_name_value(name, value):
    return "{0}: {1}\n".format(name, value)


def write_report_line(bam, base_name, out, bed_line):
    values = bed_line.strip().split("\t")
    chr = values[0]
    start = values[1]
    end = values[2]
    range = "{0}:{1}-{2}".format(chr, start, end)
    depths = get_depth_by_bed_range(bam, range)
    length = (int(end) - int(start)) + 1
    if len(depths) > 0:
        depths_int = map(int, depths)
        depth_mean = "{0:.2f}".format(numpy.mean(depths_int))
        counts = count_by_group(depths_int)
        out.write(print_line(base_name, bed_line.strip(), length, depth_mean, map(str, counts)))
    else:
        out.write(print_line(base_name, bed_line.strip(), length, 'na', []))

# TODO check on division by 0
def run_header_report(bed, bam, genome, padding=250):
    print ''
    print '*' * 70
    print datetime.now()
    print "Printing header report"
    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'header.report'
    print output_path
    with open(output_path, 'w') as out:
        vnumber_of_mapped_reads = int(number_of_mapped_reads(bam))
        out.write(print_name_value("Number of mapped reads", vnumber_of_mapped_reads))

        vnumber_of_unmapped_reads = int(number_of_unmapped_reads(bam))
        out.write(print_name_value("Number of unmapped reads", vnumber_of_unmapped_reads))

        vnumber_of_reads = int(number_of_reads(bam))
        out.write(print_name_value("Number of reads", vnumber_of_reads))
        out.write(print_name_value("Percent mapped reads", 100 * vnumber_of_mapped_reads / vnumber_of_reads))

        vnumber_reads_on_target = int(reads_on_target(bed, bam))
        out.write(print_name_value("Number of reads on target", vnumber_reads_on_target))

        percent_reads_on_target = 100 * vnumber_reads_on_target / vnumber_of_mapped_reads
        out.write(print_name_value("Percent reads on target", percent_reads_on_target))

        vnumber_padded_reads_on_target = int(reads_on_target(get_pedded_bed_file(bed, genome, padding), bam))
        out.write(print_name_value("Percent reads on padded target",
                                   100 * vnumber_padded_reads_on_target / vnumber_of_mapped_reads))

        vtotal_aligned_base_reads = total_aligned_base_reads(bam)
        out.write(print_name_value("Total aligned base reads", vtotal_aligned_base_reads))

        target_groups, cov, vtotal_aligned_base_reads_on_target, max_depth = get_analytics_target_depth(bed, bam)
        out.write(print_name_value("Total base reads on target", vtotal_aligned_base_reads_on_target))

        percent_based_reads_on_target = 100 * vtotal_aligned_base_reads_on_target / vtotal_aligned_base_reads
        out.write(print_name_value("Percent base reads on target", percent_based_reads_on_target))
        out.write(print_name_value("Bases in targeted reference ", cov))
        out.write(print_name_value("Bases covered (at least 1x) ", target_groups[0][1]))
        out.write(print_name_value("Average base coverage depth ", vtotal_aligned_base_reads_on_target / cov))
        out.write(print_name_value("Maximum base read depth", max_depth))

        #Std.Dev base read depth
        [out.write(print_name_value("Target coverage at " + str(group[0]) + "x", 100 * group[1] / cov)) for group in
         target_groups]

    print datetime.now()
    print "Done with header report"


def run_cov_report(bed, bam):
    print ''
    print '*' * 70
    print datetime.now()

    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'report'

    print 'Creating report path: ' + output_path

    bed_sorted_path = gnu_sort(bed)

    with open(bed_sorted_path) as sorted_bed, open(output_path, 'w') as out:
        out.write('\t'.join(_header) + '\n')
        for bed_line in sorted_bed.readlines():
            write_report_line(bam, base_name, out, bed_line)


def run_cov_gene_report(gene_bed, capture_bed, bam):
    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'gene.report'

    # a_bed_sorted = gnu_sort(a_bed)
    # b_bed_sorted = gnu_sort(b_bed)
    int_bed = intersect_bed(gene_bed, capture_bed)

    with open(int_bed) as bed, open(output_path, 'w') as out:
        out.write("\t".join(_header) + '\n')
        for bed_line in bed.readlines():
            write_report_line(bam, base_name, out, bed_line)


def run_cov_exomes_report(gene_bed, exome_bed, capture_bed, bam):
    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'exomes.report'

    int_gene_bed = intersect_bed(gene_bed, capture_bed)
    int_gene_exone_bed = intersect_bed(exome_bed, int_gene_bed)

    with open(int_gene_exone_bed) as  bed, open(output_path, 'w') as out:
        out.write("\t".join(_header) + '\n')
        for bed_line in bed.readlines():
            write_report_line(bam, base_name, out, bed_line)


# input  aligned.bed, gene_exome.bed, gene.bed, bam, genome
if __name__ == '__main__':
    bed_path_gene_exome = "/home/alla/worksplace/TargetSeq/data/pim/UCSC_exons_modif_canonil.sorted.bed"
    bed_path_gene = "/home/alla/worksplace/TargetSeq/data/pim/UCSC_canonical.bed"

    #bam_path ="/home/alla/worksplace/TargetSeq/data/D10_S7_L001.sorted.realign.bam"
    #bed_path_capture = "/home/alla/worksplace/TargetSeq/data/AZ_internal_13_DS_merged.bed"
    bam_path = "/home/alla/worksplace/TargetSeq/data/pim/EOL-1-1-ready.bam"
    bed_path_capture = "/home/alla/worksplace/TargetSeq/data/pim/AgilentSureSelectV4_51M.bed"
    genome_path = "/home/alla/worksplace/TargetSeq/data/bedtools/genomes/human.hg19.genome"

    run_header_report(bed_path_capture, bam_path, genome_path, padding=250)

    # run_cov_exomes_report(  bed_path_gene, bed_path_gene_exome,bed_path_capture,bam_path)


