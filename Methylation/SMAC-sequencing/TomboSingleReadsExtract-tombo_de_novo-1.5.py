##################################
#                                #
# Last modified 2018/12/08       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import h5py
import numpy as np    
import os

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def run():

    if len(sys.argv) < 3:
        print('usage: python %s tombo.per_read_stats genome.fa outfile_prefix [-m5C-only] [-m6A-only] [-CG-only] [-CG-CG-only] [-GC-only] [-m6A-CG-only] [-m6A-GC-only] [-m6A-GC-CG-only] [-doT] [-T-only] [-excludeChr chr1[,chr2,...,chrN]]' % sys.argv[0])
        print('\tnote: by default, As in all contexts and Cs in all contexts will be printed out')
        print('\tnote: the [-m6A-CG-only option] will print out As in all contexts and Cs in CpG context')
        print('\tnote: the [-m6A-GC-only option] will print out As in all contexts and Cs in GpC context')
        print('\tnote: the [-m6A-CG-GC-only option] will print out As in all contexts and Cs in GpC or CpG context')
        print('\tnote: the [-m5C-only] option will print out Cs in all contexts')
        print('\tnote: the [-CG-only] and [-GC-only] [-CG-GC-only options only apply if they [-m5C-only] has been specified')
        print('\tnote: the [-doT] option is off by defaultt; the [-T-only] option only applies if they [-doT] option has been enabled')
        sys.exit(1)

    do5C = True
    do6A = True
    CGonly = False
    GConly = False
    CGGConly = False

    ExcludedChrs = {}
    if '-excludeChr' in sys.argv:
        for chr in sys.argv[sys.argv.index('-excludeChr') + 1].split(','):
            ExcludedChrs[chr] = 1

    if '-m5C-only' in sys.argv:
        do6A = False
        GConly = False
        CGonly = False
        print('will only output m5C positions')
        if '-CG-only' in sys.argv:
            CGonly = True
            print('will only output m5C positions in CpG context')
        if '-GC-only' in sys.argv:
            GConly = True
            print('will only output m5C positions in GpC context')
        if GConly and CGonly:
            print('incompatible options, exiting')
            sys.exit(1)
        if '-CG-GC-only' in sys.argv:
            CGGConly = True
            print('will only output m5C positions in GpC or CpG context')

    dom6AGCCGonly = False
    if '-m6A-GC-CG-only' in sys.argv:
        dom6AGCCGonly = True
        print('will only output m6A, GpC and CpG positions')

    dom6AGConly = False
    if '-m6A-GC-only' in sys.argv:
        dom6AGConly = True
        print('will only output m6A and GpC positions')

    dom6ACGonly = False
    if '-m6A-CG-only' in sys.argv:
        dom6ACGonly = True
        print('will only output m6A and CpG positions')

    if '-m6A-only' in sys.argv:
        do5C = False
        print('will only output m6A positions')

    doT = False
    if '-doT' in sys.argv:
        doT = True
        print('will output T positions too')
        if '-T-only' in sys.argv:
            do5C = False
            do6A = False
            print('will only output T positions')

    tombo = sys.argv[1]
    fasta = sys.argv[2]
    outprefix = sys.argv[3]

    GenomeDict={}
    sequence=''
    inputdatafile = open(fasta)
    for line in inputdatafile:
        if line[0]=='>':
            if sequence != '':
                GenomeDict[chr] = ''.join(sequence).upper()
            chr = line.strip().split('>')[1]
            sequence=[]
            Keep=False
            continue
        else:
            sequence.append(line.strip())
    GenomeDict[chr] = ''.join(sequence).upper()

    for chr in list(ExcludedChrs.keys()):
        if chr in GenomeDict:
           del GenomeDict[chr]

    print('finished inputting genomic sequence')

    ReadDict = {}

    fMfile = h5py.File(tombo, 'r')

    groupC = fMfile['Statistic_Blocks']
    for block in groupC.values():
        reads =  block['read_ids']
        RIDtoReadNameDict = {}
        RIDnumber = 0
        for rid in reads:
            RIDtoReadNameDict[RIDnumber] = rid
            RIDnumber += 1
        for BB in block.attrs.items():
            if BB[0] == 'chrm':
                chr = BB[1].split('|')[0]
            if BB[0] == 'strand':
                strand = BB[1]
            if BB[0] == 'start':
                start = BB[1]
        if chr in GenomeDict:
            pass
        else:
            print('chromosome not found in genome.fa file, skipping', chr)
            continue
        for (pos,ll,ridNum) in block['block_stats']:
            rid = RIDtoReadNameDict[ridNum]
            if (chr,strand,rid) in ReadDict:
                pass
            else:
                ReadDict[(chr,strand,rid)] = {}
#            seq = GenomeDict[chr][pos-1:pos+2]
            seq1 = GenomeDict[chr][pos:pos+2]
            seq2 = GenomeDict[chr][pos-1:pos+1]
            if strand == '+':
                if GConly:
                    if seq2 == 'GC':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif CGonly:
                    if seq1 == 'CG':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif CGGConly:
                    if seq1 == 'CG' or seq2 == 'GC':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6AGCCGonly:
                    if seq1 == 'CG' or seq2 == 'GC':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'A':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6AGConly:
                    if seq2 == 'GC':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'A':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6ACGonly:
                    if seq1 == 'CG':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'A':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                else:
                    if do5C and (seq1[0] == 'C'):
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif do6A and seq1[0] == 'A':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '6A+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif doT and seq1[0] == 'T':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, 'T+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
            if strand == '-':
                if GConly:
                    if seq1 == 'GC':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif CGonly:
                    if seq2 == 'CG':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif CGGConly:
                    if seq2 == 'CG' or seq1 == 'GC':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6AGCCGonly:
                    if seq2 == 'CG' or seq1 == 'GC':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'T':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '6A-, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6AGConly:
                    if seq1 == 'GC':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'T':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '6A-, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                elif dom6ACGonly:
                    if seq2 == 'CG':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C+, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif seq1[0] == 'T':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '6A-, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                else:
#                    print seq1, rid, pos, ll, chr, GenomeDict[chr][pos-5:pos+5]
                    if do5C and (seq1[0] == 'G'):
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '5C-, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif do6A and seq1[0] == 'T':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, '6A-, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll
                    elif doT and seq1[0] == 'A':
                        if pos in ReadDict[(chr,strand,rid)]:
                            print('positions already seen', (chr,strand,rid), pos, 'T-, exiting', sys.exit(1))
                        ReadDict[(chr,strand,rid)][pos] = ll

    print('finished processing file')

    reads = list(ReadDict.keys())
    reads.sort()

    outfile = open(outprefix + '.reads.tsv', 'w')

    for (chr,strand,rid) in reads:
        positions = list(ReadDict[(chr,strand,rid)].keys())
        positions.sort()
        if len(positions) == 0:
            print('skipping:', (chr,strand,rid), ReadDict[(chr,strand,rid)])
            continue
        outline = chr + '\t' + str(positions[0]) + '\t' + str(positions[-1]) + '\t' + strand + '\t' + rid + '\t' + '.'
        Ps = ''
        LLs = ''
        for pos in positions:
            Ps = Ps + str(pos) + ','
            LLs = LLs + "{0:.2f}".format(ReadDict[(chr,strand,rid)][pos]) + ','
        outline = outline + '\t' + Ps[0:-1]
        outline = outline + '\t' + LLs[0:-1]
        outfile.write(outline + '\n')

    outfile.close()

    
run()
