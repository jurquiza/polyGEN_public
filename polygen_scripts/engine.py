# Import all necessary packages

import numpy as np
import csv
import itertools
import pandas as pd
import warnings
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


# Define infrastructure for following code
class Part:
    def __init__(self,name,type,sequence):
        self.name = name
        self.sequence = sequence
        self.type = type
    description = ""
    primer_forward = ""
    primer_reverse = ""
    bridge_with_next_part = ""
    bridge_with_previous_part = ""
    primer_forward_tm = 0
    primer_reverse_tm = 0
    
def builds(parts_list):
    import copy
    builds_list = []
    builds = 1
    for i in range(builds):
        build = []
        for part in parts_list:
            build.append(copy.copy(part))
        builds_list.append(build)
    return builds_list


# Optimize the overhangs used in Golden Gate assembly
def golden_gate_optimization(parts_list,gg_overhangs):

    golden_gate_overhangs=gg_overhangs
    seq_matches = []
    for x in range(len(parts_list)-1):
        seq_matches.append([])
        for overhang in golden_gate_overhangs:
            if parts_list[x].type == 'pegRNA':
                if overhang in parts_list[x].sequence[96:]:
                    seq_matches[x].append(overhang)
                else:
                    pass
            elif parts_list[x].type == 'gRNA':
                if parts_list[x+1].type == 'gRNA' or parts_list[x+1].type == 'pegRNA':
                    if overhang in parts_list[x+1].sequence[:20]:
                        seq_matches[x].append(overhang)
                    else:
                        pass
                elif parts_list[x+1].type == 'sRNA':
                    if overhang in parts_list[x+1].sequence[:-len(tRNA)]:
                        seq_matches[x].append(overhang)
                    else:
                        pass
            elif parts_list[x].type == 'sRNA':
                if parts_list[x+1].type == 'gRNA' or parts_list[x+1].type == 'pegRNA':
                    if overhang in parts_list[x+1].sequence[:20]:
                        seq_matches[x].append(overhang)
                    else:
                        pass
                elif parts_list[x+1].type == 'sRNA':
                    if overhang in parts_list[x+1].sequence[:-len(tRNA)]:
                        seq_matches[x].append(overhang)
                    else:
                        pass
            else: # if part is tRNA
                if parts_list[x+1].type == 'gRNA' or parts_list[x+1].type == 'pegRNA':
                    if overhang in parts_list[x+1].sequence[:20]:
                        seq_matches[x].append(overhang)
                    else:
                        pass
                elif parts_list[x+1].type == 'sRNA':
                    if overhang in parts_list[x+1].sequence[:-len(tRNA)]:
                        seq_matches[x].append(overhang)
                    else:
                        pass
    combs = []
    for x in itertools.product(*seq_matches):
        combs.append(x)
    for comb in combs:
        if len(comb) == len(set(comb)):
            return comb
    # if there are no possible combinations
    return None


# Define functions to easy the main computation
def reverse_complement(sequence):
    rev_comp = ""
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        rev_comp = Watson_Crick[base] + rev_comp
    return rev_comp

def complement(sequence):
    comp = ''
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        comp += Watson_Crick[base]
    return comp

def reverse(sequence):
    rev = ''
    for base in sequence:
        rev = base + rev
    return rev


# set template sequences tRNA and gRNA scaffold
tRNA = 'AACAAAGCACCAGTGGTCTAGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA'
scaffld = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'


# Perform scarless Golden Gate assembly computation with provided parts
def scarless_gg(parts_list, primer_tm_range, max_annealing_len, bb_overlaps, additional_overhangs):
    bb_overlaps = [i.lower() for i in bb_overlaps]
    additional_overhangs = [i.lower() for i in additional_overhangs]
    
    # format parts list
    builds_list = builds(parts_list)

    # Go through parts and write all known annotations into list
    new_builds_list = []
    ftrs = []
    for unpacked_list in builds_list:
        mmry = 13
        for part in unpacked_list:
            part.sequence = part.sequence.lower()
            ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+len(part.sequence), strand=1), type=part.name))
            if part.type == 'pegRNA':
                ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+20, strand=1), type='spacer'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+96, mmry+len(part.sequence)-13, strand=1), type='RT template'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+len(part.sequence)-13, mmry+len(part.sequence), strand=1), type='PBS'))
            elif part.type == 'gRNA':
                ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+20, strand=1), type='spacer'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+20+len(scaffld), mmry+len(part.sequence), strand=1), type='tRNA'))
            elif part.type == 'sRNA':
                ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+len(part.sequence)-len(tRNA), strand=1), type='sRNA'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+len(part.sequence)-len(tRNA), mmry+len(part.sequence), strand=1), type='tRNA'))
            mmry += len(part.sequence)

        # Iterate through overhang sets with increasing size until fitting one is found
        gg_opt = None
        breakit = False
        for p in range(10,51):
            for q in range(5):
                with open('overhangsets/setsof%s.csv'%p,'r') as f:
                    reader = csv.reader(f, delimiter=",")
                    temp = list(reader)[1:]
                    
                # Only grab sets that include all existing overhangs and delete the existing from the set
                if all(i in [x.lower() for x in temp[q]] for i in additional_overhangs+bb_overlaps):
                    free_overhangs = [i for i in [x.lower() for x in temp[q]] if i not in additional_overhangs+bb_overlaps]
                else:
                    continue
                
                gg_overhangs = []
                for x in free_overhangs:
                    gg_overhangs.append(x)
                gg_opt = golden_gate_optimization(unpacked_list,gg_overhangs)
                if gg_opt is not None:
                    breakit = True
                if breakit:
                    break
            if breakit:
                break
        
        # No sets include all existing overhangs
        if gg_opt is None:
            warnings.warn('The given combination of existing overhangs is not compatible with an optimal overhang set. '
                          'Computing the best overhang set not including the given existing overhangs. There might be '
                          'interference between overhangs.')
            for p in range(10,51):
                for q in range(5):
                    with open('overhangsets/setsof%s.csv'%p,'r') as f:
                        reader = csv.reader(f, delimiter=",")
                        temp = list(reader)[1:]

                    # Only grab sets that include no existing overhangs
                    if any(i in [x.lower() for x in temp[q]] for i in additional_overhangs+bb_overlaps):
                        continue
                    else:
                        free_overhangs = [x.lower() for x in temp[q]]

                    gg_overhangs = []
                    for x in free_overhangs:
                        gg_overhangs.append(x)
                    gg_opt = golden_gate_optimization(unpacked_list,gg_overhangs)
                    if gg_opt is not None:
                        breakit = True
                    if breakit:
                        break
                if breakit:
                    break
        
        # Raise error if no solution is possible
        if gg_opt is None:
            raise ValueError('The given combination of existing overhangs does not allow for an optimal overhang set.')
                
                
        #Modify sequences and design primers
        for i in range(len(unpacked_list)):

            # If current part is first part, define forward primer with left backbone overlap and reverse primer ordinarily
            if i == 0:
                unpacked_list[i].primer_forward = 'taggtctcc' + reverse_complement(bb_overlaps[0]) + unpacked_list[i].sequence[:max_annealing_len]
                unpacked_list[i].sequence = 'taggtctcc' + reverse_complement(bb_overlaps[0]) + unpacked_list[i].sequence
                if unpacked_list[i].type == 'pegRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76-max_annealing_len:unpacked_list[i].sequence.find(gg_opt[i],len(unpacked_list[i].sequence)-35)] + gg_opt[i] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i],len(unpacked_list[i].sequence)-35):] + unpacked_list[i+1].sequence[:max_annealing_len]
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i],len(unpacked_list[i].sequence)-35):] + unpacked_list[i+1].sequence
                    unpacked_list[i].sequence = unpacked_list[i].sequence[:unpacked_list[i].sequence.find(gg_opt[i],len(unpacked_list[i].sequence)-35)] + gg_opt[i] + "tgagacccg"
                elif unpacked_list[i].type == 'gRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg"
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):]
                elif unpacked_list[i].type == 'miRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg"
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):]
                else: #part is tRNA
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg"
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):]

            # If current part is last part, define reverse primer with right backbone overlap and forward primer ordinarily
            elif i == len(unpacked_list)-1:
                unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + bb_overlaps[-1] + 'tgagacccg')
                unpacked_list[i].sequence = unpacked_list[i].sequence + bb_overlaps[-1] + 'tgagacccg'

            # If current part is not first or last part, do ordinary computation
            else:
                if unpacked_list[i].type == 'pegRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76-max_annealing_len:unpacked_list[i].sequence.find(gg_opt[i],len(unpacked_list[i].sequence)-35)] + gg_opt[i] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i],len(unpacked_list[i].sequence)-35):] + unpacked_list[i+1].sequence[:max_annealing_len]
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i],len(unpacked_list[i].sequence)-35):] + unpacked_list[i+1].sequence
                    unpacked_list[i].sequence = unpacked_list[i].sequence[:unpacked_list[i].sequence.find(gg_opt[i],len(unpacked_list[i].sequence)-35)] + gg_opt[i] + "tgagacccg"
                elif unpacked_list[i].type == 'gRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg"
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):]
                elif unpacked_list[i].type == 'miRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg"
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):]
                else: #part is tRNA
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])] + gg_opt[i] + "tgagacccg"
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i]):]
        new_builds_list.append(unpacked_list)

    #Optimize primer Tm
    for unpacked_list in new_builds_list:
        for part in unpacked_list:
            
            # If both Tms are already below the range, find the nearest Gs to create GC clamp
            if mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0]:
                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-reverse(part.primer_forward).find('g')], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-reverse(part.primer_reverse).find('g')], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                
            # If one Tm is below the range, lower the other to make them most similar and find nearest Gs to create CG clamp
            elif mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0]:
                breakit=False
                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-reverse(part.primer_forward).find('g')], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                for j in range(len(part.primer_reverse),len(part.primer_reverse)-(max_annealing_len-18),-1):
                    if abs(part.primer_forward_tm-mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)) <=5 and part.primer_reverse[j-1] == 'g':
                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                        part.primer_reverse = part.primer_reverse[:j]
                        breakit=True
                if not breakit:
                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
            elif mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0]:
                breakit=False
                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-reverse(part.primer_reverse).find('g')], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                for i in range(len(part.primer_forward),len(part.primer_forward)-(max_annealing_len-18),-1):
                    if abs(part.primer_reverse_tm-mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)) <=5 and part.primer_forward[i-1] == 'g':
                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                        part.primer_forward = part.primer_forward[:i]
                        breakit=True
                if not breakit:
                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
            
            # If both Tms are within or above the range, lower them into the range and make them most similar and find nearest Gs to create GC clamp
            else:
                breakit=False
                # Do all of the above
                for i in range(len(part.primer_forward),len(part.primer_forward)-(max_annealing_len-18),-1):
                    for j in range(len(part.primer_reverse),len(part.primer_reverse)-(max_annealing_len-18),-1):
                        if abs(mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)-mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50))<=5 and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and part.primer_reverse[j-1] == 'g' and part.primer_forward[i-1] == 'g':
                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                            part.primer_forward = part.primer_forward[:i]
                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                            part.primer_reverse = part.primer_reverse[:j]
                            breakit=True
                        if breakit:
                            break
                    if breakit:
                        break

                # Do range and Gs        
                if not breakit:
                    for i in range(len(part.primer_forward),len(part.primer_forward)-(max_annealing_len-18),-1):
                        for j in range(len(part.primer_reverse),len(part.primer_reverse)-(max_annealing_len-18),-1):
                            if mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and part.primer_reverse[j-1] == 'g' and part.primer_forward[i-1] == 'g':
                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                                part.primer_forward = part.primer_forward[:i]
                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                                part.primer_reverse = part.primer_reverse[:j]
                                breakit=True
                            if breakit:
                                break
                        if breakit:
                            break

                # Do only Gs nearest to range
                if not breakit:
                    for i in range(len(part.primer_forward)-1-(max_annealing_len-18), len(part.primer_forward)):
                        for j in range(len(part.primer_reverse)-1-(max_annealing_len-18), len(part.primer_reverse)):
                            if part.primer_reverse[j-1] == 'g' and part.primer_forward[i-1] == 'g':
                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                                part.primer_forward = part.primer_forward[:i]
                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                                part.primer_reverse = part.primer_reverse[:j]
                                breakit=True
                            if breakit:
                                break
                        if breakit:
                            break
                    
                if not breakit:
                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-reverse(part.primer_forward).find('g')], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-reverse(part.primer_reverse).find('g')], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
    return new_builds_list[0],ftrs


# Define function to create a PTG from a list of parts
def PTGbldr(inserts):
    '''
    function which takes all desired parts of PTG GG assembly and gives out the respective inserts for PTG. During the 
    process, each part is appended with the same tRNA. The unit of part and tRNA are then treated as one insert.
    
    parts must be of the form [['name', 'type', 'sequence'],[...]]
    '''
    
    PTG_parts = []
    PTG_parts.append(Part('tRNA', 'tRNA', tRNA))
    
    ## Take each coding sequence of the output and append it with a tRNA to the list
    for c,prt in enumerate(inserts):
        if prt[1] == 'pegRNA':
            PTG_parts.append(Part(prt[0], prt[1], str(prt[2])))
            PTG_parts.append(Part('tRNA', 'tRNA', tRNA))
        elif prt[1] == 'gRNA':
            PTG_parts.append(Part(prt[0], prt[1], str(prt[2]) + tRNA))
        elif prt[1] == 'sRNA':
            PTG_parts.append(Part(prt[0], prt[1], str(prt[2]) + tRNA))
    
    return PTG_parts


# Execute computation
def runall(arr, tm_range=[52,72], max_ann_len=30, bb_overlaps=['tgcc','gttt'], additional_overhangs=[]):
    PTG = PTGbldr(arr)
    outpt,feat = scarless_gg(PTG, tm_range, max_ann_len, bb_overlaps, additional_overhangs)
    return outpt,feat


#full_test = [['gRNA0', 'gRNA', 'CGATTCGGCATAACGATCCC'+scaffld.upper(), 'f'],
#             ['pegRNA0', 'pegRNA', 'CGAGGAGCTGTTCACCGGGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCAGGATTGGCACCACCCCGGTGAACAGCTC', 'f'],
#             ['sRNA0', 'sRNA', 'GCTTAGATCGGATCCAAGCTA']]

#out,ftrs = runall(full_test)
#print(out[0].name)
