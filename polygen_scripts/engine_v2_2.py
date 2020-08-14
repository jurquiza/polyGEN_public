# Import all necessary packages

import numpy as np
import csv
import itertools
import re
import pandas as pd
import warnings
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#The following is a copyright notice by the authors of iBioCAD, from which
#substantial portions of this code were adopted and modified:

#Copyright (c) 2019 Scott Weisberg

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


#Define infrastructure for the following code
class Part:
    def __init__(self,name,type,sequence):
        self.name = name
        self.sequence = sequence
        self.type = type
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


# Define functions to ease the main computation
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

def flattn(list_2d):
    flattnd = []
    for subl in list_2d:
        if type(subl) is list:
            flattnd += [s for s in subl]
        else:
            flattnd.append(subl)
    return flattnd


#Optimizes the overhangs used in Golden Gate assembly
def golden_gate_optimization(parts_list,gg_overhangs):

    # Write all variable sequences in same order in a new list
    oh_list = []
    for x in parts_list:
        if x.type == 'pegRNA':
            oh_list.append([x.sequence[:20],x.sequence[96:]])
        elif x.type == 'gRNA':
            oh_list.append(x.sequence[:20])
        elif x.type == 'smRNA':
            oh_list.append(x.sequence[:-len(tRNA)])
        else: # if part is tRNA
            oh_list.append('plc')
    
    # Starting in the middle of the variable sequences and moving outwards, find overhang combinations
    for cov in range(2,max([int(np.ceil(np.true_divide(len(prt),2))) for prt in flattn(oh_list)])):
        
        cov_list = []
        for oh in oh_list:
            if type(oh) is list:
                o_list = []
                for o in oh:
                    if 2*cov >= len(o):
                        o_list.append(o)
                    else:
                        o_list.append(o[int(np.ceil(np.true_divide(len(o),2)))-cov:int(np.ceil(np.true_divide(len(o),2)))+cov])
                cov_list.append(o_list)
            else:
                if 2*cov >= len(oh):
                    cov_list.append(oh)
                else:
                    cov_list.append(oh[int(np.ceil(np.true_divide(len(oh),2)))-cov:int(np.ceil(np.true_divide(len(oh),2)))+cov])
        
        # Find possible overhang combinations
        golden_gate_overhangs=gg_overhangs
        seq_matches = []
        for x in range(len(parts_list)-1):
            seq_matches.append([])
            for overhang in golden_gate_overhangs:
                if parts_list[x].type == 'pegRNA':
                    if overhang in cov_list[x][1]:
                        seq_matches[x].append(oh_list[x][1][:oh_list[x][1].find(overhang)+4])
                elif parts_list[x].type == 'gRNA':
                    if parts_list[x+1].type == 'gRNA':
                        if overhang in cov_list[x+1]:
                            seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                    elif parts_list[x+1].type == 'pegRNA':
                        if overhang in cov_list[x+1][0]:
                            seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                    elif parts_list[x+1].type == 'smRNA':
                        if overhang in cov_list[x+1]:
                            seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                elif parts_list[x].type == 'smRNA':
                    if parts_list[x+1].type == 'gRNA':
                        if overhang in cov_list[x+1]:
                            seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                    elif parts_list[x+1].type == 'pegRNA':
                        if overhang in cov_list[x+1][0]:
                            seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                    elif parts_list[x+1].type == 'smRNA':
                        if overhang in cov_list[x+1]:
                            seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                else: # if part is tRNA
                    if parts_list[x+1].type == 'gRNA':
                        if overhang in cov_list[x+1]:
                            seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                    elif parts_list[x+1].type == 'pegRNA':
                        if overhang in cov_list[x+1][0]:
                            seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                    elif parts_list[x+1].type == 'smRNA':
                        if overhang in cov_list[x+1]:
                            seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                  
        combs = []
        for x in itertools.product(*seq_matches):
            combs.append(x)
        for comb in combs:
            if len([c[-4:] for c in comb]) == len(set([c[-4:] for c in comb])):
                return comb
    #if there are no possible combinations
    return None


# Set template sequences tRNA and gRNA scaffold
tRNA = 'AACAAAGCACCAGTGGTCTAGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA'
scaffld = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'


# Perform scarless Golden Gate assembly computation
def scarless_gg(parts_list, primer_tm_range=[52,60], max_annealing_len=30, bb_overlaps=['ggca','gttt'], additional_overhangs=[]):
    bb_overlaps = [i.lower() for i in bb_overlaps]
    additional_overhangs = [i.lower() for i in additional_overhangs]
    
    # Format parts list
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
            elif part.type == 'smRNA':
                ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+len(part.sequence)-len(tRNA), strand=1), type='smRNA'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+len(part.sequence)-len(tRNA), mmry+len(part.sequence), strand=1), type='tRNA'))
            mmry += len(part.sequence)
        
        # Iterate through overhang sets with increasing size until fitting one is found
        gg_opt = None
        breakit = False
        for p in range(10,51):
            for q in range(5):
                with open('overhangsets/setsof%s.csv'%p,'r') as f:
                    reader = csv.reader(f, delimiter=",")
                    sets = list(reader)[1:]

                temp = []
                for s in sets:
                    if len(s) != 0:
                        temp.append(s)

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
        
        if gg_opt is None:
            raise ValueError('The given combination of existing overhangs does not allow for an optimal overhang set.')
                
        #Modify sequences and design primers
        for i in range(len(unpacked_list)):
            
            # If current part is first part, define forward primer with left backbone overlap and reverse primer ordinarily
            if i == 0:
                unpacked_list[i].primer_forward = 'taggtctcc' + reverse_complement(bb_overlaps[0]) + unpacked_list[i].sequence[:max_annealing_len]
                unpacked_list[i].sequence = 'taggtctcc' + reverse_complement(bb_overlaps[0]) + unpacked_list[i].sequence
                if unpacked_list[i].type == 'pegRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76-max_annealing_len:unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + unpacked_list[i+1].sequence[:max_annealing_len]
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + unpacked_list[i+1].sequence
                    unpacked_list[i].sequence = unpacked_list[i].sequence[:unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + "tgagacccg"
                elif unpacked_list[i].type == 'gRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                elif unpacked_list[i].type == 'smRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                else: #part is tRNA
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
            
            # If current part is last part, define reverse primer with right backbone overlap and forward primer ordinarily
            elif i == len(unpacked_list)-1:
                unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + bb_overlaps[-1] + 'tgagacccg')
                unpacked_list[i].sequence = unpacked_list[i].sequence + bb_overlaps[-1] + 'tgagacccg'
            
            # If current part is not first or last part, do ordinary computation
            else:
                if unpacked_list[i].type == 'pegRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76-max_annealing_len:unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + unpacked_list[i+1].sequence[:max_annealing_len]
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + unpacked_list[i+1].sequence
                    unpacked_list[i].sequence = unpacked_list[i].sequence[:unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + "tgagacccg"
                elif unpacked_list[i].type == 'gRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                elif unpacked_list[i].type == 'smRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                else: #part is tRNA
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
        new_builds_list.append(unpacked_list)

    #Optimize primer Tm
    for unpacked_list in new_builds_list:
        for part in unpacked_list:
            
            # If both Tms are already below the range, find the nearest Gs
            if mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0]:
                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-re.search('[g,c]', reverse(part.primer_forward)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-re.search('[g,c]', reverse(part.primer_reverse)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                
            # If one Tm is below the range, lower the other to make them most similar and find nearest Gs
            elif mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0]:
                breakit=False
                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-re.search('[g,c]', reverse(part.primer_forward)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                for j in range(len(part.primer_reverse),len(part.primer_reverse)-(max_annealing_len-18),-1):
                    if abs(part.primer_forward_tm-mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)) <=5 and part.primer_reverse[j-1] in ['g','c']:
                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                        part.primer_reverse = part.primer_reverse[:j]
                        breakit=True
                if not breakit:
                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
            elif mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0]:
                breakit=False
                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-re.search('[g,c]', reverse(part.primer_reverse)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                for i in range(len(part.primer_forward),len(part.primer_forward)-(max_annealing_len-18),-1):
                    if abs(part.primer_reverse_tm-mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)) <=5 and part.primer_forward[i-1] in ['g','c']:
                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                        part.primer_forward = part.primer_forward[:i]
                        breakit=True
                if not breakit:
                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
            
            # If both Tms are within or above the range, lower them into the range and make them most similar and find nearest Gs
            else:
                breakit=False
                # Do all of the above
                for i in range(len(part.primer_forward),len(part.primer_forward)-(max_annealing_len-18),-1):
                    for j in range(len(part.primer_reverse),len(part.primer_reverse)-(max_annealing_len-18),-1):
                        if abs(mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)-mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50))<=5 and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and part.primer_reverse[j-1] in ['g','c'] and part.primer_forward[i-1] in ['g','c']:
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
                            if mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and part.primer_reverse[j-1] in ['g','c'] and part.primer_forward[i-1] in ['g','c']:
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
                            if part.primer_reverse[j-1] in ['g','c'] and part.primer_forward[i-1] in ['g','c']:
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
                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-re.search('[g,c]', reverse(part.primer_forward)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-re.search('[g,c]', reverse(part.primer_reverse)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
    return new_builds_list[0],ftrs


def pegbldr(sequence, edits, mode='PE2'):
    '''
    function that builds the necesary pegRNA and gRNA to accurately introduce a point mutation using the prime editor 
    system.
    
    sequence is the original 5' -> 3' sequence strand that is to be edited
    edits is an array containing all desired edits in the form [[index of edit, new bases, type of edit],[...]]
    '''
    
    sequence = sequence.upper()
    plc_seq = sequence[:]
    plc_rev_comp = reverse_complement(plc_seq)
    
    out = []
    for c,edt in enumerate(edits):
        
        seq = plc_seq[:]
        rev_comp = plc_rev_comp[:]
        
        inds = [int(i) for i in str(edt[0]).split(',')]
        changes = str(edt[1]).upper().split(',')

        
                ## Define pegRNA for editing
        
        # Find all PAM motifs in forward and reverse strand in respective possible regions
        pegPAMrgn_forw = seq[:inds[0]+6]                         # Define the region where PAM could possibly lie in forw strand
        pegPAMs_forw = re.finditer(r'(?=(.GG))', pegPAMrgn_forw) # Find all PAMs in region
        pegPAMrgn_rev = rev_comp[:len(seq)-inds[-1]+6]
        pegPAMs_rev = re.finditer(r'(?=(.GG))', pegPAMrgn_rev)
        pegPAMs_forw = list(pegPAMs_forw)
        pegPAMs_rev = list(pegPAMs_rev)
        pegPAMs_forw = [i.start() for i in pegPAMs_forw]
        pegPAMs_rev = [i.start() for i in pegPAMs_rev]
        
        # Check if usable PAMs are present
        if pegPAMs_forw == pegPAMs_rev == []:
            raise ValueError("There are no usable PAM motifs around the edit")
        
        # Check if there are PAMs in only one strand
        elif pegPAMs_forw == []:
            pegPAM_rev = pegPAMs_rev[np.argmin([abs(i-(len(seq)-2-inds[-1])) for i in pegPAMs_rev])]
            pegPAM = pegPAM_rev
            pegPAM_strand = 'r'
            seq = rev_comp[:]
            rev_comp = sequence[:]
            inds = [len(seq)-i-1 for i in inds] # Recalculate mutation indices for rev strand. -1 because len(x) == ind(x[-1])+1
            changes = [complement(i) for i in changes]
        elif pegPAMs_rev == []:
            pegPAM_forw = pegPAMs_forw[np.argmin([abs(i-1-inds[0]) for i in pegPAMs_forw])]
            pegPAM = pegPAM_forw
            pegPAM_strand = 'f'
        
        # If there are PAMs in both strands, choose the one closest to edit
        else:
            forw_min = min([abs(i-1-inds[0]) for i in pegPAMs_forw])
            rev_min = min([abs(i-(len(seq)-2-inds[-1])) for i in pegPAMs_rev])
            pegPAM_forw = pegPAMs_forw[np.argmin([abs(i-1-inds[0]) for i in pegPAMs_forw])]
            pegPAM_rev = pegPAMs_rev[np.argmin([abs(i-(len(seq)-2-inds[-1])) for i in pegPAMs_rev])]
            if forw_min <= rev_min: # Of the closest PAM of each strand, which is closest to mutation? -1 because len(x) == ind(x[-1])+1
                pegPAM = pegPAM_forw
                pegPAM_strand = 'f'
            else:
                pegPAM = pegPAM_rev
                pegPAM_strand = 'r'
                seq = rev_comp[:]
                rev_comp = sequence[:]
                inds = [len(seq)-i-1 for i in inds]              # Recalculate mutation indices for rev strand. -1 because len(x) == ind(x[-1])+1
                changes = [complement(i) for i in changes]
        
        if inds[-1]-pegPAM > 30:
                warnings.warn("There is no PAM motif in +/- 30 nt proximity of edit " + str(c))

        if pegPAM-20 < 0:
            raise ValueError("The provided sequence does not cover enough area around the edit")
        else:
            pegspacer = seq[pegPAM-20:pegPAM]                    # Spacer should be 20 nt in length and end at PAM

        
        # Calculate RT templates depending on type of edit
        if edt[2] == 'mut':                                      # Check if edit is point mutation
            
            pre_RT_len = max([13, inds[-1]-(pegPAM-3)])          # Set default length of RT-template to 13 (recommended by Anzalone et al. 2019) or until edit if further
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start() # From default length find next D
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]      # Retrieve RT-template
            RT_templ = [i for i in RT_templ]
            
            for c_pm,pm in enumerate(inds):                      # If several point mutations, go through each separately
                RT_templ[int(pm)-(pegPAM-3)] = changes[c_pm]     # Include edit in RT-template
                
        elif edt[2] == 'ins':
            
            pre_RT_len = max([13, inds[-1]-(pegPAM-3)+6])        # template should have additional length 5' of insert to ensure binding
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start()
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]
            RT_templ = [i for i in RT_templ]
            
            for c_pm,pm in enumerate(inds):
                RT_templ = RT_templ[:int(pm)-(pegPAM-3)] + [l for l in changes[c_pm]] + RT_templ[int(pm)-(pegPAM-3):]
                
        elif edt[2] == 'del':
            
            len_deltns = 0
            for c_pm,pm in enumerate(inds):
                len_deltns += int(changes[c_pm])-pm
            pre_RT_len = max([13, int(changes[-1])-(pegPAM-3)+6+len_deltns])
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start()
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]
            RT_templ = [i for i in RT_templ]
            
            mmry = 0
            for c_pm,pm in enumerate(inds):
                del RT_templ[pm-(pegPAM-3)-mmry:int(changes[c_pm])-(pegPAM-3)-mmry]
                mmry += int(changes[c_pm])-pm
                
        RT_templ = ''.join(RT_templ)
        
        RT_templ = reverse_complement(RT_templ)        # RT-template must be in opposite direction
                                                               
        PBS = reverse_complement(seq[pegPAM-16:pegPAM-3]) # PBS must be in opposite direction
        
        pegRNA = pegspacer + scaffld + RT_templ + PBS
        
        out.append(['pegRNA'+str(c), 'pegRNA', pegRNA, pegPAM_strand])
        
        ## Define gRNA for PE3
        if mode == 'PE3':
            
            old = []
            seq = [i for i in seq]
            for i in range(len(inds)):
                old.append(seq[inds[i]])
                seq[inds[i]] = changes[i]
            seq = ''.join(seq)
            
            gPAMrgn = reverse_complement(seq[pegPAM:])     # gRNA must bind to other strand somewhere after pegPAM
            gPAMs = re.finditer(r'(?=(.GG))', gPAMrgn)     # find all PAMs in that region
            gPAMs = np.array([i.start() for i in gPAMs])
            gPAM = gPAMs[abs(np.array(gPAMs)-(len(gPAMrgn)-44)).argmin()] # Find PAM closest to 47 nt upstream of pegPAM nick
                                                                          # (the gRNA nick should be 50 nt upstream of pegPAM nick)

            gspacer = gPAMrgn[gPAM-20:gPAM]                # Define spacer

            gRNA = gspacer

            out.append(['gRNA'+str(c), 'gRNA', gRNA, pegPAM_strand])

            seq = [i for i in seq]
            for i in range(len(inds)):
                seq[inds[i]] = old[i]
            seq = ''.join(seq)
    
    return out


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
            PTG_parts.append(Part(prt[0], prt[1], str(prt[2]) + scaffld + tRNA))
        elif prt[1] == 'smRNA':
            PTG_parts.append(Part(prt[0], prt[1], str(prt[2]) + tRNA))
    
    return PTG_parts


# Execute computation
def runall(arr, tm_range=[52,72], max_ann_len=30, bb_overlaps=['tgcc','gttt'], additional_overhangs=[]):
    PTG = PTGbldr(arr)
    outpt,feat = scarless_gg(PTG, tm_range, max_ann_len, bb_overlaps, additional_overhangs)

    oligos = []
    full_sequence = ''
    for c,o in enumerate(outpt):
        oligos.append(o.primer_forward)
        oligos.append(o.primer_reverse)
        if c == 0:
            full_sequence += o.sequence[:-13]
        elif c == len(outpt)-1:
            full_sequence += o.sequence[9:]
        else:
            full_sequence += o.sequence[9:-13]

    return outpt,full_sequence# ,oligos,feat

