# Import all necessary packages

import numpy as np
import csv
import itertools
import re
import pandas as pd
import json
import copy
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC



class InvalidUsage(Exception):
    '''
    Handler for user-interpretable errors
    
    :param message: A user-interpretable error message
    :type message: str
    :param status_code: The status code of the server, defaults to 400
    :type status_code: int, optional
    :param payload: Additional information on the error location, defaults to None
    :type payload: dict, optional
    '''

    def __init__(self, message, status_code=400, payload=None):
        '''inits InvalidUsage'''
        
        Exception.__init__(self)
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        '''
        Turns error object into dictionary
        
        :return: A dictionary containing the error message under key 'message' and the payload under the respective keys
        :rtype: dict
        '''
        
        rv = dict(self.payload or ())
        rv['message'] = self.message
        return rv

#Copyright (c) 2019 Scott Weisberg
class Part:
    '''
    Object used for defining sequence characteristics
    
    :param name: Sequence name
    :type name: str
    :param type: Sequence type
    :type type: type
    :param sequence: A nucleotide sequence consisting of [ACGTacgt]
    :type sequence: str
    '''
    
    def __init__(self,name,type,sequence):
        '''inits Part'''
        
        self.name = name
        self.sequence = sequence
        self.type = type
    
    def to_json(self):
        '''transforms data to json format'''
        
        return {
            "name": self.name,
            "sequence": self.sequence,
            "type": self.type,
            "primer_forward": self.primer_forward,
            "primer_reverse": self.primer_reverse,
            "primer_forward_tm": self.primer_forward_tm,
            "primer_reverse_tm": self.primer_reverse_tm
        }
    
    @classmethod
    def from_json(cls, json_data):
        return cls(json_data['name'],
                   json_data['type'],
                   json_data['sequence']
                   )
        
    primer_forward = ""
    primer_reverse = ""
    primer_forward_tm = 0
    primer_reverse_tm = 0
    
class Polycistron:
    
    
    def __init__(self):
        '''inits Polycistron'''
    
    sequence = ''
    parts = []
    features = []
    warning = None


# Define functions to ease the main computation
def polyToJson(poly):
    '''transforms data to json format''' 
    
    p = copy.deepcopy(poly)
    
    ftrs = [vars(feat) for feat in p.features]
    for feat in ftrs:
        feat['location'] = vars(feat['location'])
        
    return {
        "sequence": p.sequence,
        "parts": [p.to_json() for p in p.parts],
        "features": ftrs,
        "warning": p.warning
    }

def reverse_complement(sequence):
    '''
    Finds the reverse complement of a provided sequence
    
    :param sequence: The nucleic acid sequence consisting of [ACGTacgt]
    :type sequence: str
    
    :return: The reverse complement of the provided sequence
    :rtype: str
    '''
    
    rev_comp = ""
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        rev_comp = Watson_Crick[base] + rev_comp
    return rev_comp

def complement(sequence):
    '''
    Finds the complement of a provided sequence
    
    :param sequence: A nucleic acid sequence consisting of [ACGTacgt]
    :type sequence: str
    
    :return: The complement of the provided sequence
    :rtype: str
    '''
    
    comp = ''
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        comp += Watson_Crick[base]
    return comp

def reverse(sequence):
    '''
    Finds the reverse of a provided sequence
    
    :param sequence: A nucleic acid sequence consisting of [ACGTacgt]
    :type sequence: str
    
    :return: The reverse of the provided sequence
    :rtype: str
    '''
    
    rev = ''
    for base in sequence:
        rev = base + rev
    return rev

def flattn(list_2d):
    '''
    Flattens a 2D list of lists into one 1D list
    
    :param list_2d: 2D list
    :type list_2d: list
    
    :return: 1D list with all sublists of the input concatenated
    :rtype: list
    '''
    
    flattnd = []
    for subl in list_2d:
        if type(subl) is list:
            flattnd += [s for s in subl]
        else:
            flattnd.append(subl)
    return flattnd

def Diff(li1, li2):
    '''
    Finds all elements of two lists, which are unique to one of the two
    
    :param li1: First list
    :type li1: list
    :param li2: Second list
    :type li2: list
    
    :return: Differential list
    :rtype: list
    '''
    
    return (list(list(set(li1)-set(li2)) + list(set(li2)-set(li1))))
    
    
# Set template sequences tRNA and gRNA scaffold
tRNA = 'aacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca'
scaffld = 'gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc'
DR = 'aatttctactgttgtagat'


#Optimize the overhangs used in Golden Gate assembly
def golden_gate_optimization(parts_list, free_overhangsets, poltype_opt='ptg'):
    '''
    Finds linkers from optimal linker sets in a provided parts list
    
    :param parts_list: An array of objects of class part
    :type parts_list: list
    :param free_overhangsets: An array of arrays containing 4 bp linkers
    :type free_overhangsets: list
    :param poltype_opt: The type of polycistronic architecture to use. Must be one of 'ptg' or 'ca', defaults to 'ptg'
    :type poltype_opt: str, optional
    
    :return: A tuple of sequences indicating the linker between two adjacent parts. The sequence spans from the start of the respective variable sequence to the end of the determined 4 bp linker. 
    :rtype: tuple
    '''
    
    for p in parts_list:
        p.sequence = p.sequence.lower()
        
    # Write all variable sequences in same order in a new list
    oh_list = []
    if poltype_opt=='ptg':
        for x in parts_list:
            if x.type == 'pegRNA':
                oh_list.append([x.sequence[:x.sequence.find(scaffld)],x.sequence[x.sequence.find(scaffld)+len(scaffld):]])
            elif x.type == 'gRNA':
                oh_list.append(x.sequence[:x.sequence.find(scaffld)])
            elif x.type == 'smRNA':
                oh_list.append(x.sequence[:-len(tRNA)])
            else: # if part is tRNA
                oh_list.append('plc')
    elif poltype_opt=='ca':
        for x in parts_list:
            oh_list.append(x.sequence[x.sequence.find(DR)+len(DR):])
    
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
        if poltype_opt=='ptg':
            for golden_gate_overhangs in free_overhangsets:
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
        elif poltype_opt=='ca':
            for s in free_overhangsets:
                golden_gate_overhangs=s
                seq_matches = []
                for x in range(len(parts_list)-1):
                    seq_matches.append([])
                    for overhang in golden_gate_overhangs:
                        if overhang in cov_list[x]:
                            seq_matches[x].append(oh_list[x][:oh_list[x].find(overhang)+4])
                            
        #Copyright (c) 2019 Scott Weisberg
        combs = []
        for x in itertools.product(*seq_matches):
            combs.append(x)
        for comb in combs:
            if len([c[-4:] for c in comb]) == len(set([c[-4:] for c in comb])):
                return comb
    #if there are no possible combinations
    return None


#Copyright (c) 2019 Scott Weisberg
# Perform scarless Golden Gate assembly computation with provided parts
def scarless_gg(parts_list, tm_range=[55,65], max_ann_len=30, bb_overlaps=['ggca','aaac'], additional_overhangs=[], poltype='ptg', enzm='bsai'):
    '''
    Uses a list of desired parts and additional arguments to compute a corresponsing PTG. Returns a list of newly computed parts and their primers which can be used to generate the PTG.
    
    :param parts_list: An array of objects of class part, to be included in a PTG
    :type parts_list: list
    :param tm_range: An array of integers of length 2. Temperature range to aim for during primer optimization, defaults to [55,65]
    :type tm_range: list, optional
    :param max_ann_len: The maximal annealing length of the static part of the primer, defaults to 30
    :type max_ann_len: int, optional
    :param bb_overlaps: Linkers of the destination plasmid flanking the final PTG, defaults to ['ggca','aaac']
    :type bb_overlaps: list, optional
    :param additional_overhangs: Additional linkers in the destination plasmid, defaults to []
    :type additional_overhangs: list, optional
    :param poltype: Type of polycistronic architecture to use. Must be one of 'ptg' or 'ca', defaults to 'ptg'
    :type poltype: str, optional
    :param enzm: Type II restriction enzyme to use for the Golden Gate assembly. Defaults to 'bsai'
    :type enzm: str, optional
    
    :return: Returns three objects. (1) A list of newly computed parts of class Part, (2) A list of features of class SeqRecord, (3) An error or warning message
    :rtype: list, list, str
    '''
    
    polycistron = Polycistron()
    polycistron.features = [] # must overwrite with empty list because features list would accumulate across runs in same session through append command
    
    bb_overlaps = [i.lower() for i in bb_overlaps]
    additional_overhangs = [i.lower() for i in additional_overhangs]
    enzms={'bsai': ['gaggtctcg', 'cgagacctc'], 'bsmbi': ['tgcgtctca', 'tgagacgca'], 'btgzi': ['ctgcgatggagtatgtta', 'taacatactccatcgcag'], 'bpii': ['ttgaagactt', 'aagtcttcaa']} #templates found in pUU080 (bsai), pUPD2 (bsmbi), Ortega-Escalante et al. 2018 (btgzi), Kun (bpii)

    
    # Go through parts and write all known annotations into list   
    mmry = len(enzms[enzm][0])+4
    polycistron.sequence = enzms[enzm][0] + reverse_complement(bb_overlaps[0])
    for part in parts_list:
        part.sequence = part.sequence.lower()
        polycistron.sequence += part.sequence
        polycistron.features.append(SeqFeature(FeatureLocation(mmry, mmry+len(part.sequence), strand=1), type=part.type))
        if part.type == 'pegRNA':
            polycistron.features.append(SeqFeature(FeatureLocation(mmry, mmry+part.sequence.find(scaffld), strand=1), type='spacer'))
            polycistron.features.append(SeqFeature(FeatureLocation(mmry+part.sequence.find(scaffld)+len(scaffld), mmry+len(part.sequence)-13, strand=1), type='RT template'))
            polycistron.features.append(SeqFeature(FeatureLocation(mmry+len(part.sequence)-13, mmry+len(part.sequence), strand=1), type='PBS'))
        elif part.type == 'gRNA':
            polycistron.features.append(SeqFeature(FeatureLocation(mmry, mmry+part.sequence.find(scaffld), strand=1), type='spacer'))
            polycistron.features.append(SeqFeature(FeatureLocation(mmry+part.sequence.find(tRNA), mmry+len(part.sequence), strand=1), type='tRNA'))
        elif part.type == 'smRNA':
            polycistron.features.append(SeqFeature(FeatureLocation(mmry, mmry+len(part.sequence)-len(tRNA), strand=1), type='smRNA'))
            polycistron.features.append(SeqFeature(FeatureLocation(mmry+len(part.sequence)-len(tRNA), mmry+len(part.sequence), strand=1), type='tRNA'))
        if part.type == 'crRNA':
            polycistron.features.append(SeqFeature(FeatureLocation(mmry+part.sequence.find(DR)+len(DR), mmry+len(part.sequence), strand=1), type='spacer'))
        mmry += len(part.sequence)
    polycistron.sequence += bb_overlaps[1] + enzms[enzm][1]

    # Iterate through overhang sets with increasing size until fitting one is found
    gg_opt = None
    breakit = False
    existing_overhangs = additional_overhangs+bb_overlaps
    for p in range(10,51):
        free_overhangsets = []
        for q in range(5):
            with open('overhangsets/setsof%s.csv'%p,'r') as f:
                reader = csv.reader(f, delimiter=",")
                sets = list(reader)[1:]
                
            temp = []
            for s in sets:
                if len(s) != 0:
                    temp.append(s)

            # Only grab sets that include all existing overhangs and delete the existing from the set
            st = [x.lower() for x in temp[q]]
            if all(i in st or reverse_complement(i) in st for i in existing_overhangs):
                free_overhangsets.append([i for i in [x.lower() for x in temp[q]] if i not in additional_overhangs+bb_overlaps])
            else:
                continue
           
        if free_overhangsets:
            gg_opt = golden_gate_optimization(parts_list, free_overhangsets, poltype)
        if gg_opt is not None:
            breakit = True
        if breakit:
            break
        
    # No sets include all existing overhangs
    if gg_opt is None:
            
        overhangs_sublist = []
        for l in range(len(existing_overhangs)-1,0,-1):
            overhangs_subsublist = list(itertools.combinations(existing_overhangs,l))
            overhangs_sublist += overhangs_subsublist
            
        for sublist in overhangs_sublist:
            for p in range(10,51):
                free_overhangsets = []
                for q in range(5):
                    with open('overhangsets/setsof%s.csv'%p,'r') as f:
                        reader = csv.reader(f, delimiter=",")
                        sets = list(reader)[1:]

                    temp = []
                    for s in sets:
                        if len(s) != 0:
                            temp.append(s)

                    # Only grab sets that include all linkers in the current sublist and delete the existing from the set
                    if all(i in [x.lower() for x in temp[q]] for i in sublist):
                        free_overhangsets.append([i for i in [x.lower() for x in temp[q]] if i not in additional_overhangs+bb_overlaps])
                    else:
                        continue

                if free_overhangsets:
                    gg_opt = golden_gate_optimization(parts_list, free_overhangsets, poltype)
                if gg_opt is not None:
                    breakit = True
                    
                    exist = ','.join(sublist)
                    nexist = ','.join(Diff(sublist, existing_overhangs))
                    polycistron.warning = 'The given combination of existing overhangs is not compatible with an optimal overhang set. Found the set including the largest possible fraction of existing overhangs ('+exist+'). The following overhangs were disregarded: '+nexist+'. There might be interference between overhangs.'
                if breakit:
                    break
            if breakit:
                break
        
    if gg_opt is None:
        raise InvalidUsage("No combination of optimal linkers possible for the provided existing linkers", status_code=400, payload={'pge': 'sequence.html', 'box': 'link'})


    #Modify sequences and design primers
    if poltype=='ptg':
        for i in range(len(parts_list)):
            
            # If current part is first part, define forward primer with left backbone overlap and reverse primer ordinarily
            if i == 0:
                parts_list[i].primer_forward = enzms[enzm][0] + reverse_complement(bb_overlaps[0]) + parts_list[i].sequence[:max_ann_len]
                parts_list[i].sequence = enzms[enzm][0] + reverse_complement(bb_overlaps[0]) + parts_list[i].sequence
                if parts_list[i].type == 'pegRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76-max_ann_len:parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + parts_list[i+1].sequence[:max_ann_len]
                    parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + parts_list[i+1].sequence
                    parts_list[i].sequence = parts_list[i].sequence[:parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + enzms[enzm][1]
                elif parts_list[i].type == 'gRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                elif parts_list[i].type == 'smRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                else: #part is tRNA
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
            
        # If current part is last part, define reverse primer with right backbone overlap and forward primer ordinarily
            elif i == len(parts_list)-1:
                parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + bb_overlaps[-1] + enzms[enzm][1])
                parts_list[i].sequence = parts_list[i].sequence + bb_overlaps[-1] + enzms[enzm][1]
            
        # If current part is not first or last part, do ordinary computation
            else:
                if parts_list[i].type == 'pegRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76-max_ann_len:parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + parts_list[i+1].sequence[:max_ann_len]
                    parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + parts_list[i+1].sequence
                    parts_list[i].sequence = parts_list[i].sequence[:parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + enzms[enzm][1]
                elif parts_list[i].type == 'gRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                elif parts_list[i].type == 'smRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                else: #part is tRNA
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                            
    elif poltype=='ca':
        for i in range(len(parts_list)):
            
            # If current part is first part, define forward primer with left backbone overlap and reverse primer ordinarily
            if i == 0:
                parts_list[i].primer_forward = enzms[enzm][0] + reverse_complement(bb_overlaps[0]) + DR
                parts_list[i].sequence = enzms[enzm][0] + reverse_complement(bb_overlaps[0]) + parts_list[i].sequence
                
                parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[parts_list[i].sequence.find(DR):parts_list[i].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])-4:] + DR
                parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])-4:] + parts_list[i+1].sequence
                parts_list[i].sequence = parts_list[i].sequence[:parts_list[i].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    
            # If current part is last part, define reverse primer with right backbone overlap and forward primer ordinarily
            elif i == len(parts_list)-1:
                parts_list[i].primer_reverse = reverse_complement(DR + bb_overlaps[-1] + enzms[enzm][1])
                parts_list[i].sequence = parts_list[i].sequence + bb_overlaps[-1] + enzms[enzm][1]
                
            # If current part is not first or last part, do ordinary computation
            else:
                parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[parts_list[i].sequence.find(DR):parts_list[i].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])-4:] + DR
                parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])-4:] + parts_list[i+1].sequence
                parts_list[i].sequence = parts_list[i].sequence[:parts_list[i].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                            
    for o in parts_list:
        o.primer_forward = o.primer_forward[:len(enzms[enzm][0])].lower() + o.primer_forward[len(enzms[enzm][0]):len(enzms[enzm][0])+4].upper() + o.primer_forward[len(enzms[enzm][0])+4:].lower()
        o.primer_reverse = o.primer_reverse[:len(enzms[enzm][1])].lower() + o.primer_reverse[len(enzms[enzm][1]):len(enzms[enzm][1])+4].upper() + o.primer_reverse[len(enzms[enzm][1])+4:].lower()
    
    if poltype == 'ca':
        polycistron.parts = parts_list
        return polycistron # If CA, the primer optimization step can be skipped, since the DR is very short

    polycistron.parts = parts_list

    
    ## Optimise primer Tm    
    for c,prmr in enumerate(flattn([[part.primer_forward, part.primer_reverse] for part in polycistron.parts])):
    
        prmrRest = prmr[:-max_ann_len]
        prmrRestLen = len(prmrRest)
    
        scr = [1 if i in ['C','G','g','c'] else 0 for i in prmr] # Extract score from sequence for finding GC-clamps
        optns = [prmr[:i] for i in range(len(prmr),len(prmr)-(max_ann_len-18+1), -1) if sum(scr[i-5:i]) in [1,2]] # Find all GC-clamps
        checkStatic = [optn for optn in optns if tm_range[0] <= mt.Tm_NN(optn[prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= tm_range[1]] # Find all GC-clamp options with static Tm inside provided range
        checkWhole = [optn for optn in checkStatic if tm_range[0] <= mt.Tm_NN(optn, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= tm_range[1]] # Find all GC-clamp options with whole Tm inside provided range
        
        if checkWhole: # If there are GC-clamp options with whole Tm, use the shortest one to cut costs
            fin = checkWhole[-1]
        elif checkStatic: # If there are no GC-clamps with whole Tm but such with static Tm, use the shortest, to (1) cut costs and (2) get the Tm of the whole primer as close as possible to range.
            fin = checkStatic[-1]
        elif optns and mt.Tm_NN(optns[-1][prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) > tm_range[1]: # If all GC-clamps are above the range, use lowest one
            fin = optns[-1]
        elif optns and mt.Tm_NN(optns[0][prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < tm_range[0]: # If all GC-clamps are below the range, use highest one
            fin = optns[0]
        else: # If there are no GC-clamps, find the primer with static Tm closest to bottom of range for cutting costs
            countUp = prmrRestLen+18
            lastResort = prmr[prmrRestLen:countUp]
            while mt.Tm_NN(lastResort, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < tm_range[0] and prmrRestLen + len(lastResort) < len(prmr):
                countUp += 1
                lastResort = prmr[prmrRestLen:countUp]
            fin = prmr[:countUp]
        
        if c%2 == 0: # Assign the found primer and corresponding static Tm to the respective part
            polycistron.parts[int(np.floor(c/2))].primer_forward = fin
            polycistron.parts[int(np.floor(c/2))].primer_forward_tm = mt.Tm_NN(fin[prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
        else:
            polycistron.parts[int(np.floor(c/2))].primer_reverse = fin
            polycistron.parts[int(np.floor(c/2))].primer_reverse_tm = mt.Tm_NN(fin[prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
            
    return polycistron


def pegbldr(sequence, edits, mode='PE2'):
    '''
    Builds the necesary pegRNA and gRNA to accurately introduce a point mutation using the prime editor system.
    
    :param sequence: The original 5' -> 3' nucleic acid sequence consisting of [ACGTacgt] that is to be edited
    :type sequence: str
    :param edits: An array containing all desired edits in the form [[index of edit, new bases, type of edit],[...]]
    :type edits: list
    :param mode: The type of prime editing mode to design for. One of 'PE2' or 'PE3', defaults to 'PE3'
    :type mode: str, optional
    
    :return: A list of lists containing for each computed guide RNA the name, type, sequence and strand specifications
    :rtype: list
    '''
    msg = None
    
    if sequence == '':
        raise InvalidUsage("No sequence input", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})
    elif re.search(r'^[ACGTacgt]*$', sequence) is None:
        raise InvalidUsage("Invalid sequence input ", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})
    
    if edits == [['']]:
        raise InvalidUsage("No edits input", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
    for e in edits:
        if len(e) != 3:
            raise InvalidUsage("Invalid edits input syntax", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif len(e[0].split(',')) != len(e[1].split(',')):
            raise InvalidUsage("Inconsistent number of edits", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif e[2] not in ['mut', 'del', 'ins']:
            raise InvalidUsage("Invalid edit type", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif e[2] == 'mut' and (any([re.search(r'^[0-9]*$', n) is None for n in e[0].split(',')]) or any([re.search(r'^[ACGTacgt]*$', b) is None for b in e[1].split(',')])):
            raise InvalidUsage("Invalid specifications for edit type \'mut\'", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif e[2] == 'del' and any([re.search(r'^[0-9]*$', n) is None for n in e[0].split(',')+e[1].split(',')]):
            raise InvalidUsage("Invalid specifications for edit type \'del\'", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif e[2] == 'ins' and (any([re.search(r'^[0-9]*$', n) is None for n in e[0].split(',')]) or any([re.search(r'^[ACGTacgt]*$', b) is None for b in e[1].split(',')])):
            raise InvalidUsage("Invalid specifications for edit type \'ins\'", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
    
    sequence = sequence.upper()
    plc_seq = sequence[:]
    plc_rev_comp = reverse_complement(plc_seq)
    
    out = []
    for c,edt in enumerate(edits):
        
        seq = plc_seq[:]
        rev_comp = plc_rev_comp[:]
        
        inds = [int(i) for i in str(edt[0]).split(',')]
        changes = str(edt[1]).upper().split(',')
        if edt[2] == 'del':
            changes = [str(int(chng)+1) for chng in changes] # deletion stretches should include the ending index for intuitive usability

        
                ## Define pegRNA for editing
        
        allInd = inds + [int(chng) for chng in changes if edt[2] == 'del']
        maxInd = max(allInd)
        minInd = min(allInd)
        
        # Find all PAM motifs in forward and reverse strand in respective possible regions
        pegPAMrgn_forw = seq[:minInd+6]                       # Define the region where PAM could possibly lie in forw strand
        pegPAMs_forw = re.finditer(r'(?=(.GG))', pegPAMrgn_forw) # Find all PAMs in region
        pegPAMrgn_rev = rev_comp[:len(seq)-1-maxInd+6]
        pegPAMs_rev = re.finditer(r'(?=(.GG))', pegPAMrgn_rev)
        pegPAMs_forw = list(pegPAMs_forw)
        pegPAMs_rev = list(pegPAMs_rev)
        pegPAMs_forw = [i.start() for i in pegPAMs_forw]
        pegPAMs_rev = [i.start() for i in pegPAMs_rev]
        
        # Check if usable PAMs are present
        if pegPAMs_forw == pegPAMs_rev == []:
            raise InvalidUsage("There are no usable PAM motifs around the edit", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})
        
        # Check if there are PAMs in only one strand
        elif pegPAMs_forw == []:
            pegPAM_rev = pegPAMs_rev[np.argmin([abs(i-(len(seq)-2-max(inds))) for i in pegPAMs_rev])]
            pegPAM = pegPAM_rev
            pegPAM_strand = 'r'
            seq = rev_comp[:]
            rev_comp = sequence[:]
            inds = [len(seq)-i-1 for i in inds] # Recalculate mutation indices for rev strand. -1 because len(x) == ind(x[-1])+1
            changes = [complement(i) for i in changes]
        elif pegPAMs_rev == []:
            pegPAM_forw = pegPAMs_forw[np.argmin([abs(i-1-min(inds)) for i in pegPAMs_forw])]
            pegPAM = pegPAM_forw
            pegPAM_strand = 'f'
        
        # If there are PAMs in both strands, choose the one closest to edit. inf is used to make sure the list is not empty but will never be chosen.
        else:
            forw_min = min([abs(i-1-inds[0]) for i in pegPAMs_forw] + [np.inf])
            rev_min = min([abs(i-(len(seq)-2-inds[-1])) for i in pegPAMs_rev] + [np.inf])
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
                if edt[2] == 'mut':
                    inds = [len(seq)-i-1 for i in inds]              # Recalculate mutation indices for rev strand. -1 because len(x) == ind(x[-1])+1
                    changes = [complement(i) for i in changes]
                elif edt[2] == 'ins':
                    inds = [len(seq)-i for i in inds]
                    changes = [complement(i) for i in changes]
                else:
                    changes = [int(chng) for chng in changes] # Change type of changes to int
                    changes_plc = changes[:] # create placeholder
                    changes = [len(seq)-i-1 for i in inds] # switch changes and inds
                    inds = [len(seq)-i-1 for i in changes_plc]
                                  
        if max(inds)-pegPAM > 30:
                msg = "There was no PAM motif in +/- 30 nt proximity of edit " + str(c) + ". Used the nearest one, which was " + str(max(inds)-pegPAM) + " bp away."

        if pegPAM-20 < 0:
            raise InvalidUsage("The provided sequence does not cover enough area around the edit", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})
        else:
            pegspacer = seq[pegPAM-20:pegPAM]                    # Spacer should be 20 nt in length and end at PAM
        PBS = reverse_complement(seq[pegPAM-16:pegPAM-3]) # PBS must be in opposite direction

        
        # Calculate RT templates depending on type of edit
        if edt[2] == 'mut':                                      # Check if edit is point mutation
            
            pre_RT_len = max([13, max(inds)-(pegPAM-3)])          # Set default length of RT-template to 13 (recommended by Anzalone et al. 2019) or until edit if further
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start() # From default length find next D
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]      # Retrieve RT-template
            RT_templ = [i for i in RT_templ]
            
            for c_pm,pm in enumerate(inds):                      # If several point mutations, go through each separately
                RT_templ[int(pm)-(pegPAM-3)] = changes[c_pm]     # Include edit in RT-template
                
        elif edt[2] == 'ins':
            
            pre_RT_len = max([13, max(inds)-(pegPAM-3)+6])        # template should have additional length 5' of insert to ensure binding
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start()
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]
            RT_templ = [i for i in RT_templ]
            for c_pm,pm in enumerate(inds):
                RT_templ = RT_templ[:int(pm)-(pegPAM-3)] + [l for l in changes[c_pm]] + RT_templ[int(pm)-(pegPAM-3):]
                
        elif edt[2] == 'del':
            
            changes = [int(chng) for chng in changes]
            len_deltns = 0
            for c_pm,pm in enumerate(inds):
                len_deltns += changes[c_pm]-pm
            pre_RT_len = max([13, max(changes)-(pegPAM-3)+6+len_deltns])
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start()
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]
            RT_templ = [i for i in RT_templ]
            
            mmry = 0
            for c_pm,pm in enumerate(inds):
                del RT_templ[pm-(pegPAM-3)-mmry:changes[c_pm]-(pegPAM-3)-mmry]
                mmry += changes[c_pm]-pm
                
        RT_templ = ''.join(RT_templ)
        
        RT_templ = reverse_complement(RT_templ)        # RT-template must be in opposite direction
        
        pegRNA = pegspacer + scaffld.upper() + RT_templ + PBS
        
        out.append(['pegRNA'+str(c), 'pegRNA', pegRNA, pegPAM_strand])

        ## Define gRNA for PE3
        if mode == 'PE3':
            
            nuseq = [i for i in seq]
            mmry = 0

            for c2,ind in enumerate(inds):
                if edt[2] == 'mut':
                    nuseq[ind] = changes[c2]
                elif edt[2] == 'ins':
                    nuseq = nuseq[:ind] + [b for b in changes[c2]] + nuseq[ind:]
                elif edt[2] == 'del':
                    del nuseq[ind:changes[c2]]
            nuseq = ''.join(nuseq)
            
            gPAMrgn = reverse_complement(nuseq[pegPAM-3:])     # gRNA must bind to other strand somewhere downstream of pegPAM
            gPAMs = re.finditer(r'(?=(.GG))', gPAMrgn)     # find all PAMs in that region
            gPAMs = np.array([i.start() for i in gPAMs])
            gPAM = gPAMs[abs(np.array(gPAMs)-(len(gPAMrgn)-44)).argmin()] # Find PAM closest to 47 nt downstream of pegPAM nick
                                                                          # (the gRNA nick should be 50 nt downstream of pegPAM nick)

            gspacer = gPAMrgn[gPAM-20:gPAM]                # Define spacer

            gRNA = gspacer

            out.append(['gRNA'+str(c), 'gRNA', gRNA, pegPAM_strand])
    
    return out,msg


def PTGbldr(name, inserts, poltype='ptg'):
    '''
    Takes all desired parts of PTG GG assembly and gives out the respective inserts for PTG. During the 
    process, each part is appended with the same tRNA. The unit of part and tRNA are then treated as one insert.
    
    :param inserts: List of parts of the form [['type', 'sequence'],[...]]
    :type inserts: list
    :param poltype: Type of polycistronic architecture to use. Must be one of 'ptg' or 'ca', defaults to 'ptg'
    :type poltype: str, optional
    '''
    
    if poltype=='ptg':
        PTG_parts = []
        PTG_parts.append(Part(name+'_0', 'tRNA', tRNA))
        c = 1
    
        ## Take each coding sequence of the output and append it with a tRNA to the list
        for prt in inserts:
            if prt[0] not in ['pegRNA', 'gRNA', 'smRNA']:
                raise InvalidUsage("crRNAs should be processed in CA", status_code=400, payload={'pge': 'sequence.html', 'box': 'poltype_input'})
            if prt[0] == 'pegRNA':
                PTG_parts.append(Part(name+'_'+str(c), 'pegRNA', str(prt[1])))
                PTG_parts.append(Part(name+'_'+str(c+1), 'tRNA', tRNA))
                c += 2
            elif prt[0] == 'gRNA':
                PTG_parts.append(Part(name+'_'+str(c), 'gRNA', str(prt[1]) + scaffld + tRNA))
                c += 1
            elif prt[0] == 'smRNA':
                PTG_parts.append(Part(name+'_'+str(c), 'smRNA', str(prt[1]) + tRNA))
                c += 1
    
        return PTG_parts
        
    elif poltype=='ca':
        CA_parts = []
        c = 0
        
        for c,prt in enumerate(inserts):
            if prt[0] != 'crRNA':
                raise InvalidUsage("CA can only process crRNAs", status_code=400, payload={'pge': 'sequence.html', 'box': 'poltype_input'})
            CA_parts.append(Part(name+'_'+str(c), prt[0], DR + str(prt[1])))
            c += 1
        
        CA_parts.append(Part(name+'_'+str(c), 'DR', DR))
            
        return CA_parts
