"""
    Python library for bioinformatics. It covers some of problems 
    listed in (http://rosalind.info/problems/list-view/).
    Copyright (C) 2017  Ali Mert Ceylan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import itertools
from Bio import SeqIO

def rna_codon_table(path=None):
    """
    Loads RNA Codon table from file.
    
    Returns a dictionary of RNA squences with corresponding protein strings.
    """
    rna_table = open(path or 'rna_codon_table.txt')
    rna_dict = {}
  
    rna_table_text = rna_table.readlines()

    for i in range(len(rna_table_text)):
        seq,prot = rna_table_text[i].split('      \n')[0].rstrip('\n').split(' ')
        rna_dict[seq] = prot

    return rna_dict


def _load_blosum62(path=None):
    """Loads blosum62 scoring matrix from file."""
    lines = open('blosum62.txt').readlines()
    blosum62 = {}
    
    for i in range(len(lines)):
        lines[i] = lines[i].rstrip('\n')
        lines[i] = lines[i].split()

    return lines


def blosum62(path=None):
    """
    Loads blosum62 scoring matrix as dictionaries.
    
    blosum62["something"]["somethingelse"]
    """
    table={} 
    lines = _load_blosum62(path)
    
    labels = lines[0]
    lines = lines[1:]
    for i in range(len(labels)):
        temp_dict_row = {}
        for j in range(len(labels)):
            temp_dict_row[labels[j]] = lines[i][j]
        table[labels[i]] = temp_dict_row
        
    return table


def ones(path=None):
    """Returns a scoring matrix of 1 only."""
    table={}
    lines = _load_blosum62(path)

    labels = lines[0]
    lines = lines[1:]
    for i in range(len(labels)):
        temp_dict_row = {}
        for j in range(len(labels)):
            temp_dict_row[labels[j]] = 1
        table[labels[i]] = temp_dict_row
        
    return table


def import_fasta(d):
    """Loads fasta file from directory."""

    records = []
    for record in SeqIO.parse(d, "fasta"):
        records.append(record)
    
    return records


def parse_fasta(r):
    """Parse fasta file to list."""

    sequences = []
    for i in range(len(r)):
        sequences.append(r[i].seq)
    
    return sequences


def transcribe_dna_to_rna(s):
    """
    http://rosalind.info/problems/rna/
    """
    transcribed_rna = ''
    for i in range(len(s)):
        bp = str(s[i])
        if(bp == 'T'):
            transcribed_rna+='U'
        else:
            transcribed_rna+=s[i]

    return transcribed_rna


def reverse_complement_dna(s):
    """
    http://rosalind.info/problems/revc/
    """
    reverse_complement = ''
    for i in range(len(s)):
        bp = str(s[len(s)-i-1])
        if(bp == 'A'):
            reverse_complement+='T'
        elif(bp == 'T'):
            reverse_complement+='A'
        elif(bp == 'G'):
            reverse_complement+='C'
        elif(bp == 'C'):
            reverse_complement+='G'
        else:
            print('Something went wrong.')
            break

    return reverse_complement


def point_mutations(seq1, seq2):
    """
    http://rosalind.info/problems/hamm/
    """
    if(len(seq1) != len(seq2)):
        print('Size mismatch')
        #TODO throw somekind of exception here.

    pms = 0
    for i in range(len(seq1)):
        if(seq1[i] != seq2[i]):
            pms+=1

    return pms


def translate_rna(s):
    """
    http://rosalind.info/problems/prot/    
    """
    rna_dict = rna_codon_table()
    _prot_list = []
    for i in range(0,len(s)-3,3):
        _prot_list.append(rna_dict[s[i:i+3]])

    return _prot_list


def find_motif(s, m):
    """
    http://rosalind.info/problems/subs/
    """
    sl = len(m)
    match_points = []
    for i in range(len(s)-sl):
        if(m in s[i:(i+sl)]):
            match_points.append(str(i+1))  
    
    return match_points


def _lcsm(seq1, seq2):
    """
    http://rosalind.info/problems/lcsq/
    """
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    
    gridlen1 = seq1_len+1
    gridlen2 = seq2_len+1
    
    dist_grid = [[0 for x in range(gridlen2)] for y in range(gridlen1)]
    trace_grid = [['--' for x in range(gridlen2)] for y in range(gridlen1)]
    
    for i in range(gridlen1):
        dist_grid[i][0] = 0

    for i in range(gridlen2):
        dist_grid[0][i] = 0
        
    for i in range(1,gridlen1):
        for j in range(1,gridlen2):
            if(seq1[i-1]==seq2[j-1]):
                dist_grid[i][j] = dist_grid[i-1][j-1]+1
                trace_grid[i][j]='↖'#cross back
            else:
                dist_grid[i][j] = max(dist_grid[i-1][j], dist_grid[i][j-1])
                if(dist_grid[i][j]==dist_grid[i-1][j]):
                    trace_grid[i][j]='↑ '#up back
                elif(dist_grid[i][j]==dist_grid[i][j-1]):
                    trace_grid[i][j]='←'#left back

    return [dist_grid, trace_grid]


def _lcsm_glob(seq1, seq2, penalty_dict={'sigma': None, 'score':None}):
    """
    http://rosalind.info/problems/lcsm/
    http://rosalind.info/problems/glob/
    """

    seq1_len = len(seq1)
    seq2_len = len(seq2)
    
    gridlen1 = seq1_len+1
    gridlen2 = seq2_len+1

    sigma =  penalty_dict['sigma'] or 0
    score = penalty_dict['score'] or ones()

    dist_grid = [[0 for x in range(gridlen2)] for y in range(gridlen1)]
    trace_grid = [['--' for x in range(gridlen2)] for y in range(gridlen1)]
    
    for i in range(gridlen1):
        dist_grid[i][0] = -1*i*sigma
        # dist_grid[i][0] = i*sigma

    for i in range(gridlen2):
        dist_grid[0][i] = -1*i*sigma
        # dist_grid[0][i] = i*sigma
        
    for i in range(1,gridlen1):
        for j in range(1,gridlen2):
            maxl = [dist_grid[i-1][j]-sigma, dist_grid[i][j-1]-sigma]
            # if(seq1[i-1]==seq2[j-1]):
                # maxl.append(dist_grid[i-1][j-1]+int(score[seq1[i-1]][seq2[j-1]]))
            # else:
            maxl.append(dist_grid[i-1][j-1]+int(score[seq1[i-1]][seq2[j-1]]))

            dist_grid[i][j] = max(maxl)#max(dist_grid[i-1][j]-sigma, dist_grid[i][j-1]-sigma)
            if(dist_grid[i][j]==dist_grid[i-1][j]-sigma):
                trace_grid[i][j]='↑ '#up back
            elif(dist_grid[i][j]==dist_grid[i][j-1]-sigma):
                trace_grid[i][j]='←' #left back
            elif(dist_grid[i][j]==dist_grid[i-1][j-1]+int(score[seq1[i-1]][seq2[j-1]])):
                trace_grid[i][j]='↖' #cross back
            else:
                print('An error occured.')
            
    return [dist_grid, trace_grid]


def extract_lcsm_seqs(seq1, seq2, trace_grid):
    """For alignment"""
    lcsm_seq1=''
    lcsm_seq2=''
    i = len(seq1)
    j = len(seq2)
    while(i>=0 and j>=0):
        if(trace_grid[i][j]=='↑ '):
            lcsm_seq1+=seq1[i-1]
            lcsm_seq2+='-'
            i-=1
        elif(trace_grid[i][j]=='←'):
            lcsm_seq1+='-'
            lcsm_seq2+=seq2[j-1]
            j-=1
        else:
            if(trace_grid[i][j]=='--'):
                if(i>0 and j==0):
                    lcsm_seq1+=seq1[i-1]
                    lcsm_seq2+='-'
                elif(j>0 and i==0):
                    lcsm_seq1+='-'
                    lcsm_seq2+=seq2[j-1]
            else:
                lcsm_seq1+=seq1[i-1]
                lcsm_seq2+=seq2[j-1]
            i-=1
            j-=1
    return [lcsm_seq1[::-1], lcsm_seq2[::-1]]

def extract_lcsm_seq(seq1, seq2, trace_grid):
    """http://rosalind.info/problems/lcsm/"""
    lcsm_seq=''
    i = len(seq1)
    j = len(seq2)
    
    while(i > 0 and j > 0):
        if(trace_grid[i][j]=='↑ '):
            i-=1
        elif(trace_grid[i][j]=='←'):
            j-=1
        else:
            if(trace_grid[i][j]=='--'):
                break
            lcsm_seq+=seq1[i-1]
            i-=1
            j-=1
    return [lcsm_seq[::-1]]

def _extract_lcsm_sub_string(seq1, seq2, trace_grid, min=None):
    """http://rosalind.info/problems/lcsm/"""
    lcsm_list = []
    lcsm_seq=''
    i = len(seq1)
    j = len(seq2)
    
    for a in range(len(seq1),0,-1):
        for b in range(len(seq2),0,-1):
            i = a*1
            j = b*1
            lcsm_seq=''
            
            while(trace_grid[i][j]=='↖'):
                lcsm_seq+=seq1[i-1]
                i-=1
                j-=1
            
            if(len(lcsm_seq) > (min or 0)):
                lcsm_list.append(lcsm_seq[::-1])
                lcsm_seq=''
    return lcsm_list


def _get_subseq_grid(_my_records):
    """
    http://rosalind.info/problems/lcsm/
    http://rosalind.info/problems/lcsq/
    """

    subseq_grid = [[[] for x in range(len(_my_records))] for y in range(len(_my_records))]
    
    for i in range(len(_my_records)):
        for j in range(len(_my_records)):
            if(i!=j):
                dist_grid, trace_grid = _lcsm(_my_records[i],_my_records[j])
                # subseq_grid[i][j]= _extract_lcsm_sub_string(_my_records[i],_my_records[j],trace_grid)
                subseq_grid[i][j]= extract_lcsm_seq(_my_records[i], _my_records[j], trace_grid)
    return subseq_grid


def _get_subset_dict(_my_records, _subseq_grid):
    """
    http://rosalind.info/problems/lcsm/
    http://rosalind.info/problems/lcsq/
    """

    _subseq_dict_grid = [[{} for x in range(len(_my_records))] for y in range(len(_my_records))]
    subseq_dict = {}
    
    for i in range(len(_my_records)):
        for j in range(len(_my_records)):
            for k in _subseq_grid[i][j]:#it wasnt here
                if(k in _subseq_dict_grid[i][j]):
                    _subseq_dict_grid[i][j][k]+=1
                else:
                    _subseq_dict_grid[i][j][k]=1

    for i in range(len(_my_records)):
        for j in range(len(_my_records)):
            if(_subseq_dict_grid[i][j]):
                for k in _subseq_dict_grid[i][j].keys():
                    if(k in subseq_dict):
                        subseq_dict[k]+=1
                    else:
                        subseq_dict[k]=1
                        
    return subseq_dict

def longest_common_substring(_my_records):
    """http://rosalind.info/problems/lcsm/"""
    
    _lcss = ''
    
    dist_grid, trace_grid = _lcsm(_my_records[0],_my_records[1])
    _some_sub_sequences = _extract_lcsm_sub_string(_my_records[0],_my_records[1],trace_grid)
    _some_sub_sequences.sort(key = len)
                
    for s in _some_sub_sequences:
        cnt = 0
        for r in _my_records:
            if(s in r and len(s)>len(_lcss)):
                cnt+=1
                _lcss=s
            if(cnt == len(_my_records)):
                print('Break')
                return _lcss

    return _lcss        

def perm(s):
    """http://rosalind.info/problems/perm/"""
    
    l = range(1, len(s)+1)
    v = itertools.permutations(l)

    pl = []
    perms = []
    cnt=0
    for item in v:
        sn=''
        for i in item:
            sn+=s[int(i-1)]
        pl.append(sn)
        perms.append(item)
        cnt+=1

    return cnt, perms, pl

def locate_restriction_sites(_length_min, _length_max, _dna_string):
    """http://rosalind.info/problems/revp/"""
    
    rs = []
    reverse_complement_dna_string = reverse_complement_dna(_dna_string)
    for i in range(len(_dna_string)):
        for length in range(_length_min, _length_max):        
            if(length%2==0):
                _length = int(length/2)
            
                if(i >= len(_dna_string)-(_length)):
                    break
            
                if(reverse_complement_dna(_dna_string[i:i+_length]) == _dna_string[i+_length:i+2*_length]):
                    rs.append('%d, %d'%(i+1, length))
    
    return rs


def remove_introns(original_dna_string, _introns):
    """http://rosalind.info/problems/splc/"""
    _dna_string = original_dna_string[:]

    for i in range(len(_introns)):
        idx = _dna_string.find(_introns[i])
        while(idx != -1):
            _dna_string = _dna_string[0:idx] + _dna_string[idx+len(_introns[i]):(len(_dna_string))]
            idx=-1

    return _dna_string


def find_spliced_motif(s, t):
    """http://rosalind.info/problems/sseq/"""
    index_list = [0]

    for c in t:
        i = s.find(c, index_list[len(index_list)-1]+2, len(s))
        if(i > -1):
            index_list.append(i)

    return index_list[1:]


def find_shared_spliced_motif(ss):
    """http://rosalind.info/problems/lcsq/"""
    subseq_grid = _get_subseq_grid(ss)
    subseq_dict = _get_subset_dict(ss, subseq_grid)

    return subseq_dict.keys()
    

def extract_lcsm_seq(seq1, seq2, trace_grid):
    lcsm_seq=''
    i = len(seq1)
    j = len(seq2)
    while(i>0 and j>0):
            if(trace_grid[i][j]=='↑ '):
                i-=1
            elif(trace_grid[i][j]=='←'):
                j-=1
            else:
                if(trace_grid[i][j]=='--'):
                    break
                lcsm_seq+=seq1[i-1]
                i-=1
                j-=1
    return [lcsm_seq[::-1]]


def edit_distance(seq1, seq2):
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    
    gridlen1 = seq1_len+1
    gridlen2 = seq2_len+1
    
    dist_grid = [[0 for x in range(gridlen2)] for y in range(gridlen1)]
    trace_grid = [['--' for x in range(gridlen2)] for y in range(gridlen1)]
    
    for i in range(gridlen1):
        dist_grid[i][0] = i

    for i in range(gridlen2):
        dist_grid[0][i] = i
        
    for i in range(1,gridlen1):
        for j in range(1,gridlen2):
            cl = [dist_grid[i-1][j]+1, dist_grid[i][j-1]+1]
            if(seq1[i-1]==seq2[j-1]):
                dist_grid[i][j] = dist_grid[i-1][j-1]
                trace_grid[i][j]='↖'#cross back
            else:
                cl.append(dist_grid[i-1][j-1]+1)
                dist_grid[i][j]=min(cl)
                
                if(dist_grid[i][j]==dist_grid[i-1][j]+1):
                    trace_grid[i][j]='↑ '#up back
                elif(dist_grid[i][j]==dist_grid[i][j-1]+1):
                    trace_grid[i][j]='←'#left back
                else:
                    trace_grid[i][j]='↖'#cross back

    return [dist_grid[-1][-1], dist_grid, trace_grid]
