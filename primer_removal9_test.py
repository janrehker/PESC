import subprocess, sys, os, gzip, io, time, itertools, re, math
from sets import Set


#input1 = sys.argv[1]  # (realigned) mate sorted bam
input2 = sys.argv[1]  # primer bed (chr\tpos\t0/1(+- strand)\tprimer sequence(5'->3')
#output1 = sys.argv[3]

# output goes to stdout!! pipe into samtools, e.g. python primer_removal2.py *.bam primer.bed | samtools view -bh -o output.bam

#annotates size of the read group into read name in bam


alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases


# check for paired, output paircheck: 1 first in pair, 2 second in pair, 3 both (should not happen)
pair_check = 0
def check_pair(flag):
    if flag & 64:
        return '1'
    if flag & 128:
        return '2'
    if (flag & 128) and (flag & 64):
        return '3'

strand = ''
def check_strand(flag):
    if flag & 16:
        return 'r_reverse' ## read reverse strand
    if flag & 32:
        return 'm_reverse' ## mate reverse strand

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def roundup(x):
    return int(math.ceil(x / 1000.0)) * 1000



prim_list = {}
    

with open(input2, 'r') as prim_file:
    for line in prim_file:
        if line != '':
            linesplit = line.strip('\n').split('\t')
            chroms = linesplit[0]
            POS = roundup(int(linesplit[1]))
            if str(chroms + '_' + str(POS)) in prim_list:
                prim_list[str(chroms + '_' + str(POS))].append(linesplit[1:])
            if str(chroms + '_' + str(POS)) not in prim_list:
                prim_list[str(chroms + '_' + str(POS))]=[linesplit[1:]]

# instead of adding lists, better convert each list into dictionary of pos as key, list of strand, sequence as value 

#print prim_list['chr5_138150000']


name = ''
flag = ''
chrom = ''
pos = ''
CIGAR = ''
sequence = ''


list_out = []

#with open('error.txt', 'a') as error_out:
#    output = subprocess.Popen(['samtools', 'view', '-h', input1], stdout=subprocess.PIPE, stderr=error_out).communicate()[0]
#split = output.split('\n')

count = 0

for i in sys.stdin:
    remove = 0
    count = count  + 1
    deconstruct = i.split('\t')
    list_M = 1
    if i != '' and i[0:5] == 'M0072':  # lib1 is library name
        name = deconstruct[0]
        flag = int(deconstruct[1])
        chrom = deconstruct[2]
        pos = int(deconstruct[3])
        CIGAR = deconstruct[5]
        RNEXT = deconstruct[6] # mate name (should be = )
        mate = int(deconstruct[7]) # mate position
        frag_len = deconstruct[8] # length of fragment
        frag_len = (abs(int(frag_len)))
        sequence = deconstruct[9]
        sequence_length = len(sequence)
        MD = deconstruct[12]
        CIGAR = deconstruct[5]
        reCIGAR_c = (re.split('[0-9]+', CIGAR, flags=re.IGNORECASE))[1:]  ### grep letters
        reCIGAR_v = re.split('[a-z]+', CIGAR, flags=re.IGNORECASE)[0:-1]  ### grep values
        list_CIG = listz = [list(ig) for ig in zip(reCIGAR_v, reCIGAR_c)]
        countMs = 0
        count_S1 = 0
        count_S2 = 0
        for l in list_CIG:
            if 'S' in l and countMs == 0:
                count_S1 = int(l[0])
            if 'S' in l and countMs > 0:
                count_S2 = int(l[0])
            if 'M' in l:
                countMs = countMs + int(l[0])
            if 'D' in l:
                countMs = countMs + int(l[0])
            if 'I' in l:
                countMs = countMs - int(l[0])
                
        POS_1 = roundup(int(pos - 150))
        POS_2 = roundup(int(pos + 150))
        MATE_1 = roundup(int(mate - 150))
        MATE_2 = roundup(int(mate + 150))
        remove = 0
        drawer1_1 = chrom + '_' + str(POS_1)
        drawer1_2 = chrom + '_' + str(POS_2)
        drawer2_1 = chrom + '_' + str(MATE_1)
        drawer2_2 = chrom + '_' + str(MATE_2)
        chrom_list2 = []
        #print 'hello--------------------'
        #print drawer1
        #print pos
        if drawer1_1 in prim_list or drawer1_2 in prim_list or drawer2_1 in prim_list or drawer2_2 in prim_list:
            #print 'its in!'
            chrom_list1_1 = []
            chrom_list1_2 = []
            chrom_list2_1 = []
            chrom_list2_2 = []
            if drawer1_1 in prim_list:
                chrom_list1_1 = prim_list[drawer1_1]
            if drawer1_2 in prim_list:
                chrom_list1_2 = prim_list[drawer1_2]
            if drawer2_1 in prim_list:
                chrom_list2_1 = prim_list[drawer2_1]
            if drawer2_2 in prim_list:
                chrom_list2_2 = prim_list[drawer2_2]
            chrom_list = chrom_list1_1 + chrom_list1_2 + chrom_list2_1 + chrom_list2_2
            chrom_list.sort()
            chrom_list1 = list(chrom_list for chrom_list,_ in itertools.groupby(chrom_list))  ### remove overlaps
            #if name == 'lib1:4796210':
            #print chrom_list1
            
            for l in chrom_list1:
                if 150 >= int(l[0]) - pos >= -150 or 150 >= int(l[0]) - mate >= -150: # narrow it down. Primer has to be either close to read position in read1 or close to the mate position if we observe read2
                    #print l
                    chrom_list2.append(l)
        
            
        #print chrom_list2
        list_M = 1
        cut = 0
        sub = 0
        #print 'first'
        for pr in chrom_list2:
            if 'H' not in CIGAR and CIGAR != '*' and cut == 0 and pr != '':
                prim_split = pr
                if (check_pair(flag)) == '1':  #read1
                    #print 'first'
                    if 1 >= pos - (int(prim_split[0]) - len(prim_split[2].strip())) >= -1 and check_strand(flag) == 'm_reverse' and prim_split[1] == '0':
                        #print 'first found'
                        #print pr
                        #print CIGAR, pos
                        length_cut_read1 = len(prim_split[2].strip())  ## length of sequence that needs to be removed from left M and added to left S
                        list_left_S = ''
                        list_M = 0
                        cut = 0
                        seen_M = 0
                        seen_left_S = 0
                        for l in list_CIG:
                            if 'S' in l and seen_M == 0: # must be the left S
                                list_left_S = str(int(l[0]) + length_cut_read1) # this would be the corrected S if left M is spanning at least the length of length_cut_read2
                                seen_left_S = 1
                            if 'M' in l and seen_M == 0:  # left M, we want this one!
                                seen_M = 1
                                if int(l[0]) >= length_cut_read1: # M is long enough, we can cut it!
                                    list_M = int(l[0]) - length_cut_read1
                                    if list_M >= 0:
                                        cut = 1 # cut it, please
                                    else:
                                        remove = 1   # better keep that?
                        seen_M = 0
                        if cut == 1: # let's do it
                            sub = 1
                            for t in list_CIG:
                                if 'S' in t and seen_M == 0:
                                    t[0] = str(list_left_S)
                                if 'M' in t and seen_M == 0:  
                                    seen_M = 1
                                    t[0] = str(list_M)
                                    
                        CIGAR = ''.join(itertools.chain(*list_CIG))
                        if cut == 1 and seen_left_S == 0:
                            CIGAR = str(length_cut_read1) + 'S' + CIGAR
                        if cut == 1:
                            pos = pos + length_cut_read1
                        #print CIGAR, pos
                        #if name == 'M00727:24:000000000-CD2KV:1:1111:10111:4159':
                        #    print pos, length_cut_read1, str(prim_split[2].strip())


                    if sequence[0 + count_S1:(len(prim_split[2].strip()) + count_S1)] != prim_split[2].strip() and countMs <= len(prim_split[2].strip('\n')): #read equal/smaller than primer, remove read
                        #print 'remove read ' + CIGAR + ' ' + str(len(prim_split[3].strip('\n')))
                        remove = 1
                        
                    if 1 >= (int(mate) + frag_len) - (int(prim_split[0]) + len(prim_split[2].strip())) >= -1 and check_strand(flag) == 'r_reverse' and prim_split[1] == '1': #primer on the right end of the read: pos(alignment start) + Ms(alignment length) == prim-position + len(primer)
                        #print 'second found'
                        length_cut_read1 = mate + frag_len - (int(prim_split[0]))
                        list_right_S = ''
                        list_M = 0
                        cut = 0
                        seen_M = 0
                        seen_right_S = 0
                        for l in list_CIG:
                            if 'S' in l and seen_M > 0: # must be the right S
                                list_right_S = str(int(l[0]) + length_cut_read1) # this would be the corrected S if left M is spanning at least the length of length_cut_read2
                                seen_right_S = 1
                            if 'M' in l:
                                seen_M = seen_M + 1
                                if int(l[0]) >= length_cut_read1: # M is long enough, we can cut it!
                                    list_M = int(l[0]) - length_cut_read1
                                    if list_M >= 0:
                                        cut = 1 # cut it, please
                                    else:
                                        remove = 1
                                    seen_M2 = seen_M
                        if seen_M2 != seen_M: #last cut M was not the last M that appeared the last one was shorter than the part that needs cutting
                            cut =  0
                        seen_M = 0
                        if cut == 1: # let's do it
                            sub = 2
                            for t in list_CIG:
                                if 'M' in t:  
                                    seen_M = seen_M + 1
                                    if seen_M == seen_M2:  # is it the right M?
                                        t[0] = str(list_M)
                                if 'S' in t and seen_M == seen_M2:
                                    t[0] = str(list_right_S)
                        CIGAR = ''.join(itertools.chain(*list_CIG))
                        if cut == 1 and seen_right_S == 0:
                            CIGAR = CIGAR + str(length_cut_read1) + 'S'


                    if sequence[-(len(prim_split[2].strip()) - 1 - count_S2):] != reverse_complement(prim_split[2].strip()) and countMs <= len(prim_split[2].strip('\n')):
                        remove = 1

                if (check_pair(flag)) == '2':  #read2, second in pair
                    if RNEXT == '=' and 1 >= int(mate) - (int(prim_split[0]) - len(prim_split[2].strip())) >= -1 and check_strand(flag) == 'r_reverse' and prim_split[1] == '0':
                        #print 'third found'
                        if pos < (int(prim_split[0]) + 1):
                            length_cut_read2 = int(prim_split[0]) + 1 - pos  ## length of sequence that needs to be removed from left M and added to left S
                            list_left_S = ''
                            list_M = 0
                            cut = 0
                            seen_M = 0
                            seen_left_S = 0
                            for l in list_CIG:
                                if 'S' in l and seen_M == 0: # must be the left S
                                    list_left_S = str(int(l[0]) + length_cut_read2) # this would be the corrected S if left M is spanning at least the length of length_cut_read2
                                    seen_left_S = 1
                                if 'M' in l and seen_M == 0:  # left M, we want this one!
                                    seen_M = 1
                                    if int(l[0]) >= length_cut_read2: # M is long enough, we can cut it!
                                        list_M = int(l[0]) - length_cut_read2
                                        if list_M >= 0:
                                            cut = 1 # cut it, please
                                        else:
                                            remove = 1
                            seen_M = 0
                            if cut == 1: # let's do it
                                sub = 3
                                for t in list_CIG:
                                    if 'S' in t and seen_M == 0:
                                        t[0] = str(list_left_S)
                                    if 'M' in t and seen_M == 0:  
                                        seen_M = 1
                                        t[0] = str(list_M)
                            CIGAR = ''.join(itertools.chain(*list_CIG))
                            if cut == 1 and seen_left_S == 0:
                                CIGAR = str(length_cut_read2) + 'S' + CIGAR
                            if cut == 1:
                                pos = pos + length_cut_read2

                            ### remove read if alignment gets too small??? ####
                        

                    if RNEXT == '=' and 1 >= (pos + int(frag_len)) - (int(prim_split[0]) + 1 + len(prim_split[2].strip())) >= -1 and check_strand(flag) == 'm_reverse' and prim_split[1] == '1':
                        #print 'fourth found'
                    
                        if pos + countMs - 1 > (int(prim_split[0]) + 1): #overlap: correct CIGAR on the right side
                            length_cut_read2 = pos + countMs - (int(prim_split[0]))
                            list_right_S = ''
                            list_M = 0
                            cut = 0
                            seen_M = 0
                            seen_right_S = 0
                            for l in list_CIG:
                                if 'S' in l and seen_M > 0: # must be the right S
                                    list_right_S = str(int(l[0]) + length_cut_read2) # this would be the corrected S if left M is spanning at least the length of length_cut_read2
                                    seen_right_S = 1
                                if 'M' in l:
                                    seen_M = seen_M + 1
                                    if int(l[0]) >= length_cut_read2: # M is long enough, we can cut it!
                                        list_M = int(l[0]) - length_cut_read2
                                        if list_M >= 0:
                                            cut = 1 # cut it, please
                                        else:
                                            remove = 1
                                        seen_M2 = seen_M
                            if seen_M2 != seen_M: #last cut M was not the last M that appeared the last one was shorter than the part that needs cutting
                                cut =  0
                            seen_M = 0
                            if cut == 1: # let's do it
                                sub = 4
                                for t in list_CIG:
                                    if 'M' in t:  
                                        seen_M = seen_M + 1
                                        if seen_M == seen_M2:  # is it the right M?
                                            t[0] = str(list_M)
                                    if 'S' in t and seen_M == seen_M2:
                                        t[0] = str(list_right_S)
                            CIGAR = ''.join(itertools.chain(*list_CIG))
                            if cut == 1 and seen_right_S == 0:
                                CIGAR = CIGAR + str(length_cut_read2) + 'S'

        reCIGAR_c = (re.split('[0-9]+', CIGAR, flags=re.IGNORECASE))[1:]  ### grep letters
        reCIGAR_v = re.split('[a-z]+', CIGAR, flags=re.IGNORECASE)[0:-1]  ### grep values
        list_CIG = listz = [list(ig) for ig in zip(reCIGAR_v, reCIGAR_c)]  ## reconstruct cigar
        count_CIGAR = 0


        deconstruct[3] = str(pos)
        deconstruct[5] = CIGAR
    if list_M > 0:
        deconstruct = '\t'.join(deconstruct)
        #list_out.append(deconstruct)
    #if count == 200:
        #streamout = ''.join(list_out)
        sys.stdout.write(str(deconstruct))
        list_out = []
        count  = 0

#streamout = '\n'.join(list_out)
#sys.stdout.write(deconstruct)
#list_out = []


