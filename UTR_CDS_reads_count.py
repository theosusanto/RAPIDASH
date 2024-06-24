'''
count the number of reads mapped to 5'UTR, first 15 codons of CDS, CDS excluding first 15 and last 5 codons, last 5 codons on CDS as well as 3'UTR
'''

from sys import argv
import numpy
import pickle

script, inputfile, outputfile = argv

# Initialize the dictionary
# The dictionary will have four numbers: reads in the 5'UTR, the first 15 amino acid of CDS, rest of CDS and then 3'UTR
reads_count={}  


with open(inputfile, 'r') as in_handle:
    for entry in in_handle:
        entry.rstrip('\n')
        features = entry.split('\t')
        [gene,sequence,position,relative_start,relative_stop] = [features[0], features[1],features[2],int(features[3]),int(features[4])]
        read_length=len(sequence)
        print(f'gene: <{gene}>, position: <{position}>')
        #first assign the A site based on the formula
        
        if position == '5UTR':
            
            print('found 5UTR')
            
            if not gene in reads_count:
                reads_count[gene]=numpy.zeros(5)
                reads_count[gene][0] +=1
            elif gene in reads_count:
                reads_count[gene][0] +=1
        
        elif position =='CDS' and relative_start <=45: # first 15 amino acid, see how ribosome pile up may correlate with TE or anything
            
            print('found CDS, start before 45')
            
            if not gene in reads_count:
                reads_count[gene]=numpy.zeros(5)
                reads_count[gene][1] +=1
            elif gene in reads_count:
                reads_count[gene][1] +=1
        
        elif position =='CDS' and relative_start >45 and relative_stop<=-15: # middle of CDS excluding first 15 codons and last 5 codons
            
            print('found CDS, neither first nor last few bases')
            
            if not gene in reads_count:
                reads_count[gene]=numpy.zeros(5)
                reads_count[gene][2] +=1
            elif gene in reads_count:
                reads_count[gene][2] +=1
        
        elif position =='CDS' and relative_stop > -15: # CDS last 5 codons
            
            print('found CDS, near end')
            
            if not gene in reads_count:
                reads_count[gene]=numpy.zeros(5)
                reads_count[gene][3] +=1
            elif gene in reads_count:
                reads_count[gene][3] +=1     
        
        elif position =='3UTR':
            
            print('found 3UTR')
            
            if not gene in reads_count:
                reads_count[gene]=numpy.zeros(5)
                reads_count[gene][4] +=1
            elif gene in reads_count:
                reads_count[gene][4] +=1 
        
        else:
            print('could not categorize position')

print(f'length of dictionary: {len(reads_count)}')

with open(outputfile, 'w') as out_handle:
    for gene in reads_count.keys():
        print(f'writing {gene} and {reads_count[gene]} to file')
        to_write=f"{gene}\t{reads_count[gene][0]}\t{reads_count[gene][1]}\t{reads_count[gene][2]}\t{reads_count[gene][3]}\t{reads_count[gene][4]}\n"
        out_handle.write(to_write)
# gene_name, 5'UTR, first 15 codon, main CDS, last 5 codon, 3'UTR

out_handle.close()
