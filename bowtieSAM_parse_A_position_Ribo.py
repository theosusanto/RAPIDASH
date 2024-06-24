#!/usr/bin/env python

"""
What it does: Parse bowtie output .SAM file
Input: bowtie output (.sam) file
Output:
gene ID, sequence, 5'UTR/CDS/3'UTR, position of A site relative to start site, position of A site relative to stop
"""


from sys import argv
script, bowtieSAMfile, UTR_CDS_length, output = argv


in_handle1 = open(bowtieSAMfile, "rU")
in_handle2 = open(UTR_CDS_length, "rU")
out_handle = open(output, "w")

count_total=0
count_mapped=0
count_UTR5=0
count_CDS=0
count_UTR3=0
count_multiple=0

dictionary={}
for entry in in_handle2:
    features = entry.split('\t')
    [name, UTR5, CDS, UTR3] = [features[0], features[1], features[2], features[3]]
    dictionary[name] = [int(UTR5), int(CDS), int(UTR3)]



for entry in in_handle1:
    if entry.find("@") == 0:     
        continue     #remove header lines, which starts with '@'
    elif entry.find("XS:i:") >= 0:
        count_multiple+=1
        continue     #remove reads with multiple alignments
    else:
        an_alignment = entry.rstrip('\n')      #load one alignment at a time
        count_total+=1
        alignmentFeature = an_alignment.split('\t')
    #alignmentFeature is a list with each column of bowtie output as one item
        title = alignmentFeature[0]         #column 1 is the title (without the '@' in front)
        sum_of_flags = alignmentFeature[1]   #column 2 is sum of all applicable flags
        gene = alignmentFeature[2].split("|")[0]     #column 3 is the name of gene it mapped to
        start = int(alignmentFeature[3]) 	     #column 4 is 1-based offset into the forward strand (leftmost)
        mapping_quality = alignmentFeature[4]
        sequence = alignmentFeature[9] 		 #column 10 is the read seq (reverse-complemented if aligned to - stand)
        Phred = alignmentFeature[10]        #column 11 is the Phred quality socre of the read
        '''
        now begin to parse the bowtie results
        three scenarios: 
        sum_of_flags = 4, meaning no alignment
        sum_of_flags = 0, meaning align to the + strand
        sum_of_flags = 16, meaning align to the - strand
         '''
        read_length = len(sequence)
        if read_length<=25:
            A_position = start + 4 + 9
        elif read_length==26:
            A_position = start + 4 + 10
        elif read_length==27:
            A_position = start + 4 + 10
        elif read_length==28:
            A_position = start + 4 + 2
        elif read_length==29:
            A_position = start + 4 + 10
        elif read_length==30:
            A_position = start + 4 + 2
        elif read_length==31:
            A_position = start + 4 + 2
        elif read_length==32:
            A_position = start + 4 + 11
        elif read_length==33:
            A_position = start + 4 + 11
        elif read_length>=34:
            A_position = start + 4 + 3

        if sum_of_flags=="0":
            count_mapped+=1
            [UTR5,CDS,UTR3] = dictionary[gene]
            if A_position <= UTR5:
                position="5UTR"
                count_UTR5+=1
                to_write1 = "%s\t%s\t%s\t%i\t%i\n" % (gene,sequence,position,A_position,A_position-UTR5)
                out_handle.write(to_write1)
            elif A_position >UTR5 and A_position<=(UTR5+CDS):
                position="CDS"
                count_CDS+=1
                to_write2 = "%s\t%s\t%s\t%i\t%i\n" % (gene,sequence,position,A_position-UTR5,A_position-UTR5-CDS)
                out_handle.write(to_write2)
            elif A_position>(UTR5+CDS) and start<=(UTR5+CDS+UTR3+17):
                position="3UTR"
                count_UTR3+=1
                to_write3 = "%s\t%s\t%s\t%i\t%i\n" % (gene,sequence,position,A_position-UTR5-CDS,A_position-UTR5-CDS-UTR3)
                out_handle.write(to_write3)
            else:
                print(gene, A_position, sequence)
            #out_handle.write(to_write)
            
           

in_handle1.close()
in_handle2.close()
out_handle.close()


print(f"Bowtie parser results: There are {count_total} total alignments, and {count_mapped} mapped reads. {count_UTR5} in 5'UTR, {count_CDS} in CDS and {count_UTR3} in 3'UTR. {count_multiple} reads were multiply aligned\n")

