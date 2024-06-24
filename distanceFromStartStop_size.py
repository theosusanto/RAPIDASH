from sys import argv
script, mapped_file, output = argv

in_handle = open(mapped_file, "rU")
out_handle = open(output, "w")

#Initialize the dictionary
UTR5_dic={}
CDS_stop={}
for i in range(25,39):
    UTR5_dic[i]={}
    CDS_stop[i]={}
    for j in range(-27,1):
        UTR5_dic[i][j] = 0
        CDS_stop[i][j] = 0


for entry in in_handle:
    entry.rstrip('\n')
    features = entry.split('\t')
    [gene,sequence,position,relative_start,relative_stop] = [features[0], features[1],features[2],int(features[3]),int(features[4])]
    length = len(sequence)  #5'first nuc is addition by RT however the assignment of A site also have this 
    if length>=25 and length<=38:
        if position == '5UTR':
            if relative_stop <=0 and relative_stop >= -27:
                UTR5_dic[length][relative_stop] += 1
        elif position == 'CDS':
            if relative_start>100 and relative_stop>= -27 and relative_stop <=0:
                CDS_stop[length][relative_stop] += 1
                        
            
for m in range(25,39):
    for n in range(-27,1):
        to_write = "%i\t%i\tUTR\t%i\tCDS\t%i\n" % (m,n,UTR5_dic[m][n], CDS_stop[m][n])
        out_handle.write(to_write)

    
in_handle.close()
out_handle.close()