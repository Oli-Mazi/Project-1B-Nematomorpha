# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 09:44:01 2025

@author: okmaz
"""
#####first pathway to files containing orthogroups are processed into a list
import os
files=sorted(os.listdir("C:/Users/okmaz/OneDrive/Desktop/Project 1B/Nematomorpha project/MSAseq"))
print(files)

####Aformentioned orthogroup files are turned into a dictionary, where each header is paired with its sequence
seqdict={}
for f in files: 
    with open("C:/Users/okmaz/OneDrive/Desktop/Project 1B/Nematomorpha project/MSAseq/"+ f, "r") as file:
        for line in file:
            line=line.strip()
            if line[0]==">":
                head=line
                seqdict[head]=""
            else: 
                seqdict[head]+=line  
   
                
####Next each orthogroup file has its sequences extracted and turned into a string
####where each string is made up of only one species. This is done via using aformentioned dictionarys key content
CaEl=""
NeMu=""
AcAu=""
PaVa=""
ChFu=""
ChFo=""
GoAq=""
GoMo=""
GoPa=""
 
              
keylist=list(seqdict.keys())

for key in keylist:
    if key[-4:]=="NeMu": 
        NeMu+=seqdict[key]
    if key[-4:]=="AcAu": 
        AcAu+=seqdict[key] 
    if key[-4:]=="PaVa": 
        PaVa+=seqdict[key]
    if key[-4:]=="ChFu":
        ChFu+=seqdict[key] 
    if key[-4:]=="ChFo": 
        ChFo+=seqdict[key] 
    if key[-4:]=="GoAq": 
        GoAq+=seqdict[key] 
    if key[-4:]=="GoMo": 
        GoMo+=seqdict[key]
    if key[-4:]=="GoPa": 
        GoPa+=seqdict[key] 
    if key[-2:]==".1":
        CaEl+=seqdict[key]


####Now the strings are turned into a fasta format style
####by first turning the strings into a dictionary with the key being the fasts file header

seqftreedict={">N.munidae": NeMu, ">A.australiensis":AcAu,">P.varius":PaVa, ">C.fukuii":ChFu, ">C.formosanus": ChFo, ">G.aquaticus": GoAq, ">G.montsenyensis": GoMo, ">G.paranensis": GoPa, ">C.elegans": CaEl}

with open("C:/Users/okmaz/OneDrive/Desktop/Project 1B/Nematomorpha project/MSAseq/sumMSA.fa", "w") as file:
    for head, seq in seqftreedict.items():
        file.write(head+"\n")
        file.write(seq+"\n")

####Next partion files are made by using the primary dictionary containing the seperated sequences
####the length of each sequence can be calculated
####Thus the concatenated gene start and end positions within in the concatenated sequence to be calculaed

sizelist=[]
for i in keylist: 
    sizelist.append(len(seqdict[i]))
rsizelist=[]
for i in range(1,19):
    p=(i*9)-1
    rsizelist.append(sizelist[p])

endpos=[]
for l in range(0,19):
    endpos.append(sum(rsizelist[0:l]))

endpos=endpos[1:]

startpos=[1]
for l in range(0,18):
    sp=endpos[l]+1
    startpos.append(sp)

startpos=startpos[0:18]

####The lists with start and end positions can be combined into one file, in the Rax-ML partion file style

with open("C:/Users/okmaz/OneDrive/Desktop/Project 1B/Nematomorpha project/MSAseq/posMSA.partition", "w") as file:
    for n in range(0,18):
        line="AA, part"+str((n+1))+"="+str(startpos[n])+"-"+str(endpos[n])+"\n"
        file.write(line)
