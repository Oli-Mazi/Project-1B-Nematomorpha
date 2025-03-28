# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 09:39:34 2025

@author: okmaz
"""
####import and make a list of files containing the orthogroups (each file contains one orthogroup)
import os
files=os.listdir("C:/Users/okmaz/OneDrive/Desktop/Project 1B/Nematomorpha project/Tree-seq-file")
print(files)

####Below code goes through each file, adds the content line by line into a new file
#### but at the same time removes "*"at the end of predicted amino acid sequences
for n in files: 
    with open("C:/Users/okmaz/OneDrive/Desktop/Project 1B/Nematomorpha project/Tree-seq-file/"+ n, "r") as file:
        with open("C:/Users/okmaz/OneDrive/Desktop/Project 1B/Nematomorpha project/cleanseqfile/"+ "g"+ n , "w") as newfile:
            for line in file:
                line=line.strip("\n")
                if line[-1]=="*":
                    nline=line.strip("*")
                    newfile.write(nline+"\n")
                else:
                    nline=line
                    newfile.write(nline+"\n")