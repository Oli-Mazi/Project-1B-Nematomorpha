######Data selection
NeMu=read.table(file.choose(), sep="\t")
ChFu=read.table(file.choose(), sep="\t")
ChFo=read.table(file.choose(), sep="\t")
PaVa=read.table(file.choose(), sep="\t")
GoMo=read.table(file.choose(), sep="\t")
GoAq=read.table(file.choose(), sep="\t")
GoPa=read.table(file.choose(), sep="\t")
AcAu=read.table(file.choose(), sep="\t")

#####Finding all types of labels for colour coordination

looking1=data.frame(unique(NeMu$V3))
looking2=data.frame(unique(ChFu$V3))
looking3=data.frame(unique(ChFo$V3))
looking4=data.frame(unique(PaVa$V3))
looking5=data.frame(unique(GoMo$V3))
looking6=data.frame(unique(GoAq$V3))
looking7=data.frame(unique(GoPa$V3))
looking8=data.frame(unique(AcAu$V3))

#####Find unique repetitive element classifications
repeatlabels=c(looking1$unique.NeMu.V3.,looking2$unique.ChFu.V3.,looking3$unique.ChFo.V3.,looking4$unique.PaVa.V3.,looking5$unique.GoMo.V3.,looking6$unique.GoAq.V3.,looking7$unique.GoPa.V3.,looking8$unique.AcAu.V3.)

####unique labels were manually assigned colours, see below
colassigner=as.data.frame(unique(repeatlabels))

colassigner=read.csv(file.choose())

####Data format reassignment
set=list(NeMu = NeMu, ChFu = ChFu, ChFo = ChFo, PaVa=PaVa, GoMo = GoMo, GoAq= GoAq, GoPa=GoPa, AcAu= AcAu)


for (n in names(set)){
  newname=paste0(n, "Chrom")
    newdata=data.frame(Chrom=set[[n]]$V1, Start=set[[n]]$V4, End=set[[n]]$V5, Name=set[[n]]$V3)
  assign(newname, newdata)}



####Colour assignment based on TE type 

set2=list(NeMuChrom = NeMuChrom, ChFuChrom = ChFuChrom, ChFoChrom = ChFoChrom, PaVaChrom=PaVaChrom, GoMoChrom = GoMoChrom, GoAqChrom= GoAqChrom, GoPaChrom=GoPaChrom, AcAuChrom= AcAuChrom)

for (c in names(set2)){
  newdata=merge(set2[[c]], colassigner, by.x="Name", by.y="type", all.x=TRUE)
  assign(c,newdata)
}


#####Plotting TE types across sequence

####NeMu Chromplot unusable

chromPlot(bands=NeMuChrom,annot1=NeMuChrom, plotRndchr=T)


######ChFu Chromplot unusable

chromPlot(bands=ChFuChrom,annot1=ChFuChrom, plotRndchr=T)

####ChFo chromplot usable, with filtering of small fragments out, and estimation of shown nucleotides

ChFoChromsubset=subset(ChFoChrom, Chrom%in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8"))

chromPlot(bands=ChFoChromsubset,annot1=ChFoChromsubset, plotRndchr=T, bin=5000, title= "C. formosanus")

legend(x="bottomright", title="TE Type", legend=c("LINE", "Unknown LINE", "LTR-RT", "Unknown LTR-RT", "Low Complexity", "Helitron", "DNA Tranposon", "Non-LTR-RT", "TIR", "SINE", "Long Terminal Repeats", "Repeats", "Target Site Duplications"), fill=c("purple", "plum2", "turquoise", "blue", "navajowhite2", "green", "red", "yellow", "brown", "orange", "deeppink1", "pink", "black", cex=0.5))

ChFocount=as.list(unique(ChFoChromsubset$Chrom))
ntshown=0
for (i in ChFocount){
  chrdata=subset(ChFoChromsubset, Chrom== i)
  maxval=max(chrdata$End)
  ntshown=ntshown+maxval
  print(paste("ChFo shown:", ntshown))
}

#####PaVa chromplot usable, and estimation of shown nucleotides

chromPlot(bands=PaVaChrom,annot1=PaVaChrom, plotRndchr=T, bin=5000, title="P. varius")
legend(x="bottomright", title="TE Type", legend=c("LINE", "Unknown LINE", "LTR-RT", "Unknown LTR-RT", "Low Complexity", "Helitron", "DNA Tranposon", "Non-LTR-RT", "TIR", "SINE", "Long Terminal Repeats", "Repeats", "Target Site Duplications"), fill=c("purple", "plum2", "turquoise", "blue", "navajowhite2", "green", "red", "yellow", "brown", "orange", "deeppink1", "pink", "black", cex=1.5))

PaVacount=as.list(unique(PaVaChrom$Chrom))
ntshown=0
for (i in PaVacount){
  chrdata=subset(PaVaChrom, Chrom== i)
  maxval=max(chrdata$End)
  ntshown=ntshown+maxval
  print(paste("PaVa shown:", ntshown))
}

#####GoMo chromplot usable but with heavy filtering of fragments,  and estimation of shown nucleotides

GoMoChromsubset=subset(GoMoChrom, Chrom %in% c("940783","940784","940785","940786","940787"))

chromPlot(bands=GoMoChromsubset,annot1=GoMoChromsubset, plotRndchr=T, bin=5000, title="G. montsenyensis")
legend(x="bottomright", title="TE Type", legend=c("LINE", "Unknown LINE", "LTR-RT", "Unknown LTR-RT", "Low Complexity", "Helitron", "DNA Tranposon", "Non-LTR-RT", "TIR", "SINE", "Long Terminal Repeats", "Repeats", "Target Site Duplications"), fill=c("purple", "plum2", "turquoise", "blue", "navajowhite2", "green", "red", "yellow", "brown", "orange", "deeppink1", "pink", "black", cex=1.5))


GoMocount=as.list(unique(GoMoChromsubset$Chrom))
ntshown=0
for (i in GoMocount){
  chrdata=subset(GoMoChromsubset, Chrom== i)
  maxval=max(chrdata$End)
  ntshown=ntshown+maxval
  print(paste("GoMo shown:", ntshown))
}

#####GoAq chromplot usable

GoAquni=as.data.frame(unique(GoAqChrom$Chrom))

chromPlot(bands=GoAqChrom,annot1=GoAqChrom, plotRndchr=T, bin=5000, title= "G. aquaticus")

legend(x="bottomright", title="TE Type", legend=c("LINE", "Unknown LINE", "LTR-RT", "Unknown LTR-RT", "Low Complexity", "Helitron", "DNA Tranposon", "Non-LTR-RT", "TIR", "SINE", "Long Terminal Repeats", "Repeats", "Target Site Duplications"), fill=c("purple", "plum2", "turquoise", "blue", "navajowhite2", "green", "red", "yellow", "brown", "orange", "deeppink1", "pink", "black", cex=1.5))

##### GoPa chromplot is unusable 

gopascaf=subset(GoPaChrom, End>1500000)
unigopa=as.data.frame(unique(gopascaf$Chrom))
GoPaChromsubset=subset(GoPaChrom, Chrom%in%as.vector(unigopa$`unique(gopascaf$Chrom)`))

chromPlot(bands=GoPaChromsubset,annot1=GoPaChrom, plotRndchr=T)


####AcAu chromplot usable with filtering, and estimation of shown nucleotides
acauscaf=subset(AcAuChrom, End>10000000)
uniacau=as.data.frame(unique(acauscaf$Chrom))

AcAuChromsubset=subset(AcAuChrom, Chrom %in% c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4"))


chromPlot(bands=AcAuChromsubset,annot1=AcAuChrom, plotRndchr=T, bin=5000, title= "A. australiensis")


AcAucount=as.list(unique(AcAuChromsubset$Chrom))
ntshown=0
for (i in AcAucount){
  chrdata=subset(AcAuChromsubset, Chrom== i)
  maxval=max(chrdata$End)
  ntshown=ntshown+maxval
  print(paste("AcAu shown:", ntshown))
}
