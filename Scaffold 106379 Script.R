scaffoldID <- "Scaffold106379"

scaffold <- readLines("Output_FASTAs/Scaffold106379.fa")
scaffoldID <- grep(pattern = "^>", x = scaffold, value = T)
scaffoldID <- gsub(pattern = ">", replacement = "", x = scaffoldID)
# Apply UDF to find ORFs and translate these into AA
# --> we will use the AA sequences to conduct a BLAST analysis and identify ORFs coding for Aquaporin exons
tryCatch(
  {
    for(i in 1:length(scaffoldID)){
      findORFsTranslateDNA2AA(scaffold = scaffold, scaffoldID = scaffoldID[i], MinLen = 40)
    }
  })



# Open the CSV output of the UDF containing data on predicted ORFs and associated sequences 
orf.report <- read.csv("ORFs_report/Scaffold106379_ORFs.csv")

# How many ORFs were found?
nrow(orf.report) 

# What are the longest ORFs?
orf.report[order(-orf.report$width),] 

# Look at translated ORFs (=AA sequences)
translated.orfs <- readLines("AA_ORFs/Scaffold106379_ORFs.fa")
head(translated.orfs)


###
#Build map of scaffold with ORFs
###

#Read scaffold FASTA file (line by line)
scaffold <- readLines("Output_FASTAs/Scaffold106379.fa")

#List all files with ORFs
ORFfiles <- list.files(path = "ORFs_report", pattern = ".csv", full.names = T)

#Check if folder where ORF maps will be saved exists
# if not then creates it
output_dir <- file.path(paste0(getwd(), "/ORFs_map/"))

if(dir.exists(output_dir)){
  print(paste0("Dir", output_dir, " already exists!"))
}else{
  print(paste0("Created ", output_dir))
  dir.create(output_dir)  
}

#Produce a map (in pdf format) for each scaffold
for(i in 1:length(ORFfiles)){
  print(ORFfiles[i])
  #Read file in
  ORF <- read.csv(ORFfiles[i])
  
  #Process FASTA scaffold sequence
  seq <- strsplit(scaffold[grep(paste(">", ORF$scaffoldID[1], sep=''), scaffold)+1], split='')
  
  #Separate ORFs by strand
  ORFplus <- subset(ORF, ORF$strand == "+")
  ORFneg <- subset(ORF, ORF$strand == "-")
  
  #Create plot
  pdf(paste("ORFs_map/", ORF$scaffoldID[1], "_ORFs_annotated.pdf", sep=''))
  
  #Initiate plot
  plot(x=1, y=1, xlim=c(0,length(seq[[1]])), ylim=c(0,2), type='n', bty="n", axes=F, xlab = "", ylab='')
  #Add title
  text(x=5, y=2, paste(ORF$scaffoldID[1], " (strand: + in grey and - in blue)", sep=''), adj=0, cex=.8)
  
  #Create a segment with length of scaffold
  segments(x0=0, x1=length(seq[[1]]), y0=1, y1=1, col='black', lwd=3)
  
  #Add ORFs: rectangles (grey: +, blue: -)
  rect(xleft=ORFplus$start, xright=ORFplus$end, ybottom=0.75, ytop=1.25, col='grey')
  rect(xleft=ORFneg$start, xright=ORFneg$end, ybottom=0.75, ytop=1.25, col='blue')
  text(x=(ORFneg$start + ORFneg$end)/2, y=0.7, paste(ORFneg$ORFID, " (", ORFneg$start, ":", ORFneg$end, ")", sep=''), srt=90, col='blue', adj=1, cex=0.4)
  text(x=(ORFplus$start + ORFplus$end)/2, y=1.3, paste(ORFplus$ORFID, " (", ORFplus$start, ":", ORFplus$end, ")", sep=''), srt=90, col='black', adj=0, cex=0.4)
  
  #Add x axis
  axis(side = 1)
  mtext("Sequence (bp)", pos = c(0,0.5), side=1, line=2, cex.lab=0.6,las=1)
  #Close pdf
  dev.off()
}





##
## ORF_10 (2) - not aqp
# ORF_14 (6) - AQP
# ORF_21 (14) - AQP
# ORF_23 (16) - not aqp
# ORF_4 (19) - not aqp
# ORF_5 (20) - not aqp
# ORF_9 (24) - AQP

# 6, 14, 24 are AQP (14, 6, 24 on map)

# 2, 16, 19, 20 are not AQP
##
##

ORFs2protSeq <- function(fastaAA, ORFsToCat, stopCodon){
  #Check if the FASTA file exists
  checkFile <- file.exists(fastaAA)
  if(checkFile == FALSE){
    #Stop the function and print error
    stop(call. = FALSE, paste(fastaAA, "does not exists! Please check that your working directory is set properly.", sep = " "))
  }else{
    #Open FASTA file
    fas <- readLines(fastaAA)
    #Find ORFs in Fasta file
    ToFind <- match(ORFsToCat, fas)
    #Check that all the ORFs were found in the Fasta
    checkORFs <- is.na(ToFind)
    if(length(which(checkORFs == TRUE)) > 0){
      #Stop the function and print error
      stop(call. = FALSE, paste(paste(ORFsToCat[which(checkORFs == TRUE)], collapse = " "), "couldn't be found in", fastaAA, " Please fix the error(s).", sep = " "))
    }else{
      print(paste("Producing protein sequence with", paste(ORFsToCat, collapse = ", "), sep = " "))
      #Cat ORFs into AA sequence
      protSeq <- paste(fas[ToFind+1], collapse = "")
      if(is.logical(stopCodon) == TRUE){
        #Discard stop codons from final sequence
        protSeq <- gsub("[*]", "", protSeq)  
      }
    }
    #Produce final protein sequence in FASTA format
    # header
    head <- paste(">Protein_sequence_from_", paste(gsub("[>]", "", ORFsToCat), collapse = "|"), sep= "")
    protSeqOUT <- paste(head, protSeq, sep = "\n")
    print(protSeqOUT)
  }
  return(protSeqOUT)
}



#Name of the FASTA file (make sure the path is correct)
fastaAA <- "AA_ORFs/Scaffold106379_ORFs.fa"

#Fasta headers of the ORFs you want to concatenate
# Here some ORFs names are faulty
ORFsToCat <- c(">Scaffold106379_ORF_21", ">Scaffold106379_ORF_14", ">Scaffold106379_ORF_9")

#Stop codon
stopCodon <- TRUE

#Apply UDF
ORFs2protSeq(fastaAA, ORFsToCat, stopCodon)


##
#AA sequence
##
#User copy AA sequence in the object below
# Here, associated to Scaffold106379_ORF_(21, 14, 9)
AAseq <- "nMVDVVYIQVGELAGIAVGATVMLNILVAGPSSGGSMNPVRTLGPAVAAGNYTVIWIYLLAPTLGALAGAGVYTLVKLNAVEGEEPRAVRSFRRMQLGAEFVGTFILIFAATAGPIVNQKYANSETLIGNAACAGLAVMIIILSTGHISGAHLNPSLTIAFAALRHFPWAQVPAYILAQVSASICASFALKGVFHPFMSGGVTVPTVSTGQAFALEFLITFNLLFVVTAVATDTRAVSITYDSFSFFLFLLNLMESGGDGQTPPASAPATPGTPLPLFPASRVDSLSYERKSMPRCMCLPVHAPSWGASHSCFTDFPSPDISLTRKVHRTPPVFCIF*"
#21 is
#14 is
#9 is 

##
#Create a data.frame to store data on
# - N_NPA: Number of NPA motifs
# - NPA_motif: NPA motif
# - StartNPA: Starting position of NPA motif in AA seq
##
NPAdat <- data.frame("N_NPA" = numeric(length(AAseq)), "NPA_motif" = character(length(AAseq)), "StartNPA" = numeric(length(AAseq)))

#Start by matching to find NP motif and locations along sequence
tmp <- Biostrings::matchPattern(pattern = "NP", Biostrings::AAString(AAseq)) 
tmp <- as.data.frame(IRanges::ranges(tmp))

if(nrow(tmp) > 0){
  #Print that it found NP motif
  print("NP found")
  #There is an NP motif
  #N NP
  NPAdat$N_NPA <- nrow(tmp)
  #Infer motifs
  NPAdat$NPA_motif <- paste(paste(rep("NP", nrow(tmp)), strsplit(as.vector(AAseq), split='')[[1]][as.numeric(tmp[,2]+1)], sep = ''), collapse = '/')
  #Start motifs
  NPAdat$StartNPA <- paste(tmp[,1], collapse = '/')
}else{
  #There isn't an NP motif, so we look for a PA motif
  tmp <- Biostrings::matchPattern(pattern = "PA", Biostrings::AAString(AAseq))
  tmp <- as.data.frame(IRanges::ranges(tmp))
  if(nrow(tmp) > 0){
    #Print that it found PA motif
    print("PA found")
    #N NP
    NPAdat$N_NPA <- nrow(tmp)
    #Infer motifs
    NPAdat$NPA_motif <- paste(paste(strsplit(as.vector(AAseq), split='')[[1]][as.numeric(tmp[,1]-1)], rep("PA", nrow(tmp)), sep = ''), collapse = '/')
    #Start motifs
    NPAdat$StartNPA <- paste(tmp[,1]-1, collapse = '/')
  }else{
    next
  }
}

print(paste("The following NPA motifs were found for AA sequence:", AAseq, sep=" "))


print(NPAdat)