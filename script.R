# Marta's thesis script

# Packages

## We load all the needed packages before starting the analysis

packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
              "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw",
              "BSgenome.Scerevisiae.UCSC.sacCer3")
suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))

# We load files

E_Ori <- read.table("/Users/rodrigo/Desktop/Analysis_pack/E_Rep.bed", header = F, quote = "/t")[ ,1:4]
L_Ori <- read.table("/Users/rodrigo/Desktop/Analysis_pack/L_Rep.bed", header = F, quote = "/t")[ ,1:4]
All_Ori <- read.table("/Users/rodrigo/Desktop/Analysis_pack/OriginList_Full.bed", header = TRUE, quote = "/t")
All_Ori_Link <- "/Users/rodrigo/Desktop/Analysis_pack/OriginList_Full.bed"
image <- load.image("/Users/rodrigo/Desktop/Analysis_pack/ForkPic_1.jpeg")


colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); E_Ori$mid <- round((E_Ori$chromStart + E_Ori$chromEnd)/2)
colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); L_Ori$mid <- round((L_Ori$chromStart + L_Ori$chromEnd)/2)
# run the alignment 

# We introduce the experiment name, this name will be the name our folder containing results will have. 

Pro_1 <- "Test"

# First quality check
if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "QR", ".html"))){
  
  fls = c("/Users/rodrigo/Desktop/MartaV/analysis/25_S25_R1_001.fastq.gz", "/Users/rodrigo/Desktop/MartaV/analysis/25_S25_R2_001.fastq.gz", 
          "/Users/rodrigo/Desktop/MartaV/analysis/26_S26_R1_001.fastq.gz", "/Users/rodrigo/Desktop/MartaV/analysis/26_S26_R2_001.fastq.gz", 
          "/Users/rodrigo/Desktop/MartaV/analysis/29_S29_R1_001.fastq.gz", "/Users/rodrigo/Desktop/MartaV/analysis/29_S29_R2_001.fastq.gz")
  
  names(fls) = sub(".fastq", "", basename(fls))
  
  qas = lapply(seq_along(fls),
               function(i, fls) qa(readFastq(fls[i]), names(fls)[i]),
               fls)
  qa = do.call(rbind, qas)
  rpt = report(qa, dest = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "QR", ".html"))
  
}

# We create a function to run the alignments

RunAlign_paired <- function(File_R1, File_R2, SampName){
  
  # create a temporary directory to those transitory files we won't need to store
  
  tempdir(check = TRUE)
  
  Sam <- tempfile(fileext = ".sam") #sam file
  Bam <- tempfile(fileext = ".bam") 
  nmCollate <- tempfile(fileext = ".bam")
  fixMat <- tempfile(fileext = ".bam")
  SrtBam <- tempfile(fileext = ".bam")
  
  
  Pro_1 <- Pro_1
  Pro_2 <- SampName
  
  suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1))) # We create a new directory in our computer desktop named as the named previously chosen "Pro_1"
  
  suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Bam"))) # We create another directory inside the previous one to store the bam files.
  
  AlnLog <- paste0("~/Desktop/", Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".log") # We will store the information in a .log file created in our BAM directory
  SFBam <- paste0("~/Desktop/", Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".bam") # .bam file containing the alignment information we will store in our directory
  
  #read the indexed reference genome for the alignment of sequenced data
  ref_index <- "/Users/rodrigo/Desktop/Analysis_pack/bowtie2-2.4.4-macos-x86_64/indexes/S288C_Ref" # This indexed genome has previously been created
  
  #following commands will run the alignment, check quality, sort, filter and index the resultant bam file 
  
  system(sprintf("(/Users/rodrigo/Desktop/Analysis_pack/bowtie2-2.4.4-macos-x86_64/bowtie2 -p 16  --no-discordant --fr -x %s -1 %s -2 %s -S %s) 2> %s", 
                 ref_index, File_R1, File_R2, Sam, AlnLog)) # We make the alignment using bowtie2
  
  system(sprintf("/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools view -bS -@ 15 -q 30 -f 2 %s > %s", Sam, Bam)) # We convert our sam file into bam
  
  system(sprintf("/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools collate -@ 15 -o %s %s", nmCollate, Bam)) # We group our reads by name
  
  system(sprintf("/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools fixmate -@ 15 -m %s %s", nmCollate, fixMat)) # Fixes mate errors after alignment
  
  system(sprintf("/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools sort -l 9 -@ 15 -m 1024M  -O bam -o %s %s", SrtBam, fixMat)) # We sort our reads by coordinates
   
  system(sprintf("/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools markdup -@ 15 %s %s", SrtBam, SFBam)) # We mark duplicates from a coordinate sorted file
  
  system(sprintf("/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools index -@ 15 %s", SFBam)) # We index the bam file 
  
  unlink(c(Sam, Bam, nmCollate, fixMat, SrtBam), recursive = T, force = T) # We delete the files inside the temporary directory
  
}

# We run the alignments of the files containing our data

RunAlign_paired(File_R1 = "/Users/rodrigo/Desktop/MartaV/analysis/25_S25_R1_001.fastq.gz",
                File_R2 = "/Users/rodrigo/Desktop/MartaV/analysis/25_S25_R2_001.fastq.gz",
                SampName = 'Input')

RunAlign_paired(File_R1 = "/Users/rodrigo/Desktop/MartaV/analysis/26_S26_R1_001.fastq.gz",
                File_R2 = "/Users/rodrigo/Desktop/MartaV/analysis/26_S26_R2_001.fastq.gz",
                SampName = 'BrDU')

RunAlign_paired(File_R1 = "/Users/rodrigo/Desktop/MartaV/analysis/29_S29_R1_001.fastq.gz",
                File_R2 = "/Users/rodrigo/Desktop/MartaV/analysis/29_S29_R2_001.fastq.gz",
                SampName = 'SUP')


## Get the bin 
## To get the bin, we calculate the mean of the fragment size from each bam file obtained after alignment, and then we get the mean of the means from every bam file

bamFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Bam", "/", "*", ".bam")) # bamFiles are those files containing .bam

  mFragS <- c()
  for(i in 1:length(bamFiles)) {
    binMFS <- round_any(mean(c(mFragS, mean(getPESizes(bamFiles[i], param=readParam(pe="both"))$sizes))), 100)
  }
  binSize <- binMFS
  
  # binSize is 300


## Coverage
# We define bamFiles as all the files ending with .bam and present in the folder
  
  # We create a new function to get the coverage
  BamCoverage <- function(bamFile, binSize = 300, stepSize = 10, slidingWindow = TRUE){
    # The parameters inside the function are the bamFile, the binSize that has been fixed as 300 (we know its value), stepsize that we fixed as 10 and the presence of the sliding window
    Pro_1 <- unlist(strsplit(tools::file_path_sans_ext(basename(bamFile)), split='_', fixed=TRUE))[[1]] #extract protein name
    Pro_2 <- unlist(strsplit(tools::file_path_sans_ext(basename(bamFile)), split='_', fixed=TRUE))[[2]] #extract sample name
    
    tempdir(check = TRUE)
    # We store temporary files in a temporary directory
    GenomFile <- tempfile(fileext = ".txt")
    binFile <- tempfile(fileext = ".bed")
    # We also use awk language to process the created files obtained, we select values over 0
    # idxstats prints the index information from the indexed bam created by using samtools index
    command_1 <- "/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools idxstats %s | awk 'BEGIN {OFS=\"\\t\"} {if ($2>0) print ($1,$2)}' >  %s" # We get the indexed bam file
    system(sprintf(command_1, bamFile, GenomFile))
    
    if(slidingWindow == TRUE){
      command_2 <- "/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools makewindows -g %s -w %s -s %s > %s" # Creates interval windows across a genome
      system(sprintf(command_2, GenomFile, binSize, stepSize, binFile))
    } else {
      command_2 <- "/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools makewindows -g %s -w %s > %s"
      system(sprintf(command_2, GenomFile, binSize, binFile))
    }
    
    pncFiles_watson <- tempfile(fileext = ".bed") # We store the Watson information in a file placed into a temporary file, also the crick info
    pncFiles_crick <- tempfile(fileext = ".bed")
  # Those rows containing a quality score below 30 are removed, not taken into account
    # bedtools genomecov to sotore coverage characteristics in bed format
    command_3 <- "/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s | grep -v XS:i: | /Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools view -@ 8 -b - | /Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools genomecov -5 -d -ibam stdin -strand + | awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' > %s"
    command_4 <- "/Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s | grep -v XS:i: | /Users/rodrigo/Desktop/Analysis_pack/samtools-1.13/samtools view -@ 8 -b - | /Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools genomecov -5 -d -ibam stdin -strand - | awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' > %s"
    
    system(sprintf(command_3, binFile, bamFile, paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"), pncFiles_watson))
    system(sprintf(command_4, binFile, bamFile, paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"), pncFiles_crick))
    
    suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Coverage"))) # We create a new directory to store coverage bed files, watson and crick separated
    
    finFiles_watson <- paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "watson.bed")
    finFiles_crick <- paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "crick.bed")
    # the last command to get map the results and obtain coverage bed files 
    command_5 <- "/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools map -a %s -b %s -null 0 -o sum | awk 'BEGIN {OFS=\"\\t\"} {if ($4>=0) print $1,$2,$3,\"%s\",$4}' > %s"
    
    system(sprintf(command_5, binFile, pncFiles_watson, paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"), finFiles_watson))
    system(sprintf(command_5, binFile, pncFiles_crick, paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"), finFiles_crick))
    
    unlink(c(GenomFile, binFile, pncFiles_watson, pncFiles_crick), recursive = T, force = T) # Delete those files in the temporary directory
    
  }
  # Run the function for all the bam files
  for(i in 1:length(bamFiles)){
    BamCoverage(bamFile = bamFiles[i])
  }
  

## Calculate ratio
  # Two step normalization

  CoverageFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", "*", ".bed"))
  # Define a new function
  CalculateRatio <- function(IP_coverage, Input_coverage){
    
    IP.df = read.table(IP_coverage, header = F)
    In.df = read.table(Input_coverage, header = F)
    
    suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Ratios"))) 
    
    # First normalization step: equalize INPUT and IP sizes to make the comparable
    IP_Sum <- sum(as.numeric(IP.df[,5]))
    In_Sum <- sum(as.numeric(In.df[,5]))
    corrFactor <- IP_Sum/In_Sum # To equalize we get the correction factor as the division of the IP between the INPUT coverage scores
    Ratio <- round(IP.df[,5]/In.df[,5]/corrFactor, 4) # We get the ratio from multiplying the IP coverage score by the correction factor and dividing by the input coverage score
    In.score.norm <- round(In.df[,5]*corrFactor) # The input normalization is the input coverage multiplied by the correction factor
    Ratio[!is.finite(Ratio)] <- 0
    
    strand <- cbind.data.frame(IP.df[,1], IP.df[,2], IP.df[,3], IP.df[,4], IP.df[,5], In.score.norm, Ratio)
    
    chroms <- unique(strand[,1])
    s_strand <- NULL
    for(i in 1:length(chroms)){
      chr <- strand[strand[,1]==chroms[i], ]
      x <- chr[,2]
      y <- chr$Ratio
      splineObject <- smooth.spline(x, y) # Fits a cubic smoothing spline to the supplied data, stored in column 8 in the final output
      chr$splineSmooth <- round(as.numeric(splineObject$y), 3)
      s_strand <- rbind.data.frame(s_strand, chr)
    }
    # We use Poison distribution to establish those values being statistically significant
    s_strand$ppois <- ppois(  q=s_strand[,5] - 1, 
                              lambda=s_strand$`In.score.norm`,  # The lambda from the Poison refers to the previously calculated score normalization from the input
                              lower.tail=FALSE, log=FALSE      ) # Upper tail in the distribution
    
    # We add column names to the final file
    colnames(s_strand) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'ip.score', 'in.score', 'ratio', 'smooth', 'pvalue')
# We write a new table containing the ratios and that will be stored in a new folder named Ratios
    write.table(s_strand, paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", s_strand$name[1], ".bed"),
                quote=FALSE, row.names=FALSE, sep="\t")
    
  }
  # We need to get the ratios from the files containing information about our protein, therefore we get the ratios against the INPUT, representing the rest of the DNA
  # BrdU against input
  CalculateRatio(CoverageFiles[1], CoverageFiles[3]) # crick
  CalculateRatio(CoverageFiles[2], CoverageFiles[4]) # watson
  # sup against input
  CalculateRatio(CoverageFiles[5], CoverageFiles[3]) #crick
  CalculateRatio(CoverageFiles[6], CoverageFiles[4]) # watson
  # Input against input, not necessary, just a control in this case
  CalculateRatio(CoverageFiles[3], CoverageFiles[3]) # crick
  CalculateRatio(CoverageFiles[4], CoverageFiles[4]) # watson
  


# Peak calling

  bamFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Bam", "/", "*", ".bam")) # We storing bam files from the folder as bamFiles
  
  ProcessPeaks <- function(bamFiles){ # We define a new function to obtain the peaks files
    # Inside the first function, we define a second one to call the peaks
    PeakFinder <- function(IPBam, InBam){
      # Differenciate between input, the one used to normalise and IP
      Input_combined <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools bamtobed -i %s", InBam)) ) # We convert the bam to a bed file
      # Positive are watson and negative are crick
      Input_Plus <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"+\") print $0}'", InBam)) ) # We select those rows containing a positive value in the sixth column // awk language
      Input_Minus <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"-\") print $0}'", InBam)) )  # We select those rows containing a negative value in the sixth column
      
      IP_combined <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools bamtobed -i %s", IPBam)) )
      IP_Plus <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"+\") print $0}'", IPBam)) )
      IP_Minus <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"-\") print $0}'", IPBam)) )
      # Another function to obtain peaks placed at the origin
      OriginPeaks <- function(IP_DF, Input_DF){
        
        tempdir(check = TRUE)  ## We create a temporary directories
        IP <- tempfile(fileext=".bed")
        Input <- tempfile(fileext=".bed")
        outDir <- tempdir()
        peakFile <- tempfile(fileext = ".bed") # temporary file containing peaks
        ColHeads <- "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\"" # We name the columns in the file
        # Define the tables in which we are gonna store data
        write.table(IP_DF, file = IP, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        write.table(Input_DF, file = Input, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        # We use macs2 algorithm to obtain peaks, defining t: the data file, c: the mock data file, f: the format of the input file, g: genome size, 
        # n: prefix string, outdir: where to store files, nomodel: not to build the shifting model
        system(sprintf("macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n %s --outdir %s 2> /dev/null", IP, Input, "Peak", outDir))
        allPeaks <- read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#") # List of all the peaks stores into the temporary directory outDir
        
        write.table(allPeaks, file = peakFile, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE) # We define the allPeaks table containing information
        
        Peaks_at_Origins <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'", 
                                                    peakFile, All_Ori_Link, ColHeads)), header = TRUE ) ## We use intersect to screen for overlaps between two sets
        Peaks_at_Origins <- Peaks_at_Origins[!duplicated(Peaks_at_Origins$oriName), ] # To find those duplicated elements
        # We define the origin center as the mean between the start and the end
        Peaks_at_Origins$oriCenter <- round((Peaks_at_Origins$oriStart + Peaks_at_Origins$oriEnd)/2) # We add the origin center as the intermediate point between oriStart and oriEnd
        
        allPeaks <- allPeaks[,c(1,2,3,4,5)]
        
        names(allPeaks) <- c("chrom", 'peakStart', 'peakEnd', 'peakLength', 'peakSummit') # We name the columns of the file allPeaks
        # By using antijoin we select those rows not shared between both dataframes
        ALLbutROs <- anti_join(allPeaks, Peaks_at_Origins, by = c("chrom", 'peakStart', 'peakEnd'))
        
        PeakList <- list(Peaks_at_Origins, allPeaks, ALLbutROs) #Store in a list the results from the othes files
        
        names(PeakList) <- c("ROpeaks", "ALLpeaks", "ALLbutROs") # We name the files
        
        return(PeakList) # We tell them to give us the list of peaks got 
        
        unlink(c(IP, Input, peakFile)) # We remove temporary directories
        unlink(outDir, recursive = TRUE)
        
      }
      
      Watson <- OriginPeaks(IP_DF=IP_Plus, Input_DF=Input_Plus) # We define the Watson as the data contained in the previously selected positive values, same and opposite with crick
      Crick <- OriginPeaks(IP_DF=IP_Minus, Input_DF=Input_Minus)
      Combo <- OriginPeaks(IP_DF=IP_combined, Input_DF=Input_combined) # Both w and c
      
      ResPeaks <- list(Watson, Crick, Combo) # All the peaks in the three lists
      names(ResPeaks) <- c("watsonPeaks", "crickPeaks", "comboPeaks")
      
      return(ResPeaks)
      
    }
    
    IP_Input <- PeakFinder(IPBam=bamFiles[1], InBam=bamFiles[2]) # We use the recently defined function to get the peaks from the IP and INPUT
    ## 
##    #Define overlapping peaks
    
    tempdir(check = TRUE)
    
    pF_1 <- tempfile(fileext = ".bed")
    pF_2 <- tempfile(fileext = ".bed")
    
    #Overlaps between both watson and crick strands
    # 
    IP_ColHeads <- "\"chrom\\tBWpStart\\tBWpEnd\\tBWpSummit\\tBCpStart\\tBCpEnd\\tBCpSummit\\toriName\\toriCenter\""
    
    write.table(IP_Input$watsonPeaks$ROpeaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE) #We store the peaks from watson and placed at origins
    write.table(IP_Input$crickPeaks$ROpeaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)  #We store the peaks from crick and placed at origins
    Overlapping_RO_Peaks <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$11,$12,$14,$15,$18}'",
                                                    pF_1, pF_2, IP_ColHeads)), header = TRUE ) #overlapping peaks by using intersect of the file %s, -wa: Write the original entry in A for each overlap. -wb: Write the original entry in B for each overlap. -a and -b is the file. Print separated by tabs.
    Overlapping_RO_Peaks <- Overlapping_RO_Peaks[!duplicated(Overlapping_RO_Peaks$oriName), ]
    
    Primary_RO_Peaks <- IP_Input$comboPeaks$ROpeaks # Both w and c at replication origins, combopeaks: both W and C
    
    
    ##
    IP_ColHeads <- "\"chrom\\tBWpStart\\tBWpEnd\\tBWpSummit\\tBCpStart\\tBCpEnd\\tBCpSummit\""
    
    write.table(IP_Input$watsonPeaks$ALLpeaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE) # We store this into the file inside the temporary directory
    write.table(IP_Input$crickPeaks$ALLpeaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
    Overlapping_ALL_Peaks <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$7,$8,$10}'",
                                                     pF_1, pF_2, IP_ColHeads)), header = TRUE ) # Same as twice before
    
    Primary_ALL_Peaks <- IP_Input$comboPeaks$ALLpeaks # Both w and c not only at replication origins
    
    
    ##
    IP_ColHeads <- "\"chrom\\tBWpStart\\tBWpEnd\\tBWpSummit\\tBCpStart\\tBCpEnd\\tBCpSummit\""
    
    write.table(IP_Input$watsonPeaks$ALLbutROs, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE) # All except for replication origins 
    write.table(IP_Input$crickPeaks$ALLbutROs, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
    Overlapping_ALLminusRO_Peaks <- read.table(pipe(sprintf("/Users/rodrigo/Desktop/Analysis_pack/bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$7,$8,$10}'",
                                                            pF_1, pF_2, IP_ColHeads)), header = TRUE )
    
    Primary_ALLminusRO_Peaks <- IP_Input$comboPeaks$ALLbutROs
    
    
    unlink(c(pF_1, pF_2))
    
    #
    ###Calculate Fork Locations from BrDU data
    # As BrdU is in those recently replicated regions, it is possible to stablish replication origins with its information
    SignalRanges <- function(PeakList){
      # Got by the peak list previously stored
      
      Chroms <- paste0("chr", as.roman(1:16)) # chromosomes from 1 to 16
      
      ForkPos <- NULL
      
      
      for(i in 1:16){
        
        PeakPos <- NULL
        
        i <- i
        
        Chr_Peaks <- PeakList[PeakList$chrom == Chroms[i], ] # Peak list containing a column named chrom refers to chromosome
        # Creating a loop
        
        if(length(Chr_Peaks$chrom)==0) next 
        
        
        for(y in 1:length(Chr_Peaks$chrom)){
          
          y <- y
          # Same for w and c but separetely
          WatsonLeft <- round((Chr_Peaks$BWpSummit[y] - Chr_Peaks$BWpStart[y])) # From the middle to the start, left arm
          WatsonRight <- round((Chr_Peaks$BWpEnd[y] - Chr_Peaks$BWpSummit[y])) # From the middle to the end, right arm
          
          CrickLeft <- round((Chr_Peaks$BCpSummit[y] - Chr_Peaks$BCpStart[y]))
          CrickRight <- round((Chr_Peaks$BCpEnd[y] - Chr_Peaks$BCpSummit[y]))
          # All peaks from w and c:
          Peaks <- cbind.data.frame(WatsonLeft=WatsonLeft, WatsonRight=WatsonRight, 
                                    CrickLeft=CrickLeft, CrickRight=CrickRight)
          
          
          PeakPos <- rbind.data.frame(PeakPos, Peaks)
          PeakPos[PeakPos < 0] <- 0
          
        }
        
        ForkPos <- rbind.data.frame(ForkPos, PeakPos) # The position of the fork is given by the data of the brdu location in the genome
        
      }
      return(ForkPos)
    }
    
    if(dim(Overlapping_ALL_Peaks)[1]>0){
      ##Lagging and Leading strand lengths synthesized 
      ForkPos <- SignalRanges(Overlapping_ALL_Peaks)
      
      LaggLeadSynthesis <- cbind.data.frame(lagging=ForkPos$CrickLeft+ForkPos$WatsonRight, 
                                            leading=ForkPos$WatsonLeft+ForkPos$CrickRight)
      
      
      LeadingAverage <- round(mean(LaggLeadSynthesis$leading))
      LaggingAverage <- round(mean(LaggLeadSynthesis$lagging))
      LeadingSd <- round(sd(LaggLeadSynthesis$leading))
      LaggingSd <- round(sd(LaggLeadSynthesis$lagging))
      
      ###Estimate Averaging Window
      
      #Remove extreme outliers from the data, w and c
      Left <- c(ForkPos$WatsonLeft, ForkPos$CrickLeft)
      Right <- c(ForkPos$WatsonRight, ForkPos$CrickRight)
      
      Q <- quantile(Left, probs=0.99, na.rm = T) 
      I <- IQR(Left) # interquartile range
      up  <-  Q + 1.5*I # Upper Range
      NewLeft <- Left[which(Left < up)] # those values below the outliers
    
      Q <- quantile(Right, probs=0.99, na.rm = T)
      I <- IQR(Right)
      up  <-  Q + 1.5*I # Upper Range
      NewRight <- Right[which(Right < up)]
      
      #Range for Averaging Window using the outliers values
      LeftLim <- max(NewLeft) #Left pan
      RightLim <- max(NewRight) #Right pan
      
      #Round the window
      LeftSide <- round(LeftLim/1000+0.5)*1000
      RightSide <- round(RightLim/1000+0.5)*1000
      # Stablish average window if both limits are the same or different
      if(LeftSide==RightSide){
        AveragingWindow <- LeftSide
      } else {
        AveragingWindow <- max(c(LeftSide, RightSide))
      } 
    }
    
    #
    StatDat <- rbind.data.frame(LeftLim, RightLim, paste0(LeadingAverage, "±", LeadingSd), paste0(LaggingAverage, "±", LaggingSd), 
                                AveragingWindow, binSize, stepSize) # binds information summary of the parameters obtained
    
    rownames(StatDat) <- c('Left', 'Right', "LeadSynthesis", "LaggSynthesis", "AveragingWindow", 'bin', 'slide')
    colnames(StatDat) <- " "
    
    
    #
    # Create the files containing these results
    dir.create(paste0("~/Desktop/", Pro_1, "/", "Peaks"), showWarnings = FALSE)
    # origin peaks
    write.table(Primary_RO_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Primary_RO_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(Overlapping_RO_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Overlapping_RO_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
    # all peaks
    write.table(Primary_ALL_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Primary_ALL_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(Overlapping_ALL_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Overlapping_ALL_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
    # not origin peaks
    write.table(Primary_ALLminusRO_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Primary_ALLbutOri_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(Overlapping_ALLminusRO_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Overlapping_ALLbutOri_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
    
    # summary info
    write.table(StatDat, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Synthesis.bed"), quote=FALSE, sep="\t")
    
  }
  
  ProcessPeaks(bamFiles) # Call the functions by using the bamfiles


###
# we call the files depending on their type as different ways to use these names later in every function we define
peakFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", "*", "_ALL_Peaks.bed"))
watsonFiles <-  Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", "*", "watson.bed"))
crickFiles <-  Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", "*", "crick.bed"))

BrdU_watson <- read.table(watsonFiles[1], header = T)
BrdU_crick <- read.table(crickFiles[1], header = T)
# Create new datasets from the brdu data
BrdU_coverage <- cbind.data.frame(chrom=BrdU_watson$chrom, chromStart=BrdU_watson$chromStart, chromEnd=BrdU_watson$chromEnd, wat.score=BrdU_watson$ip.score, cri.score=BrdU_crick$ip.score)
BrdU_ratio <- cbind.data.frame(chrom=BrdU_watson$chrom, chromStart=BrdU_watson$chromStart, chromEnd=BrdU_watson$chromEnd, wat.score=BrdU_watson$ratio, cri.score=BrdU_crick$ratio)

BrdU_Peaks <- read.table(peakFiles[2], header = T)


#10- Plot IP profiles
# The next step after obtaining files containing data and information is to use this info to create plots that allow the researcher to see the representation and interpret it 

stepSize <- 10

BrdU_Plot_function <- function(CoverageFile, peakFile, DataType, SamplType){
#   
  CoverageFile = CoverageFile
  peakFile = peakFile
#   
  scoreVals <- c(CoverageFile$wat.score, CoverageFile$cri.score) # The scores are stored in the coverage files separated between watson and crick
#   
  Q <- quantile(scoreVals, probs=c(.01, .99), na.rm = FALSE) # Quantile of the score values defined by watson and crick scores
  I <- IQR(scoreVals) # Interquartile range
  up  <-  Q[2]+1.5*I # Upper Range  
  low <- Q[1]-1.5*I # Lower Range
  Scores <- scoreVals[which(scoreVals < up & scoreVals > low)] # Select those scores between the upper an lower range
  thd <- round(mean(Scores)+(12*sd(Scores))) # To establish limits in the plotting later
  Ylim_sc <- c(-thd, thd) # This will be the ylim from the y axis in the plot
  
  wcVals <- log2(CoverageFile$wat.score/CoverageFile$cri.score) # The value of the watson over crick is the logarithm of the watson coverage score between the crick
  wcVals[!is.finite(wcVals)] <- 0
#   
  Q <- quantile(wcVals, probs=c(.01, .99), na.rm = FALSE) 
  I <- IQR(wcVals)
  up  <-  Q[2]+1.5*I # Upper Range  
  low <- Q[1]-1.5*I # Lower Range
  wc <- wcVals[which(wcVals < up & wcVals > low)]
  wcs <- round(mean(wc)+(3*sd(wc)))
  Ylim_wc <- c(-wcs, wcs)
#   
  yAxis_reads <- c(-0.3*thd, -0.6*thd, -0.9*thd, 0*thd, 0.3*thd, 0.6*thd, 0.9*thd) # Y axis for the reads
  yAxis_wc <- c(-0.3*wcs, -0.6*wcs, -0.9*wcs, 0*wcs, 0.3*wcs, 0.6*wcs, 0.9*wcs)  # Y axis for the w over c
#   # Chromosome names
  Coverage_chr <- CoverageFile[CoverageFile$chrom == seqnames(Scerevisiae)[k], ]  
  All_Ori_chr <- All_Ori[All_Ori$chrom == seqnames(Scerevisiae)[k], ]
  Peaks2Plot_chr <- peakFile[peakFile$chrom == seqnames(Scerevisiae)[k], ]
#   
  Coverage <-  Coverage_chr[Coverage_chr$chromStart>=S & Coverage_chr$chromStart<=E, ] 
  Ori_chr <- All_Ori_chr[All_Ori_chr$chromStart>=S & All_Ori_chr$chromStart<=E, ]
  PeakReg <- Peaks2Plot_chr[Peaks2Plot_chr$peakStart>=S & Peaks2Plot_chr$peakStart<=E, ]
#   
  CovWat <- Coverage$wat.score # watson score as the coverage of watson
  CovCri <- Coverage$cri.score
#   
#   ###
  steps <- 10 # We set step as 10 as previously done with stepsize
#   
  if(stepSize == steps){ # If step is the same as stepsize we mantain all the values from the covs
    CovWat <- CovWat
    CovCri <- CovCri
  }
  if(stepSize > steps){
    s <- stepSize/steps # If stepsize higher than 10, we get the ratio by dividin them 
    CovWat <- as.vector(sapply(CovWat, function (x) rep(x,s))) # sapply() function takes list, vector or data frame as input and gives output in vector or matrix
    CovCri <- as.vector(sapply(CovCri, function (x) rep(x,s))) # rep- replicate using the S got by the division of the stepsize and step
  }
  if(stepSize < steps){
    s <- round(steps/stepSize)
    CovWat <- round(as.vector(tapply(CovWat, gl(length(CovWat)/s, s), mean))) # tapply() computes a measure (mean, median, min, max, etc..) or a function for each factor variable in a vector
    CovCri <- round(as.vector(tapply(CovCri, gl(length(CovCri)/s, s), mean))) # with gl() we generate factors with the pattern given by the s, dividing the number of elements in the cov of crick by it, and the s (from de division)
  }
#   
#   ### 
#   
  WCrat <- log2(CovWat/CovCri) # Watson over crick ratio by dividing them
  WCrat[!is.finite(WCrat)] <- 0
#   Smooth spline of the w, c and ratio
  Cov.wat <- round(as.numeric(smooth.spline(1:length(CovWat), CovWat)$y), 3) # smoothing the coverage watson, this minimizes the residual square sum. The arguments are the length of the data and the data column y
  Cov.cri <- round(as.numeric(smooth.spline(1:length(CovCri), CovCri)$y), 3)
#   
  WC.rat <- round(as.numeric(smooth.spline(1:length(WCrat), WCrat)$y), 3)
#  if the lengths of the smooth splined covs is less than the length per row divided by the steps (10) minus 1, the covs are changed as shown 
  if(length(Cov.wat) < (Length_per_Row/steps - 1)){
    Cov.wat <- c(Cov.wat, rep(NA, (Length_per_Row/steps - 1)-length(Cov.wat)) ) # By replicating tehe previously mentioned operation minus the length of the cov in watson
  }
#   
  if(length(Cov.cri) < (Length_per_Row/steps - 1)){
    Cov.cri <- c(Cov.cri, rep(NA, (Length_per_Row/steps - 1)-length(Cov.cri)) )
  }
#   
  if(length(WC.rat) < (Length_per_Row/steps - 1)){
    WC.rat <- c(WC.rat, rep(NA, (Length_per_Row/steps - 1)-length(WC.rat)) )
  }
#   
#   #plotting parameters:
  par(mar = c(0,0,0,0))
  suppressWarnings( # Plotting coverage watson
    plot(Cov.wat, type='h', ylim=Ylim_sc, col =  rgb(100,0,0,alpha=180, maxColorValue=255), ylab=' ', # h: vertical lines
         xlab=' ', xaxt='n', yaxt='n', lwd=0.07, bty = 'n',  cex.lab=1, las = 2, xaxs='i')
  ) 
  lines(Cov.cri*(-1), type='h', col =  rgb(0,100,0,alpha=180, maxColorValue=255), lwd=0.07) # coverage crick
  # different colors, different rgb
#   
  segments(x0 = c(Coverage$chromStart - S)/steps, x1 = c(Coverage$chromEnd - S)/steps, y0 = 0, y1= 0, lwd = 0.5)
#   
#  #plot left y axis
  axis(side = 2, at = yAxis_reads, labels = round(yAxis_reads), line = 0, tick = TRUE, lwd.ticks = 1.5, las = 2, cex.axis = 0.8)
  abline(h=yAxis_reads, lwd=0.05, col =  rgb(112,128,144,alpha=225, maxColorValue=255)) # dark col, tick axis drawn
#   
#   #draw the peaks
  if(DataType=="strandedRatio"){
    if(length(PeakReg$peakStart) > 0){
      for(i in 1:length(PeakReg$peakStart)){
        segments(x0 = c(PeakReg$peakStart[i] - S)/steps, x1 = c(PeakReg$peakEnd[i] - S)/steps, 
                 y0 = par('usr')[4]-(thd*0.1), y1= par('usr')[4]-(thd*0.1), lwd = 4, col = "red", xpd = TRUE) # par('usr') to legend insertion
        draw.circle(x = (PeakReg$peakSummit[i] - S)/steps, y = par('usr')[4]-(thd*0.1), radius = 50, border = "yellow", lwd=2, col="blue")
      }
#       
      if((PeakReg$peakSummit[length(PeakReg$peakStart)] - S)/steps<9000){
        mtext(side=3, line=-0.70, at=9900, adj=1, cex=0.5, "macs2-peaks", col = 'red')
      } 
    }
  }
#   
#   #draw Replication origins as circles
#   
  if(length(Ori_chr$chromStart) > 0){
    for(i in 1:length(Ori_chr$chromStart)){
     draw.circle(x = ((Ori_chr$chromStart[i]+Ori_chr$chromEnd[i])/2 - S)/steps, y = 0, radius = 60, border = "purple", lwd=2, col="yellow")
    }
  }
#   
#   #put chromosome name as we are plotting all chromosomes
#   
  if(SamplType=="BrdU" & DataType=="readCoverage"){
    title(main = paste(Pro_1, " ", SamplType, " Chromosome", gsub("[[:punct:]]*chr[[:punct:]]*", "", seqnames(Scerevisiae)[k])), col="gray", adj = 0, cex.main=1.5, line = 0, outer = TRUE)
  }
#   
#   #put the datatype, both readCoverage and normalised Ratio
  if(DataType=="readCoverage"){
    mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "read-Coverage")
  }
  if(DataType=="strandedRatio"){
    if(length(PeakReg$peakStart) == 0){
      mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
    } 
    if(length(PeakReg$peakStart) > 0){
      if((PeakReg$peakSummit[1] - S)/steps<1500){
        mtext(side=1, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
      } else {
        mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
      }
    } 
  }
#   
#   #put origin names
#   
  S1 <- Starts[seq(1, length(Starts), 3)]
  S2 <- Starts[seq(2, length(Starts), 3)]
  S3 <- Starts[seq(3, length(Starts), 3)]
#   
  if(length(Ori_chr$chromStart) > 0){
#     
    if(DataType=="readCoverage"){
#       
      Draw_name <- function(Ss, x){
#         
        colors_ars <- rep("gray1", length(Ori_chr$name))
        colors_ars[which(Ori_chr$stat=='early')] <- "red"
        colors_ars[which(Ori_chr$stat=='late')] <- "blue"
#         
        Ori_ticks <- round((c(Ori_chr$chromStart+Ori_chr$chromEnd)/2 - S)/steps) # establish the origins by the mean and dividing by the steps, as done all along the process
#         
        GetDistances <- function(x){ # define the function to get the distances ploted
          x <- sort(x)
          Distances <- c()
          for(i in 1:length(x)-1){
            Distances <- c(Distances, x[i + 1] - x[i])
          }
          return(Distances)
        }
#         
        Ori_Dists <- GetDistances(Ori_ticks) # distances for the origins, calculated by the mean and the division by the step
#         
        DistIndex <- which(Ori_Dists<500)
#         
        CloseOriIndices <- unique(sort(c(DistIndex, DistIndex+1)))
#         
        for(i in 1:length(Ss)){
          if(S == Ss[i]){
#             
            if(length(Ori_ticks)>0){
#               
              if(length(CloseOriIndices)==0){
                Ori_ticks_distal <- Ori_ticks
                Ori_name_distal <- Ori_chr$name
                color_distal <- colors_ars
              } else {
                Ori_ticks_distal <- Ori_ticks[-CloseOriIndices]
                Ori_name_distal <- Ori_chr$name[-CloseOriIndices]
                color_distal <- colors_ars[-CloseOriIndices]
              }
#               
              if(Ori_ticks_distal[1]>300){
                axis(side = 3, at = Ori_ticks_distal, labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                axis(side = 3, at = Ori_ticks_distal, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#                 # intToUtf8(9650) use arrows
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal)), at = Ori_ticks_distal, side = 3, line=0, col = "purple", cex = 1)
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal)), at = Ori_ticks_distal, side = 3, line=0.075, col = "yellow", cex = 0.65)
#                 
                text(x = Ori_ticks_distal, y = grconvertY(x, from = "ndc"), labels = Ori_name_distal, xpd = NA, srt = 0, col = color_distal, cex = 0.9 )
              } else {
                if(length(Ori_ticks_distal)>1){
                  axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#                   
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0, col = "purple", cex = 1)
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0.075, col = "yellow", cex = 0.65)
#                   
                  text(x = Ori_ticks_distal[-1], y = grconvertY(x, from = "ndc"), labels = Ori_name_distal[-1], xpd = NA, srt = 0, col = color_distal[-1], cex = 0.9 )
#                   ###
                  axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#                   
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
#                   
                  text(x = Ori_ticks_distal[1], y = grconvertY(x+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
                } else {
                  axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#                   
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
#                   
                  text(x = Ori_ticks_distal[1], y = grconvertY(x+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
                }
#                 
              }
              
            }
#             
            if(length(CloseOriIndices)>1){
              Consec_Oris <- split(CloseOriIndices, cumsum(c(1, diff(CloseOriIndices) != 1)))
#               
              if(length(Consec_Oris)>=1){
                for(j in 1:length(Consec_Oris)){
                  for(s in 0:(length(Consec_Oris[[j]])-1)){
#                     
                    Pos <- Ori_ticks[Consec_Oris[[j]]][s+1]
                    Nam <- Ori_chr$name[Consec_Oris[[j]]][s+1]
                    Col <- colors_ars[Consec_Oris[[j]]][s+1]
                    
                    tckL <- seq(0.05, 0.75, length.out = 5)
                    
                    axis(side = 3, at = Pos, labels = F, line = 0, tick = TRUE, tck=-(tckL[s+1]), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    axis(side = 3, at = Pos, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#                     
                    lineL <- seq(0, 6.20, length.out = 5)
#                     
                    mtext(rep(intToUtf8(9650), length(Pos)), at = Pos, side = 3, line=(lineL[s+1])+0, col = "purple", cex = 0.90)
                    mtext(rep(intToUtf8(9650), length(Pos)), at = Pos, side = 3, line=(lineL[s+1])+0.075, col = "yellow", cex = 0.585)
#                     
                    texL <- seq(0, 0.0725, length.out = 5)
#                     
                    text(x = Pos, y = grconvertY(texL[s+1]+x, from = "ndc"), labels = Nam, xpd = NA, srt = 0, col = Col, cex = 0.80 )
                  }
                }
              }
            }
#             
          }
        }
#         
      }
#       
      Draw_name(S1, 0.92)
#       
      Draw_name(S2, 0.615)
#       
      Draw_name(S3, 0.31)
    }
  }
#   
#   
#   
#   
  #plot x axis with by adding the chromosomal coordinates
  if(E>seqlengths(Scerevisiae)[k]){
    axis_ticks <-  seq(0, (seqlengths(Scerevisiae)[k] - S)/steps, ((seqlengths(Scerevisiae)[k] - S)/(steps*(seqlengths(Scerevisiae)[k] - S)/10000))) 
  } else {
    axis_ticks <-  seq(0, (E - S)/steps, ((E - S)/(steps*(Length_per_Row/10000))))
  }
  axis_ticks[1] <- 1
  axis(side = 1, at = axis_ticks, labels = rep(NA, 11), line = 0, tick = TRUE, tck=0.03, lwd.ticks = 1.5)
  axis(side = 1, at = axis_ticks, labels = rep(NA, 11), line = 0, tick = TRUE, tck=-0.03, lwd.ticks = 1.5)
  title(ylab='readDensity', col="gray", cex.lab=1.25, line = 0, outer = T)
#   
  if(DataType=="strandedRatio" & SamplType=="BrdU"){
    axis(side = 1, at = axis_ticks, labels = round((S/1000)+axis_ticks*(steps/1000)), line = 0, tick = F)
    title(xlab="Chromosomal Coordinates (Kbp)", col="gray", cex.lab=1.25, line = 0, outer = T)
  }
#   

  box("figure", col="forestgreen") 
  
}

  raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_global_profiles.pdf"), width = 9, height = 11, units = "in", res = 400)
  par(oma=c(2,2,2,2))
  for(k in 1:16){
    #     
    #define plot layout for global profiles
    LayOut.Dims <- function(x, y){
      xx <- c()
      y <- y
      for(i in 1:length(x)){
        xx <- c(xx, rep(x[i], y))
      }
      return(xx)
    }
    Plot.Nums <- c(1:6)
    Cols <- 18
    #     
    mat <- matrix(LayOut.Dims(Plot.Nums, Cols),
                  length(Plot.Nums),Cols,byrow=TRUE)
    #     
    vec <- c(1,4,7,10)
    new_mat <- matrix(0,nrow=length(Plot.Nums)+4,ncol=Cols)
    new_mat[-vec,] <- mat  
    #     
    fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(2,3,3,2,3,3,2,3,3,1)), ]
    #     
    layout(fin_mat, c(1,1), c(1,1), TRUE)
    #     
    #     
    k <- k
    #     
    #define plot intervals by chromosomes
    Length_per_Row <- 100000
    Plotting_Rows <- c(seq(0, seqlengths(Scerevisiae)[k]+Length_per_Row, Length_per_Row))
    Starts <- Plotting_Rows[-length(Plotting_Rows)]
    Ends <- Plotting_Rows[-1]
    #     
    #plot global profiles, all the profiles in all the chromosomes
    for(i in 1:length(Plotting_Rows[-length(Plotting_Rows)])){
      #       
      S <- Starts[i]
      E <- Ends[i]
      #       
      BrdU_Plot_function(CoverageFile = BrdU_coverage, peakFile = BrdU_Peaks, DataType="readCoverage", SamplType = "BrdU")
      BrdU_Plot_function(CoverageFile = BrdU_ratio, peakFile = BrdU_Peaks, DataType="strandedRatio", SamplType = "BrdU")
      #       
    }
  }
  dev.off() # close the document
  #   
  
  
  ## To plot just a region of a chromosome and not all the profiles:
Local_Profile <- function(CoverageFile, peakFile, DataType, SamplType, ChromCoords){ # We create a new function, the same as before but adding ChromCoords
#   # We define the chromosome names form the seqnames of s cerevisiae
  Coverage_chr <- CoverageFile[CoverageFile$chrom == seqnames(Scerevisiae)[k], ]
  All_Ori_chr <- All_Ori[All_Ori$chrom == seqnames(Scerevisiae)[k], ]
  Peaks2Plot_chr <- peakFile[peakFile$chrom == seqnames(Scerevisiae)[k], ]
#   # We indicate the regions included for the coverage, origins and peaks, the S and the E will be defined later
  Coverage <-  Coverage_chr[Coverage_chr$chromStart>=S & Coverage_chr$chromStart<=E, ]
  Ori_chr <- All_Ori_chr[All_Ori_chr$chromStart>=S & All_Ori_chr$chromStart<=E, ]
  PeakReg <- Peaks2Plot_chr[Peaks2Plot_chr$peakStart>=S & Peaks2Plot_chr$peakStart<=E, ]
#   # We take the scores of the two strands as coverage data
  CovWat <- Coverage$wat.score
  CovCri <- Coverage$cri.score
#   # We set the step as before done with the BrdU plot
  steps <- 10
#   
  if(stepSize == steps){
    CovWat <- CovWat
    CovCri <- CovCri
  }
  if(stepSize > steps){
    s <- stepSize/steps
    CovWat <- as.vector(sapply(CovWat, function (x) rep(x,s)))
    CovCri <- as.vector(sapply(CovCri, function (x) rep(x,s)))
  }
  if(stepSize < steps){
    s <- round(steps/stepSize)
    CovWat <- round(as.vector(tapply(CovWat, gl(length(CovWat)/s, s), mean)))
    CovCri <- round(as.vector(tapply(CovCri, gl(length(CovCri)/s, s), mean)))
  }
#   # This part is the same as before
#   
  scoreVals <- c(CovWat, CovCri)
#   
  Q <- quantile(scoreVals, probs=c(.01, .99), na.rm = FALSE)
  I <- IQR(scoreVals)
  up  <-  Q[2]+1.5*I # Upper Range  
  low <- Q[1]-1.5*I # Lower Range
  Scores <- scoreVals[which(scoreVals < up & scoreVals > low)]
  thd <- round(mean(Scores)+(12*sd(Scores)))
  Ylim_sc <- c(-thd, thd)
#   
  wcVals <- log2(CoverageFile$wat.score/CoverageFile$cri.score) # Same calculation by using the logarithm 
  wcVals[!is.finite(wcVals)] <- 0
#   
  Q <- quantile(wcVals, probs=c(.01, .99), na.rm = FALSE)
  I <- IQR(wcVals)
  up  <-  Q[2]+1.5*I # Upper Range  
  low <- Q[1]-1.5*I # Lower Range
  wc <- wcVals[which(wcVals < up & wcVals > low)]
  wcs <- round(mean(wc)+(3*sd(wc)))
  Ylim_wc <- c(-wcs, wcs)
#   
  yAxis_reads <- c(-0.3*thd, -0.6*thd, -0.9*thd, 0*thd, 0.3*thd, 0.6*thd, 0.9*thd)
  yAxis_wc <- c(-0.3*wcs, -0.6*wcs, -0.9*wcs, 0*wcs, 0.3*wcs, 0.6*wcs, 0.9*wcs)
#   
  WCrat <- log2(CovWat/CovCri)
  WCrat[!is.finite(WCrat)] <- 0
#   
  Cov.wat <- round(as.numeric(smooth.spline(1:length(CovWat), CovWat)$y), 3) # Again using the smooth.spline
  Cov.cri <- round(as.numeric(smooth.spline(1:length(CovCri), CovCri)$y), 3)
#   
  WC.rat <- round(as.numeric(smooth.spline(1:length(WCrat), WCrat)$y), 3)
#   
  Length_per_Row <- E - S
#   
  if(length(Cov.wat) < (Length_per_Row/steps - 1)){
    Cov.wat <- c(Cov.wat, rep(NA, (Length_per_Row/steps - 1)-length(Cov.wat)) )
  }
#   
  if(length(Cov.cri) < (Length_per_Row/steps - 1)){
    Cov.cri <- c(Cov.cri, rep(NA, (Length_per_Row/steps - 1)-length(Cov.cri)) )
  }
#   
  if(length(WC.rat) < (Length_per_Row/steps - 1)){
    WC.rat <- c(WC.rat, rep(NA, (Length_per_Row/steps - 1)-length(WC.rat)) )
  }
#   
  #plot
  par(mar = c(0,0,0,0))
  suppressWarnings(
    plot(Cov.wat, type='h', ylim=Ylim_sc, col =  rgb(100,0,0,alpha=180, maxColorValue=255), ylab=' ', 
         xlab=' ', xaxt='n', yaxt='n', lwd=0.07, bty = 'n',  cex.lab=1, las = 2, xaxs='i')
  ) 
  lines(Cov.cri*(-1), type='h', col =  rgb(0,100,0,alpha=180, maxColorValue=255), lwd=0.07)
#   
  segments(x0 = c(Coverage$chromStart - S)/steps, x1 = c(Coverage$chromEnd - S)/steps, y0 = 0, y1= 0, lwd = 0.5)
#   
  #plot left y axis
  axis(side = 2, at = yAxis_reads, labels = round(yAxis_reads), line = 0, tick = TRUE, lwd.ticks = 1.5, las = 2, cex.axis = 0.8)
  abline(h=yAxis_reads, lwd=0.05, col =  rgb(112,128,144,alpha=225, maxColorValue=255))
#   
  #draw the peaks
  if(DataType=="strandedRatio"){
    if(length(PeakReg$peakStart) > 0){
      for(i in 1:length(PeakReg$peakStart)){
        segments(x0 = c(PeakReg$peakStart[i] - S)/steps, x1 = c(PeakReg$peakEnd[i] - S)/steps, 
                 y0 = par('usr')[4]-(thd*0.1), y1= par('usr')[4]-(thd*0.1), lwd = 4, col = "red", xpd = TRUE)
        draw.circle(x = (PeakReg$peakSummit[i] - S)/steps, y = par('usr')[4]-(thd*0.1), radius = 50, border = "yellow", lwd=2, col="blue")
      }
#       
      if((PeakReg$peakSummit[length(PeakReg$peakStart)] - S)/steps<9000){
        mtext(side=3, line=-0.70, at=9900, adj=1, cex=0.5, "macs2-peaks", col = 'red')
      } 
    }
  }
#   
  #draw Replication origins
#   
  if(length(Ori_chr$chromStart) > 0){
    for(i in 1:length(Ori_chr$chromStart)){
      draw.circle(x = ((Ori_chr$chromStart[i]+Ori_chr$chromEnd[i])/2 - S)/steps, y = 0, radius = 60, border = "purple", lwd=2, col="yellow")
    }
  }
#   
  #put chromosome name
#   
  if(SamplType=="BrdU" & DataType=="readCoverage"){
    title(main = paste0(Pro_1, "_", SamplType, "_Chromosome", gsub("[[:punct:]]*chr[[:punct:]]*", "", seqnames(Scerevisiae)[k]), ":", S, "-", E), col="gray", adj = 0, cex.main=1.5, line = 0, outer = TRUE)
  }
#   
#   
  # #put sample name
   if(DataType=="readCoverage"){
     mtext(side=3, line=0.75, at=-50, adj=1, cex=1, SamplType)
   }
#   
  #put the datatype
  if(DataType=="readCoverage"){
    mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "read-Coverage")
  }
  if(DataType=="strandedRatio"){
    if(length(PeakReg$peakStart) == 0){
     mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
    } 
    if(length(PeakReg$peakStart) > 0){
      if((PeakReg$peakSummit[1] - S)/steps<1500){
        mtext(side=1, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
      } else {
        mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
      }
    } 
  }
#   
  #put origin names
#   
  if(length(Ori_chr$chromStart) > 0){
#     
    colors_ars <- rep("gray1", length(Ori_chr$name))
    colors_ars[which(Ori_chr$stat=='early')] <- "red"
    colors_ars[which(Ori_chr$stat=='late')] <- "blue"
#     
    Ori_ticks <- round((c(Ori_chr$chromStart+Ori_chr$chromEnd)/2 - S)/steps)
#     
    GetDistances <- function(x){
      x <- sort(x)
      Distances <- c()
      for(i in 1:length(x)-1){
        Distances <- c(Distances, x[i + 1] - x[i])
      }
      return(Distances)
    }
#     
    Ori_Dists <- GetDistances(Ori_ticks)
#     
    DistIndex <- which(Ori_Dists<500)
#     
    CloseOriIndices <- unique(sort(c(DistIndex, DistIndex+1)))
#     
#     
    if(DataType=="readCoverage"){
#       
      if(length(Ori_ticks)>0){
#         
        if(length(CloseOriIndices)==0){
          Ori_ticks_distal <- Ori_ticks
          Ori_name_distal <- Ori_chr$name
          color_distal <- colors_ars
        } else {
          Ori_ticks_distal <- Ori_ticks[-CloseOriIndices]
          Ori_name_distal <- Ori_chr$name[-CloseOriIndices]
          color_distal <- colors_ars[-CloseOriIndices]
        }
#         
        if(Ori_ticks_distal[1]>300){
          axis(side = 3, at = Ori_ticks_distal, labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
          axis(side = 3, at = Ori_ticks_distal, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#           
          mtext(rep(intToUtf8(9650), length(Ori_ticks_distal)), at = Ori_ticks_distal, side = 3, line=0, col = "purple", cex = 1)
          mtext(rep(intToUtf8(9650), length(Ori_ticks_distal)), at = Ori_ticks_distal, side = 3, line=0.075, col = "yellow", cex = 0.65)
#           
          text(x = Ori_ticks_distal, y = grconvertY(0.715, from = "ndc"), labels = Ori_name_distal, xpd = NA, srt = 0, col = color_distal, cex = 0.9 )
        } else {
          if(length(Ori_ticks_distal)>1){
            axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
            axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#             
            mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0, col = "purple", cex = 1)
            mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0.075, col = "yellow", cex = 0.65)
#             
            text(x = Ori_ticks_distal[-1], y = grconvertY(0.715, from = "ndc"), labels = Ori_name_distal[-1], xpd = NA, srt = 0, col = color_distal[-1], cex = 0.9 )
#             ###
            axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
            axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#             
            mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
            mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
#             
            text(x = Ori_ticks_distal[1], y = grconvertY(0.715+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
           } else {
            axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
            axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#             
            mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
            mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
#             
            text(x = Ori_ticks_distal[1], y = grconvertY(0.715+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
          }
#           
        }
#         
      }
#       
      if(length(CloseOriIndices)>1){
        Consec_Oris <- split(CloseOriIndices, cumsum(c(1, diff(CloseOriIndices) != 1)))
#         
        if(length(Consec_Oris)>=1){
          for(j in 1:length(Consec_Oris)){
            for(s in 0:(length(Consec_Oris[[j]])-1)){
#               
              Pos <- Ori_ticks[Consec_Oris[[j]]][s+1]
              Nam <- Ori_chr$name[Consec_Oris[[j]]][s+1]
              Col <- colors_ars[Consec_Oris[[j]]][s+1]
#               
              tckL <- seq(0.05, 0.75, length.out = 5)
#               
              axis(side = 3, at = Pos, labels = F, line = 0, tick = TRUE, tck=-(tckL[s+1]), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              axis(side = 3, at = Pos, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#               
              lineL <- seq(0, 6.20, length.out = 5)
#               
              mtext(rep(intToUtf8(9650), length(Pos)), at = Pos, side = 3, line=(lineL[s+1])+0, col = "purple", cex = 0.90)
              mtext(rep(intToUtf8(9650), length(Pos)), at = Pos, side = 3, line=(lineL[s+1])+0.075, col = "yellow", cex = 0.585)
#               
              texL <- seq(0, 0.0725, length.out = 5)
#               
              text(x = Pos, y = grconvertY(texL[s+1]+0.715, from = "ndc"), labels = Nam, xpd = NA, srt = 0, col = Col, cex = 0.80 )
            }
          }
        }
      }
#       
#       
#       
#       
    } else {
#       
      if(DataType=="strandedRatio"){
        axis(side = 3, at = Ori_ticks, labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
        axis(side = 3, at = Ori_ticks, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
#         
        mtext(rep(intToUtf8(9650), length(Ori_ticks)), at = Ori_ticks, side = 3, line=0, col = "purple", cex = 1)
        mtext(rep(intToUtf8(9650), length(Ori_ticks)), at = Ori_ticks, side = 3, line=0.075, col = "yellow", cex = 0.65)
      }
    }
#     
  }
#   
#   
  #plot x axis
  if(E>seqlengths(Scerevisiae)[k]){
    axis_ticks <-  seq(0, (seqlengths(Scerevisiae)[k] - S)/steps, ((seqlengths(Scerevisiae)[k] - S)/(steps*(seqlengths(Scerevisiae)[k] - S)/10000))) 
  } else {
    axis_ticks <-  seq(0, (E - S)/steps, ((E - S)/(steps*(Length_per_Row/10000))))
  }
  axis_ticks[1] <- 1
  axis(side = 1, at = axis_ticks, labels = rep(NA, 11), line = 0, tick = TRUE, tck=0.03, lwd.ticks = 1.5)
  axis(side = 1, at = axis_ticks, labels = rep(NA, 11), line = 0, tick = TRUE, tck=-0.03, lwd.ticks = 1.5)
  title(ylab='readDensity', col="gray", cex.lab=1.25, line = 0, outer = T)
#   
  if(DataType=="strandedRatio"){
    axis(side = 1, at = axis_ticks, labels = round((S/1000)+axis_ticks*(steps/1000)), line = 0, tick = F)
    title(xlab="Chromosomal Coordinates (Kbp)", col="gray", cex.lab=1.25, line = 0, outer = T)
  }
#   
 
  box("figure", col="forestgreen") 
}
# 
ChromCoords = "4:50000-150000" # We establish manually the chromosome coordinates

  k <- as.numeric(unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[1]]) # Indicate the first number in the chromcoords as the number of k to establish it as the chromosome to represent
  SE <- unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[2]] # We select the part after the : of the chromcoords
  S <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[1]]) # In this part, the strating point of the coordinates is the S, the beginning of the region represented
  E <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[2]])  # The part after the - represent the end of the region to represent
  # We name the file by using the coordinates indicated (k, S and E)
  raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "Chrom", as.roman(k), "_",  as.character(S), "-", as.character(E), ".pdf"), width = 9, height = 11, units = "in", res = 400)
  #   Establish the parameters to plotting
  par(oma=c(2,2,2,2))
  LayOut.Dims <- function(x, y){ # Function to indicate the design dimensions and disposition of the plots
    xx <- c()
    y <- y
    for(i in 1:length(x)){
      xx <- c(xx, rep(x[i], y))
    }
    return(xx)
  }
  Plot.Nums <- c(1:2); Cols <- 12
  mat <- matrix(LayOut.Dims(Plot.Nums, Cols), length(Plot.Nums),Cols,byrow=TRUE)
  vec <- c(1,3,5)
  new_mat <- matrix(0,nrow=length(Plot.Nums)+3,ncol=Cols)
  new_mat[-vec,] <- mat  
  fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(5,3,1,3,5)), ]
  layout(fin_mat, c(1,1), c(1,1), TRUE)
  
  Local_Profile(CoverageFile = BrdU_coverage, peakFile = BrdU_Peaks, DataType="readCoverage", SamplType = "BrdU", ChromCoords) # We plot both read coverage and stranded ratio
  Local_Profile(CoverageFile = BrdU_ratio, peakFile = BrdU_Peaks, DataType="strandedRatio", SamplType = "BrdU", ChromCoords)
  
  dev.off()
  
# Different conditions whether we want to plot global or local profiles:
  
 if(PlotIPprofile == "Global"){
#   
   raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_global_profiles.pdf"), width = 9, height = 11, units = "in", res = 400)
   par(oma=c(2,2,2,2))
   for(k in 1:16){
#     
     #define plot layout for global profiles
     LayOut.Dims <- function(x, y){
       xx <- c()
       y <- y
       for(i in 1:length(x)){
         xx <- c(xx, rep(x[i], y))
       }
       return(xx)
     }
     Plot.Nums <- c(1:6)
     Cols <- 18
#     
     mat <- matrix(LayOut.Dims(Plot.Nums, Cols),
                   length(Plot.Nums),Cols,byrow=TRUE)
#     
     vec <- c(1,4,7,10)
     new_mat <- matrix(0,nrow=length(Plot.Nums)+4,ncol=Cols)
     new_mat[-vec,] <- mat  
#     
     fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(2,3,3,2,3,3,2,3,3,1)), ]
#     
     layout(fin_mat, c(1,1), c(1,1), TRUE)
     #layout.show(n)
#     
#     
     k <- k
#     
     #define plot intervals by chromosomes
     Length_per_Row <- 100000
     Plotting_Rows <- c(seq(0, seqlengths(Scerevisiae)[k]+Length_per_Row, Length_per_Row))
     Starts <- Plotting_Rows[-length(Plotting_Rows)]
     Ends <- Plotting_Rows[-1]
#     
     #plot global profiles
     for(i in 1:length(Plotting_Rows[-length(Plotting_Rows)])){
#       
       S <- Starts[i]
       E <- Ends[i]
#       
       BrdU_Plot_function(CoverageFile = BrdU_coverage, peakFile = BrdU_Peaks, DataType="readCoverage", SamplType)
       BrdU_Plot_function(CoverageFile = BrdU_ratio, peakFile = BrdU_Peaks, DataType="strandedRatio", SamplType)
#       
     }
   }
   dev.off()
#   
 }
# 
 if(PlotIPprofile == "Local"){
#   
   k <- as.numeric(unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[1]])
   SE <- unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[2]]
   S <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[1]])
   E <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[2]])
   
   raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "Chrom", as.roman(k), "_",  as.character(S), "-", as.character(E), ".pdf"), width = 9, height = 11, units = "in", res = 400)
#   
   par(oma=c(2,2,2,2))
   LayOut.Dims <- function(x, y){
     xx <- c()
     y <- y
     for(i in 1:length(x)){
       xx <- c(xx, rep(x[i], y))
     }
     return(xx)
   }
   Plot.Nums <- c(1:2); Cols <- 12
   mat <- matrix(LayOut.Dims(Plot.Nums, Cols), length(Plot.Nums),Cols,byrow=TRUE)
   vec <- c(1,3,5)
   new_mat <- matrix(0,nrow=length(Plot.Nums)+3,ncol=Cols)
   new_mat[-vec,] <- mat  
   fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(5,3,1,3,5)), ]
   layout(fin_mat, c(1,1), c(1,1), TRUE)
   
   Local_Profile(CoverageFile = BrdU_coverage, peakFile = BrdU_Peaks, DataType="readCoverage", SamplType, ChromCoords)
   Local_Profile(CoverageFile = BrdU_ratio, peakFile = BrdU_Peaks, DataType="strandedRatio", SamplType, ChromCoords)
   
   dev.off()
   
 }
 
 if(PlotIPprofile == "Both"){
   
   raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_global_profiles.pdf"), width = 9, height = 11, units = "in", res = 400)
   par(oma=c(2,2,2,2))
   for(k in 1:16){
     
     #define plot layout for global profiles
     LayOut.Dims <- function(x, y){
       xx <- c()
       y <- y
       for(i in 1:length(x)){
         xx <- c(xx, rep(x[i], y))
       }
       return(xx)
     }
     Plot.Nums <- c(1:6)
     Cols <- 18
     
     mat <- matrix(LayOut.Dims(Plot.Nums, Cols),
                   length(Plot.Nums),Cols,byrow=TRUE)
#     
     vec <- c(1,4,7,10)
     new_mat <- matrix(0,nrow=length(Plot.Nums)+4,ncol=Cols)
     new_mat[-vec,] <- mat  
#     
     fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(2,3,3,2,3,3,2,3,3,1)), ]
#     
     layout(fin_mat, c(1,1), c(1,1), TRUE)
     #layout.show(n)
     
     
     k <- k
     
     #define plot intervals by chromosomes
     Length_per_Row <- 100000
     Plotting_Rows <- c(seq(0, seqlengths(Scerevisiae)[k]+Length_per_Row, Length_per_Row))
     Starts <- Plotting_Rows[-length(Plotting_Rows)]
     Ends <- Plotting_Rows[-1]
#     
     #plot global profiles
     for(i in 1:length(Plotting_Rows[-length(Plotting_Rows)])){
       
       S <- Starts[i]
       E <- Ends[i]
       
       BrdU_Plot_function(CoverageFile = BrdU_coverage, peakFile = BrdU_Peaks, DataType="readCoverage", SamplType)
       BrdU_Plot_function(CoverageFile = BrdU_ratio, peakFile = BrdU_Peaks, DataType="strandedRatio", SamplType)
       #############################################################################
       ## Coverage_chr <- CoverageFile[CoverageFile$chrom == seqnames(Scerevisiae)[k], ] PROBLEM WHEN LOOKING FOR [k]
     }
   }
   dev.off()
#   
   # plot local profile
   
   k <- as.numeric(unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[1]])
   SE <- unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[2]]
   S <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[1]])
   E <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[2]])
#   
   raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "Chrom", as.roman(k), "_",  S, "-", E, ".pdf"), width = 9, height = 11, units = "in", res = 400)
#   
   par(oma=c(2,2,2,2))
   LayOut.Dims <- function(x, y){
     xx <- c()
     y <- y
     for(i in 1:length(x)){
       xx <- c(xx, rep(x[i], y))
     }
     return(xx)
   }
   Plot.Nums <- c(1:2); Cols <- 12
   mat <- matrix(LayOut.Dims(Plot.Nums, Cols), length(Plot.Nums),Cols,byrow=TRUE)
   vec <- c(1,3,5)
   new_mat <- matrix(0,nrow=length(Plot.Nums)+3,ncol=Cols)
   new_mat[-vec,] <- mat  
   fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(5,3,1,3,5)), ]
   layout(fin_mat, c(1,1), c(1,1), TRUE)
   
   Local_Profile(CoverageFile = BrdU_coverage, peakFile = BrdU_Peaks, DataType="readCoverage", SamplType, ChromCoords)
   Local_Profile(CoverageFile = BrdU_ratio, peakFile = BrdU_Peaks, DataType="strandedRatio", SamplType, ChromCoords)
   
   dev.off()
   
 }
 
 if(PlotIPprofile == "None"){
   
   message("No IP profiles plotted...")
#   
 }
# 

  # AVERAGE ENRICHMENT
# plot average enrichment profile around the replication origins

peakFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", "*", ".bed"))

LeftLim <- as.numeric(read.table(peakFiles[7], row.names = 1)["Left", ]) # Number 7 is synthesis bed, we obtain each value from the file
RightLim <- as.numeric(read.table(peakFiles[7], row.names = 1)["Right", ])
LeadingAverage <-  read.table(peakFiles[7], row.names = 1)["LeadSynthesis", ]
LaggingAverage <-  read.table(peakFiles[7], row.names = 1)["LaggSynthesis", ]
binSize <- as.numeric(read.table(peakFiles[7], row.names = 1)["bin", ])
stepSize <- as.numeric(read.table(peakFiles[7], row.names = 1)["slide", ])
AveragingWindow <- as.numeric(read.table(peakFiles[7], row.names = 1)["AveragingWindow", ])

# To obtain the average enrichment profiles we define a new function
Average_Enrichment <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = TRUE){
  
  watsonRatio <- read.table(watsonRatio, header = T)
  crickRatio <- read.table(crickRatio, header = T)
  PeakList <- PeakList
  
  Window = AveragingWindow
  
  if(Normalise2Input == TRUE){ # We indicate different V values depending on the normalization, this vlues will be later columns in IP files
    V <- 7 # The function is defined with normalisetoinput as true, so the number will be 7
  } else {
    V <- 5
  } 
  
  ###
  
  PeakList$AvBstart <- PeakList$mid - Window # The start and end of the averaging profiles will be defined by the window
  PeakList$AvBend <- PeakList$mid + Window
  
  chrS <- paste0("chr", as.roman(1:16)) # We name the chromosomes from 1 to 16 as roman numbers (I to XVI)
  
  #watson
  IP_R <- watsonRatio
  IP_W <- NULL
  for(i in 1:length(chrS)){
    
    i <- i
    
    OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ] # PeakList containing chromsomes as the ones established from I to XVI
    IP_ROs <- IP_R[IP_R$chrom == chrS[i], ] # Same to the watsonratio
    
    if(length(OriList_ROs$chrom)==0) next 
    
    IP_C <- NULL
    for(y in 1:length(OriList_ROs$chrom)){
      
      y <- y
      
      IP_Z <- IP_ROs[IP_ROs$chromStart>=OriList_ROs$AvBstart[y] & IP_ROs$chromStart<=OriList_ROs$AvBend[y], ] # The intervals based on the peaklist chromosomes
      
      if(length(IP_Z[,V])==round(2*Window/stepSize)){ # The V varies if normalizing to input or not
        Ratios <- IP_Z[,V] # Ratios are directly taken from IP_Z if the length of it is the same as two times the window between 10 (stepsize)
      } 
      
      if(length(IP_Z[,V]) < round(2*Window/stepSize)){ # Two different conditions in case IP_Z length is less than two times window/stepsize
        if(length(1:(length(IP_Z[,V])/2)) < round(2*Window/stepSize)/2){ # Total length of tee IP_Z from 1 is less than the 2*window/stepSize
          Ratios <- c(rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)), IP_Z[,V] ) # Ratios are replicates of the vector 0, 2*window/stepsize minus the lenght of IP_Z start of chromsomes, and the IP_Z column V
        }
        if(length((length(IP_Z[,V])/2+1):length(IP_Z[,V])) < round(2*Window/stepSize)/2){ # Total length of tee IP_Z from half IP_Z plus 1 to IP_Z column 7 is less than the 2*window/stepSize
          Ratios <- c(IP_Z[,V], rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)) ) # Ratios are the combination of the IP_Z and the same as in the other case: replicates of the vector 0, 2*window/stepsize minus the lenght of IP_Z start of chromsomes
        }
      } 
      
      if(length(IP_Z[,V]) > round(2*Window/stepSize)){ # Last case is if the IP_Z column V is bigger than 2*window/stepsize
        Ratios <- IP_Z[,V][-c((round(2*Window/stepSize)+1):length(IP_Z[,V]))] # Ratios are the opposite combination 
      }
      
      Rat <- matrix(Ratios, ncol = 1) # Ratios as matrix
      
      IP_C <- cbind(IP_C, Rat)
    }
    IP_W <- cbind(IP_W, IP_C) # The Watson ratios
  }
  IP_W <- as.data.frame(IP_W)
  colnames(IP_W) <- c(1:length(PeakList$chrom))
  
  #crick
  IP_R <- crickRatio # Same for the crick strand
  IP_Cr <- NULL
  for(i in 1:length(chrS)){
    
    OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
    IP_ROs <- IP_R[IP_R$chrom == chrS[i], ]
    
    if(length(OriList_ROs$chrom)==0) next 
    
    IP_C <- NULL
    for(y in 1:length(OriList_ROs$chrom)){
      
      IP_Z <- IP_ROs[IP_ROs$chromStart>=OriList_ROs$AvBstart[y] & IP_ROs$chromStart<=OriList_ROs$AvBend[y], ]
      
      if(length(IP_Z[,V])==round(2*Window/stepSize)){
        Ratios <- IP_Z[,V]
      } 
      
      if(length(IP_Z[,V]) < round(2*Window/stepSize)){
        if(length(1:(length(IP_Z[,V])/2)) < round(2*Window/stepSize)/2){
          Ratios <- c(rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)), IP_Z[,V] )
        }
        if(length((length(IP_Z[,V])/2+1):length(IP_Z[,V])) < round(2*Window/stepSize)/2){
          Ratios <- c(IP_Z[,V], rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)) )
        }
      } 
      
      if(length(IP_Z[,V]) > round(2*Window/stepSize)){
        Ratios <- IP_Z[,V][-c((round(2*Window/stepSize)+1):length(IP_Z[,V]))]
      }
      
      Rat <- matrix(Ratios, ncol = 1)
      
      IP_C <- cbind(IP_C, Rat)
    }
    IP_Cr <- cbind(IP_Cr, IP_C)
  }
  IP_Cr <- as.data.frame(IP_Cr)
  colnames(IP_Cr) <- c(1:length(PeakList$chrom))
  
  ###
  rowStat <- function(DF){ # New function to establish the statistic parameters we want to obtain as the median, q.25, q.75, mean, sd and standard error
    
    Quantiles <- apply(DF, 1, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    Ses <- apply(DF, 1, std.error)
    Means <- apply(DF, 1, mean)
    Sds <- apply(DF, 1, sd)
    
    return(list(q.25 = Quantiles["25%",],
                Median = Quantiles["50%",],
                q.75 = Quantiles["75%",],
                Std.err = Ses,
                Mean = Means,
                Sd = Sds))
  }
  # To get the average enrichment for watson and crick we get the statistical parameters
  watsonAverage <- cbind.data.frame(watson.q25 = rowStat(IP_W)$q.25, watson.median = rowStat(IP_W)$Median, watson.q75 = rowStat(IP_W)$q.75, watson.se = rowStat(IP_W)$Std.err, watson.mean = rowStat(IP_W)$Mean, watson.sd = rowStat(IP_W)$Sd)
  crickAverage <- cbind.data.frame(crick.q25 = rowStat(IP_Cr)$q.25, crick.median = rowStat(IP_Cr)$Median, crick.q75 = rowStat(IP_Cr)$q.75, crick.se = rowStat(IP_Cr)$Std.err, crick.mean = rowStat(IP_Cr)$Mean, crick.sd = rowStat(IP_Cr)$Sd)
  
  TwoStrands <- cbind.data.frame(watsonAverage, crickAverage)
  
  rownames(TwoStrands) <- paste0("bin", 1:length(TwoStrands[,1]))
  
  return(TwoStrands)
}
# We get the average enrichment for early and late origins in the three samples 
IP_AvE_Early <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = TRUE)
Input_AvE_Early <- Average_Enrichment(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = E_Ori, Normalise2Input = TRUE)
SUP_AvE_Early <- Average_Enrichment(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = E_Ori, Normalise2Input = TRUE)

IP_AvE_Late <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = TRUE)
Input_AvE_Late <- Average_Enrichment(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = L_Ori, Normalise2Input = TRUE)
SUP_AvE_Late <- Average_Enrichment(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = L_Ori, Normalise2Input = TRUE)

# W/C AVERAGE BY log2
# We define a new function to get the average enrichment of the watson strand over the crick using the same arguments
WatsonOverCrick_Average <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = TRUE){
  # Same as before
  Watson <- read.table(watsonRatio, header = T)
  Crick <- read.table(crickRatio, header = T)
  PeakList <- PeakList
  
  Window = AveragingWindow
  
  if(Normalise2Input == TRUE){ # Same V value as in the previous defined function
    V <- 7
  } else {
    V <- 5
  } 
  
  ###
  
  # As before:
  PeakList$AvBstart <- PeakList$mid - Window
  PeakList$AvBend <- PeakList$mid + Window

  chrS <- paste0("chr", as.roman(1:16))
  
  IP_T <- NULL
  for(i in 1:length(chrS)){
    
    i <- i
    # IN this case we put W and C together in the same loop and not separatedly as before
    OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
    Crick_ROs <- Crick[Crick$chrom == chrS[i], ]
    Watson_ROs <- Watson[Watson$chrom == chrS[i], ]
    
    if(length(OriList_ROs$chrom)==0) next 
    
    IP_C <- NULL
    
    for(y in 1:length(OriList_ROs$chrom)){
      
      y <- y
      
      Crick_S <- Crick_ROs[Crick_ROs$chromStart>=OriList_ROs$AvBstart[y] & Crick_ROs$chromStart<=OriList_ROs$AvBend[y], ] # Crick replication origins by using origin list
      Watson_S <- Watson_ROs[Watson_ROs$chromStart>=OriList_ROs$AvBstart[y] & Watson_ROs$chromStart<=OriList_ROs$AvBend[y], ]
      
      IP_Z <- log2(Watson_S[,V] / Crick_S[,V]); IP_Z[!is.finite(IP_Z)] <- 0 # The w over c got by the logarithm of the division of w by c
      # Conditions as before depending on the IP_Z value compared to two times the window divided by the stepsize
      if(length(IP_Z)==round(2*Window/stepSize)){
        Ratios <- IP_Z
      } 
      
      if(length(IP_Z) < round(2*Window/stepSize)){
        if(length(1:(length(IP_Z)/2)) < round(2*Window/stepSize)/2){
          Ratios <- c(rep(0, (round(2*Window/stepSize))-length(Crick_S$chromStart)), IP_Z )
        }
        if(length((length(IP_Z)/2+1):length(IP_Z)) < round(2*Window/stepSize)/2){
          Ratios <- c(IP_Z, rep(0, (round(2*Window/stepSize))-length(Crick_S$chromStart)) )
        }
      } 
      
      if(length(IP_Z) > round(2*Window/stepSize)){
        Ratios <- IP_Z[-c((round(2*Window/stepSize)+1):length(IP_Z))]
      }
      
      Rat <- matrix(Ratios, ncol = 1)
      
      IP_C <- cbind(IP_C, Rat)
      
    }
    
    IP_T <- cbind(IP_T, IP_C)
    IP_T <- as.data.frame(IP_T) # Final dataframe is the ratios
  }
  colnames(IP_T) <- c(1:length(PeakList$chrom))
  
  ##extract bin stats by defining a new function containing meduan, mean, sd, standard error, q.25 and q.75
  rowStat <- function(DF){
    
    
    Quantiles <- apply(DF, 1, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    Ses <- apply(DF, 1, std.error)
    Means <- apply(DF, 1, mean)
    Sds <- apply(DF, 1, sd)
    
    return(list(q.25 = Quantiles["25%",],
                Median = Quantiles["50%",],
                q.75 = Quantiles["75%",],
                Std.err = Ses,
                Mean = Means,
                Sd = Sds))
  }
  # The W over C average enrichment is obtained by getting from the IP_T the statistical parameters
  watsonOvercrickAvg <- cbind.data.frame(q25 = rowStat(IP_T)$q.25, median = rowStat(IP_T)$Median, q75 = rowStat(IP_T)$q.75, se = rowStat(IP_T)$Std.err, mean = rowStat(IP_T)$Mean, sd = rowStat(IP_T)$Sd)
  
  rownames(watsonOvercrickAvg) <- paste0("bin", 1:length(watsonOvercrickAvg[,1]))
  
  return(watsonOvercrickAvg)
  
}
##
# We obtain the Watson over Crick average enrichment of all the three samples and in early and late origins
IP_WoC_Early <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = TRUE)
Input_WoC_Early <- WatsonOverCrick_Average(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = E_Ori, Normalise2Input = TRUE)
SUP_WoC_Early <- WatsonOverCrick_Average(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = E_Ori, Normalise2Input = TRUE)

IP_WoC_Late <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = TRUE)
Input_WoC_Late <- WatsonOverCrick_Average(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = L_Ori, Normalise2Input = TRUE)
SUP_WoC_Late <- WatsonOverCrick_Average(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = L_Ori, Normalise2Input = TRUE)

# Bias at individual peaks, not average between all the peaks
Bias_at_individual_peaks <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = TRUE){ # Same arguments in the new function
  
  watsonRatio <- read.table(watsonRatio, header = T)
  crickRatio <- read.table(crickRatio, header = T)
  PeakList <- PeakList
  
  Window = AveragingWindow
  
  PeakList$AvBstart <- PeakList$mid - Window
  PeakList$AvBend <- PeakList$mid + Window
  
  chrS <- paste0("chr", as.roman(1:16))
  
  IP_T <- NULL
  
  for(i in 1:length(chrS)){
    
    i <- i
    
    IP_ROs_pos <- watsonRatio[watsonRatio$chrom == chrS[i], ] # watson as positive
    IP_ROs_neg <- crickRatio[crickRatio$chrom == chrS[i], ] # Crick as negative 
    Chr_Peaks <- PeakList[PeakList$chrom == chrS[i], ]
    
    if(length(Chr_Peaks$chrom)==0) next 
    
    ###
    
    IP_C <- NULL
    
    for(y in 1:length(Chr_Peaks$chrom)){
      
      y <- y
      # All time differenciating between w and c as positive and negative
      IP_Z_pos <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$AvBstart[y] & IP_ROs_pos$chromEnd<=Chr_Peaks$AvBend[y], ]
      IP_Z_neg <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$AvBstart[y] & IP_ROs_neg$chromEnd<=Chr_Peaks$AvBend[y], ]
      # Need to add new concept: right and left as we are dividing the ROs on four sections: WL WR CR CL
      IP_Z_CL <- sum(IP_Z_pos[,5][1:(Window/stepSize)], na.rm = TRUE)
      IP_Z_CR <- sum(IP_Z_pos[,5][(Window/stepSize):((2*Window)/stepSize - 1)], na.rm = TRUE)
      IP_Z_WL <- sum(IP_Z_neg[,5][1:(Window/stepSize)], na.rm = TRUE)
      IP_Z_WR <- sum(IP_Z_neg[,5][(Window/stepSize):((2*Window)/stepSize - 1)], na.rm = TRUE)
      
      In_Z_CL <- sum(IP_Z_pos[,6][1:(Window/stepSize)], na.rm = TRUE)
      In_Z_CR <- sum(IP_Z_pos[,6][(Window/stepSize):((2*Window)/stepSize - 1)], na.rm = TRUE)
      In_Z_WL <- sum(IP_Z_neg[,6][1:(Window/stepSize)], na.rm = TRUE)
      In_Z_WR <- sum(IP_Z_neg[,6][(Window/stepSize):((2*Window)/stepSize - 1)], na.rm = TRUE)
      # Leading and lagging strands are formed by W and C: 
      IP_Z_CL_WR <- IP_Z_CL + IP_Z_WR; IP_Z_CL_WR[is.na(IP_Z_CL_WR)]<-0 # Lagging strand
      IP_Z_WL_CR <- IP_Z_WL + IP_Z_CR; IP_Z_WL_CR[is.na(IP_Z_WL_CR)]<-0 # Leading strand
      # Binomial test to stablish the significance of the peaks 
      IP_Z_Pvalue <- binom.test(round(c(IP_Z_CL_WR, IP_Z_WL_CR)))$p.value
      
      In_Z_CL_WR <- In_Z_CL + In_Z_WR; In_Z_CL_WR[is.na(In_Z_CL_WR)]<-0
      In_Z_WL_CR <- In_Z_WL + In_Z_CR; In_Z_WL_CR[is.na(In_Z_WL_CR)]<-0
      # Ratio lagging leading by division:
      In_Z_LaLe <- In_Z_CL_WR/In_Z_WL_CR; In_Z_LaLe[!is.finite(In_Z_LaLe)] <- 0
      IP_Z_LaLe <- IP_Z_CL_WR/IP_Z_WL_CR; IP_Z_LaLe[!is.finite(IP_Z_LaLe)] <- 0
      
      
      if(Normalise2Input=="NO"){
        IP_Z_LagLead <- log2(IP_Z_LaLe); IP_Z_LagLead[!is.finite(IP_Z_LagLead)] <- 0 # If not normalising: the logarithm of the IP
      } else {
        IP_Z_LagLead <- log2(IP_Z_LaLe/In_Z_LaLe); IP_Z_LagLead[!is.finite(IP_Z_LagLead)] <- 0 # If normalising, the logarithm of the division of the IP by the input
      }
      
      IP_Z <- cbind.data.frame(round(IP_Z_CL_WR), round(IP_Z_WL_CR), IP_Z_Pvalue, IP_Z_LagLead) # All the data, for leading, lagging and statistical values
      IP_C <- rbind.data.frame(IP_C, IP_Z)
    }
    IP_C <- cbind.data.frame(Chr_Peaks$chrom, Chr_Peaks$name, IP_C)
    colnames(IP_C) <- c('chrom', 'name', "Lagg.sum", "Lead.sum", "p_value", "Bias") # Adding column names
    
    status <- rep('null', length(Chr_Peaks$chrom))
    signif <- rep('null', length(Chr_Peaks$chrom))
    status[which(IP_C$Bias>=0)] <- 'lagging_bias'; status[which(IP_C$Bias<0)] <- 'leading_bias' # Tha bias is positive: lagging, the bias is negative: leading
    signif[which(IP_C$p_value<=10e-6)] <- 'significant'; signif[which(IP_C$p_value>10e-6)] <- 'not_signif' # Establishing significance of the values by the p value 
    
    IP_C$description <- paste0(signif, "_", status)
    
    IP_T <- rbind(IP_T, IP_C)
    IP_T <- as.data.frame(IP_T)
  }
  return(IP_T)
}
##
IP_IdB_Early <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = TRUE)
Input_IdB_Early <- Bias_at_individual_peaks(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = E_Ori, Normalise2Input = TRUE)
SUP_IdB_Early <- Bias_at_individual_peaks(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = E_Ori, Normalise2Input = TRUE)

IP_IdB_Late <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = TRUE)
Input_IdB_Late <- Bias_at_individual_peaks(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = L_Ori, Normalise2Input = TRUE)
SUP_IdB_Late <- Bias_at_individual_peaks(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = L_Ori, Normalise2Input = TRUE)

## Plotting results
PlotEnrichments <- function(DataFile, PlotHeader){
  
  Watson <- round(as.numeric(smooth.spline(1:length(DataFile$watson.median), DataFile$watson.median)$y), 2) # Medians of both W and C
  Crick <- round(as.numeric(smooth.spline(1:length(DataFile$crick.median), DataFile$crick.median)$y), 2)
  
  Wat25 <- round(as.numeric(smooth.spline(1:length(DataFile$watson.q25), DataFile$watson.q25)$y), 2) # Watson quartile 25 and 75
  Wat75 <- round(as.numeric(smooth.spline(1:length(DataFile$watson.q75), DataFile$watson.q75)$y), 2)
  
  Cri25 <- round(as.numeric(smooth.spline(1:length(DataFile$crick.q25), DataFile$crick.q25)$y), 2) # Crick quartile 25 and 75
  Cri75 <- round(as.numeric(smooth.spline(1:length(DataFile$crick.q75), DataFile$crick.q75)$y), 2)
  
  Y <- max(round(abs(range(c(Wat75, Cri75*(-1))))+0.5)) # We set the Y that will establish the Y lims in the plot
  
  plot(Watson,
       ylim = c(-Y, Y),
       main = PlotHeader,
       ylab = "Average Enrichment", cex.main=1, xlab = "Distance from BrDU peakSummit (Kbp)", # Axis tittles 
       xaxt = "n", col = 'brown3', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs = 'i', yaxs = 'i')
  
  lines(Wat25, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2)) # Watson as red
  lines(Wat75, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
  
  polygon(x = c(1:length(Watson), rev(1:length(Watson))), # We draw the figure with the limits establish around the Watson line representation
          y = c(Wat25, rev(Wat75)), 
          col = adjustcolor("red", alpha.f = 0.2), border = NA)
  
  lines((Crick)*(-1), lwd=2, col = 'cornflowerblue', type = 'l') # Crick as blue
  lines((Cri25)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
  lines((Cri75)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
  
  polygon(x = c(1:length(Crick*(-1)), rev(1:length(Crick*(-1)))), # Figure between Crick limits
          y = c(Cri25*(-1), rev(Cri75*(-1))), 
          col = adjustcolor("cornflowerblue", alpha.f = 0.2), border = NA)
  
  rect(0, -Y, (length(Watson))/2 - (LeftLim/stepSize), Y, col = adjustcolor("grey", alpha.f = 0.65), border = NA) # We draw grey zones (rectangles) at both sides. This rectangles represent the length of W minus the left limit divided by 10, the result is a little bit more than a quarter axis
  rect((length(Watson))/2 + (RightLim/stepSize), -Y, length(Watson), Y, 
       col = adjustcolor("grey", alpha.f = 0.65), border = NA)
  
  text((AveragingWindow)*2-500, Y-1, labels = "Watson", cex = 0.9, col = 'brown3')
  text((AveragingWindow)*2-500, -Y+1, labels = "Crick", cex = 0.9, col = 'cornflowerblue')
  
  abline(h=0,lwd=0.4); abline(v=(length(Watson))/2,lwd=0.4) # We set the axis placing the origin in the center of the averagin window + and - (-3,3)
  axisLabels <- seq(-AveragingWindow,
                    +AveragingWindow,
                    length.out = 9)
  axisLabels[c(2,4,6,8)] <- NA
  At <- (AveragingWindow/stepSize)*seq(0,2,0.25); At[1] <- 1
  axis(1, at=At, labels = signif(axisLabels/1000, 2))
  
}

# Second function to plot Watson over Crick

PlotAverages <- function(DataFile, PlotHeader){
  # As before, we plot the median, the quartile 25 nd the q 75. In the other case we were plotting both W and C, in this case we plot the log2 got before, so we do not need to add two datasets
  Med <- round(as.numeric(smooth.spline(1:length(DataFile$median), DataFile$median)$y), 2)
  q25 <- round(as.numeric(smooth.spline(1:length(DataFile$q25), DataFile$q25)$y), 2)
  q75 <- round(as.numeric(smooth.spline(1:length(DataFile$q75), DataFile$q75)$y), 2)
  
  Y <- max(round(abs(range(q75))+0.5))
  
  plot(Med,
       ylim = c(-Y, +Y),
       main = PlotHeader,
       ylab = "log2 watson/crick", cex.main=1, xlab = "Distance from BrDU peakSummit (Kbp)", # Axis tittles 
       xaxt = "n", col = 'blue', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs='i', yaxs='i')
  
  polygon(x = c(1:length(Med), rev(1:length(Med))), # The polygon here is not comprehended between the axis and the W or C, is around the line 
          y = c(q25, rev(q75)), # Limits q.25 and 75
          col = adjustcolor("blue", alpha.f = 0.2), border = NA)
  
  rect(0, -Y, (length(Med))/2 - (LeftLim/stepSize), Y, col = adjustcolor("grey", alpha.f = 0.65), border = NA) # We add the grey rectangles as before
  rect((length(Med))/2 + (RightLim/stepSize), -Y, length(Med), Y, 
       col = adjustcolor("grey", alpha.f = 0.65), border = NA)
  
  abline(h=0, lwd=0.4); abline(v=(length(Med))/2,lwd=0.4) # Same as before
  axisLabels <- seq(-AveragingWindow,
                    +AveragingWindow,
                    length.out = 9)
  axisLabels[c(2,4,6,8)] <- NA
  At <- (AveragingWindow/stepSize)*seq(0,2,0.25); At[1] <- 1
  axis(1, at=At, labels = signif(axisLabels/1000, 2))
  
}

# Plotting the Bias 
PlotIdBias <- function(DataFile, PlotHeader){
  
  DataFile$Bias[!is.finite(DataFile$Bias)] <- 0
  DataFile$Bias[is.na(DataFile$Bias)] <- 0
  
  boxplot(DataFile$Bias, ylim = c(-2,2), las = 2, ylab = 'log2 lagging/leading', # We change the plot type to boxplot
          cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
  # Not lagging nor leading
  if(length(DataFile$Bias[which(DataFile$p_value > 10e-6)])>0){ # For those containing a p value over the limit
    spreadPoints(values=DataFile$Bias[which(DataFile$p_value > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1) # We draw them in blue 
  }
  # Lagging
  if(length(DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)])>0){ 3 # For those with a bias over or equal 0 and a p value below or equal to the limit
    spreadPoints(values=DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)], position=1.0, pointCex=0.65, col="red", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1) # We draw them in red
  }
  # Leading
  if(length(DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)])>0){ # For those with a bias below the limitand a p value below or equal to the limit 
    spreadPoints(values=DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)], position=1.0, pointCex=0.65, col="green", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1) # We draw them in green
  }
  
  # We defined what were lagging, leading and indt according to their p value 
  Lagg <- length(DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)]) # Lagging are the ones with a bias over 0
  Lead <- length(DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)]) # Leading are the ones with a bias below 0
  Indt <- length(DataFile$Bias[which(DataFile$p_value > 10e-6)])
  # Add the description
  legend("bottomright", legend = c(paste0("lagg", " (", Lagg, ")"), 
                                   paste0("lead", " (", Lead, ")"), 
                                   paste0("inde", " (", Indt, ")")), 
         col = c("red", "green", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
  
  
  InName <- PlotHeader
  
  #calculate p value
  #Decision Tree
  biasP <- binom.test(round(c(Lagg+Lead, Indt)))$p.value # By using the binomial test we conclude if the bias is significative or not
  
  if(biasP <= 10e-3 & (Lagg+Lead) > Indt){
    
    biasQ <- binom.test(round(c(Lagg, Lead)))$p.value
    
    if(biasQ <= 10e-6){
      if(Lagg > Lead){
        conc <- paste0("Strong bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3)) # Lagging bias higher than leading bias
      }
      if(Lead > Lagg){
        conc <- paste0("Strong bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3)) # Leading bias higher than lagging bias
      }
    }
    # Weaker bias as the p value is not below 10e-6
    if(biasQ > 10e-6 & biasQ <= 10e-4){
      if(Lagg > Lead){
        conc <- paste0("Weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
      }
      if(Lead > Lagg){
        conc <- paste0("Weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
      }
    }
    # Weaker as it is between 10e-4 and 10e-2
    if(biasQ > 10e-4  & biasQ <= 10e-2){
      if(Lagg > Lead){
        conc <- paste0("Very weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
      }
      if(Lead > Lagg){
        conc <- paste0("Very weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
      }
    }
    # No bias
    if(biasQ > 10e-2){
      conc <- paste0("No significant strandedness in DNA synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
    }
  } 
  else 
  {
    conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
  }
  
  
  
  mtext(conc, side = 1, line = 2, cex = 0.65)
}


pdf(paste0("~/Desktop/", Pro_1, "/", Pro_1, "_Results.pdf"), width = 10, height = 12)

par(oma=c(0,0,0,0))

PlotMat <- {matrix(c(     0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                          
                          0,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,0,
                          0,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,0,
                          0,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,0,
                          0,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,0,
                          0,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,0,
                          
                          0,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,0,
                          0,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,0,
                          0,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,0,
                          0,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,0,
                          0,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,0,
                          
                          0,0,0,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,0,0,
                          0,0,0,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,0,0,
                          0,0,0,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,0,0,
                          0,0,0,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,0,0,
                          0,0,0,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,0,0), 
                   
                   20,20,byrow=TRUE)}

layout(PlotMat, c(1,1), c(1,1), TRUE)

plot(image, axes = FALSE)
PlotEnrichments(Input_AvE_Early, "Input_Early"); 
PlotEnrichments(IP_AvE_Early, "IP_Early"); PlotEnrichments(SUP_AvE_Early, "SUP_Early")
PlotAverages(Input_WoC_Early, "Input_Early"); 
PlotAverages(IP_WoC_Early, "IP_Early"); PlotAverages(SUP_WoC_Early, "SUP_Early")
PlotIdBias(Input_IdB_Early, "Input_Early"); 
PlotIdBias(IP_IdB_Early, "IP_Early"); PlotIdBias(SUP_IdB_Early, "SUP_Early")

plot(image, axes = FALSE)
PlotEnrichments(Input_AvE_Late, "Input_Late"); 
PlotEnrichments(IP_AvE_Late, "IP_Late"); PlotEnrichments(SUP_AvE_Late, "SUP_Late")
PlotAverages(Input_WoC_Late, "Input_Late"); 
PlotAverages(IP_WoC_Late, "IP_Late"); PlotAverages(SUP_WoC_Late, "SUP_Late")
PlotIdBias(Input_IdB_Late, "Input_Late"); 
PlotIdBias(IP_IdB_Late, "IP_Late"); PlotIdBias(SUP_IdB_Late, "SUP_Late")

## Additional info included in the pdf
plot(NULL, xlim=c(0,2), ylim=c(0,2), ylab=" ", xlab=" ", yaxt="n", xaxt="n")
txt <- paste0("99% BrDU Signal are within ", LeftLim, " (left) and ", RightLim, " (right)", " bps from center", "\n",
              "Leading synthesis = ", LeadingAverage, " bps;", " Lagging synthesis = ", LaggingAverage, " bps", "\n",
              ", ", "binSize = ", binSize, " bps", ", ", "slide = ", stepSize, " bps")
text(1, 1, labels = txt, cex = 1)

dev.off()




message("Analysis complete.")

message(paste0("Check results at Desktop folder - ", Pro_1))

