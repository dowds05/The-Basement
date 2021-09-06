cd /Users/lab/Desktop/Bob_Aim3/Linux/DATA # cd to directory

echo $PATH # Set the PATH
export PATH=$PATH:/Users/lab/opt/minimap2:/Users/lab/Desktop/Bob_Aim3/Linux/samtools:/Users/lab/Desktop/Bob_Aim3/Linux/bcftools:/Users/lab/Desktop/Bob_Aim3/Linux/bbmap
##################################################################################################################
##                                grep to search Amarel                     
##################################################################################################################
# Open another terminal
ssh rad267@amarel.rutgers.edu # ssh to Amarel
cd /scratch/rad267/FASTA # cd to scratch
# line count for files
wc -l search.txt
wc -l search.fasta
wc -m ref.fa
# scp search.txt from LOCAL directory to the amarel server 
cd /Users/lab/Desktop/Bob_Aim3/Linux/DATA
scp query.fa rad267@amarel.rutgers.edu:/scratch/rad267/FASTA/ # Send to Amarel

# From /scratch directory search for taxa of interest on AMAREL
srun grep -Ff query.fa -A1 tot_fasta > search.fasta # Run grep

# scp to LOCAL workstation directory & run. Note: "." = here 
scp rad267@amarel.rutgers.edu:/scratch/rad267/FASTA/search.fasta . # Download to local

##################################################################################################################
##                                BBDuk                    
##################################################################################################################
# https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
# https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk.sh
# “Duk” stands for Decontamination Using Kmers. 
# used this to install
#/Users/lab/Desktop/Bob_Aim3/Linux/bbmap/stats.sh in=/Users/lab/Desktop/Bob_Aim3/Linux/bbmap/resources/phix174_ill.ref.fa.gz

# Adapter trimming
# bbduk.sh in=reads.fq out=clean.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
/Users/lab/Desktop/Bob_Aim3/Linux/bbmap/bbduk.sh -Xmx12g in=/Users/lab/Desktop/Bob_Aim3/Linux/DATA/search.fasta out=/Users/lab/Desktop/Bob_Aim3/Linux/DATA/clean.fasta outm=/Users/lab/Desktop/Bob_Aim3/Linux/DATA/unclean.fasta ref=/Users/lab/Desktop/Bob_Aim3/Linux/DATA/barcodes.fasta ktrim=N k=24 hdist=1 maskmiddle=t fixjunk

##################################################################################################################
##                                minimap2                    
##################################################################################################################
# https://github.com/lh3/minimap2/blob/master/README.md#uguide
# minimap2 is used to ORIENT reads to the index sequence downloaded from NCBI or on Geneious
/Users/lab/Desktop/Bob_Aim3/Linux/minimap2/minimap2 -d /Users/lab/Desktop/Bob_Aim3/Linux/DATA/search.mmi /Users/lab/Desktop/Bob_Aim3/Linux/DATA/ref.fa # indexing

# Verify line numbers by dividing fasta line count by two. use bc, which can do math on command line. 
# https://www.networkworld.com/article/3268964/how-to-do-math-on-the-linux-command-line.html
bc # An arbitrary precision calculator language

# create alignment to ORIENT reads & IMPORT into Geneious 
/Users/lab/Desktop/Bob_Aim3/Linux/minimap2/minimap2 -a /Users/lab/Desktop/Bob_Aim3/Linux/DATA/search.mmi /Users/lab/Desktop/Bob_Aim3/Linux/DATA/clean.fasta > /Users/lab/Desktop/Bob_Aim3/Linux/DATA/orient.sam # orientation

# FROM Geneious, annotated/sorted/dedup (i.e. ann.fa) to create ALIGNMENT 
wc -l ann.fa
/Users/lab/Desktop/Bob_Aim3/Linux/minimap2/minimap2 -a /Users/lab/Desktop/Bob_Aim3/Linux/DATA/search.mmi /Users/lab/Desktop/Bob_Aim3/Linux/DATA/ann.fa > /Users/lab/Desktop/Bob_Aim3/Linux/DATA/aln.sam # alignment

# Run unambig.sh to remove "N, ?, K, R, M, W, Y" 
bash /Users/lab/Desktop/Bob_Aim3/Linux/Scripts/unambig.sh

# Rename Consensus.fasta to TREATMENT & IMPORT to Geneious
cat Consensus.fa > VHFX-M-12_B.fa
