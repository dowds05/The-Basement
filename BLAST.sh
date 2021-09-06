##################################################################################################################
##                                Basic Local Alignment Search Tool (BLAST)                     
##################################################################################################################
cd /Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA # cd to directory
cd /Users/lab/Desktop/Bob_Aim3/LRC_BLAST/ # cd to directory

echo $PATH # Set the PATH
export PATH=$PATH:/Users/lab/BLAST_2.10.0/usr/local/ncbi/blast/bin

##################################################################################################################
##                                grep to search Amarel                     
##################################################################################################################
ssh rad267@amarel.rutgers.edu # ssh to Amarel
cd /scratch/rad267/FASTA # cd to scratch
# line count for files
wc -l taxa.txt
wc -l search.fasta
wc -m ref.fa
# scp search.txt from LOCAL directory to the amarel server 
cd /Users/lab/Desktop/Bob_Aim3/Linux/DATA
scp taxa.txt rad267@amarel.rutgers.edu:/scratch/rad267/FASTA/ # Send to Amarel

# From /scratch directory search for taxa of interest on AMAREL
srun grep -Ff taxa.txt -A1 tot_fasta > taxa.fa # Run grep
# (OR) run as a batch job
sbatch SEARCH.sh 

squeue -u rad267 # Check the activity of the job

# scp to LOCAL workstation directory & run. Note: "." = here 
scp rad267@amarel.rutgers.edu:/scratch/rad267/FASTA/taxa.fa . # Download to local

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
/Users/lab/Desktop/Bob_Aim3/Linux/bbmap/bbduk.sh -Xmx12g in=/Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA/taxa.fa out=/Users/lab/Desktop/Bob_Aim3/Linux/DATA/clean.fasta outm=/Users/lab/Desktop/Bob_Aim3/Linux/DATA/unclean.fasta ref=/Users/lab/Desktop/Bob_Aim3/Linux/DATA/barcodes.fasta ktrim=N k=24 hdist=1 maskmiddle=t fixjunk


##################################################################################################################
##                                MegaBLAST to ID very similar sequences                     
##################################################################################################################

# Make database to query against
makeblastdb -in /Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA/taxa.fa -dbtype nucl # Make database to query against
makeblastdb -in /Users/lab/rOpDB/rOpDB.fasta -dbtype nucl # rOpDB database to query against
makeblastdb -in /Users/lab/Desktop/Bob_Aim3/LRC_BLAST/Muri_species.fasta -dbtype nucl # Muribac species database


# Run megablast
/Users/lab/BLAST_2.10.0/usr/local/ncbi/blast/bin/blastn -db /Users/lab/rOpDB/rOpDB.fasta -query /Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA/taxa.fa -out /Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA/taxa.out.csv -task megablast -word_size 36 -reward 1 -penalty -4 -gapopen 0 -gapextend 2 -evalue 1e-5 -outfmt "10 qacc sacc salltitles pident length evalue" -max_hsps 10 -max_target_seqs 10 # Run megablast


# Run blastn
/Users/lab/BLAST_2.10.0/usr/local/ncbi/blast/bin/blastn -db /Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA/taxa.fa -query /Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA/query.fa -out /Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA/query.out.csv -task blastn -max_target_seqs 10 -evalue 0.25 -perc_identity 95 -outfmt "10 qacc sacc pident length evalue" # Run blastn

# blastn -query  genomic.sequences.2.fasta -db genomic.sequences.1.fasta -task blastn -outfmt 7 -max_target_seqs 10 -evalue 0.5 -perc_identity 95 > blast.out
