#!/bin/bash

#SBATCH --clusters=amarel            # Select which system(s) to use
#SBATCH --partition=main             # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=RENAME            # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=10           # Cores per task (>1 if multithread tasks)
#SBATCH --mem=16G                     # Real memory (RAM) required per node
#SBATCH --time=00:60:00              # Total run time limit (DD-HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT file for SLURM output
#SBATCH --mail-user=rad267@rutgers.edu # Send me email on output
#SBATCH --mail-type=ALL         # Or NONE,BEGIN,END,FAIL,REQUEUE,ALL

## Run the job
srun
# Week 3

# Use this to add sample/treatment/information to BLAST outputs
# CD-S, F
awk -F "," 'BEGIN { OFS = "," } {$8="1"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_3"; $13="101";  print}' 1.out.csv > 1.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="2"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_3"; $13="102";  print}' 2.out.csv > 2.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="3"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_3"; $13="103";  print}' 3.out.csv > 3.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="4"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_3"; $13="104";  print}' 4.out.csv > 4.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="5"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_3"; $13="105";  print}' 5.out.csv > 5.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="6"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_3"; $13="106";  print}' 6.out.csv > 6.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="7"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_3"; $13="107";  print}' 7.out.csv > 7.new.csv
# CD-S, M
awk -F "," 'BEGIN { OFS = "," } {$8="8"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_3"; $13="108";  print}' 8.out.csv > 8.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="9"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_3"; $13="109";  print}' 9.out.csv > 9.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="10"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_3"; $13="110";  print}' 10.out.csv > 10.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="11"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_3"; $13="111";  print}' 11.out.csv > 11.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="12"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_3"; $13="112";  print}' 12.out.csv > 12.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="13"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_3"; $13="113";  print}' 13.out.csv > 13.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="14"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_3"; $13="114";  print}' 14.out.csv > 14.new.csv
# VHF-S, F
awk -F "," 'BEGIN { OFS = "," } {$8="15"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="115";  print}' 15.out.csv > 15.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="16"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="116";  print}' 16.out.csv > 16.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="17"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="117";  print}' 17.out.csv > 17.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="18"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="118";  print}' 18.out.csv > 18.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="19"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="119";  print}' 19.out.csv > 19.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="20"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="120";  print}' 20.out.csv > 20.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="21"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="121";  print}' 21.out.csv > 21.new.csv
# VHF-S, M
awk -F "," 'BEGIN { OFS = "," } {$8="22"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="122";  print}' 22.out.csv > 22.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="23"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="123";  print}' 23.out.csv > 23.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="24"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="124";  print}' 24.out.csv > 24.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="25"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="125";  print}' 25.out.csv > 25.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="26"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="126";  print}' 26.out.csv > 26.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="27"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="127";  print}' 27.out.csv > 27.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="28"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="128";  print}' 28.out.csv > 28.new.csv
# CD-X, F
awk -F "," 'BEGIN { OFS = "," } {$8="29"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_3"; $13="129";  print}' 29.out.csv > 29.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="30"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_3"; $13="130";  print}' 30.out.csv > 30.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="31"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_3"; $13="131";  print}' 31.out.csv > 31.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="32"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_3"; $13="132";  print}' 32.out.csv > 32.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="33"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_3"; $13="133";  print}' 33.out.csv > 33.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="34"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_3"; $13="134";  print}' 34.out.csv > 34.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="35"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_3"; $13="135";  print}' 35.out.csv > 35.new.csv
# CD-X, M
awk -F "," 'BEGIN { OFS = "," } {$8="36"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_3"; $13="136";  print}' 36.out.csv > 36.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="37"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_3"; $13="137";  print}' 37.out.csv > 37.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="38"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_3"; $13="138";  print}' 38.out.csv > 38.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="39"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_3"; $13="139";  print}' 39.out.csv > 39.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="40"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_3"; $13="140";  print}' 40.out.csv > 40.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="41"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_3"; $13="141";  print}' 41.out.csv > 41.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="42"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_3"; $13="142";  print}' 42.out.csv > 42.new.csv
# VHF-X, F
awk -F "," 'BEGIN { OFS = "," } {$8="43"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="143";  print}' 43.out.csv > 43.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="44"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="144";  print}' 44.out.csv > 44.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="45"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="145";  print}' 45.out.csv > 45.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="46"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="146";  print}' 46.out.csv > 46.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="47"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="147";  print}' 47.out.csv > 47.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="48"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="148";  print}' 48.out.csv > 48.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="49"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_3"; $13="149";  print}' 49.out.csv > 49.new.csv
# VHF-X, M
awk -F "," 'BEGIN { OFS = "," } {$8="50"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="150";  print}' 50.out.csv > 50.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="51"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="151";  print}' 51.out.csv > 51.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="52"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="152";  print}' 52.out.csv > 52.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="53"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="153";  print}' 53.out.csv > 53.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="54"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="154";  print}' 54.out.csv > 54.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="55"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="155";  print}' 55.out.csv > 55.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="56"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_3"; $13="156";  print}' 56.out.csv > 56.new.csv

# Week 4 
# CD-S, F
awk -F "," 'BEGIN { OFS = "," } {$8="57"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_4"; $13="101";  print}' 57.out.csv > 57.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="58"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_4"; $13="102";  print}' 58.out.csv > 58.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="59"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_4"; $13="103";  print}' 59.out.csv > 59.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="60"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_4"; $13="104";  print}' 60.out.csv > 60.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="61"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_4"; $13="105";  print}' 61.out.csv > 61.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="62"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_4"; $13="106";  print}' 62.out.csv > 62.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="63"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_4"; $13="107";  print}' 63.out.csv > 63.new.csv
# CD-S, M
awk -F "," 'BEGIN { OFS = "," } {$8="64"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_4"; $13="108";  print}' 64.out.csv > 64.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="65"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_4"; $13="109";  print}' 65.out.csv > 65.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="66"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_4"; $13="110";  print}' 66.out.csv > 66.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="67"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_4"; $13="111";  print}' 67.out.csv > 67.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="68"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_4"; $13="112";  print}' 68.out.csv > 68.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="69"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_4"; $13="113";  print}' 69.out.csv > 69.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="70"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_4"; $13="114";  print}' 70.out.csv > 70.new.csv
# VHF-S, F
awk -F "," 'BEGIN { OFS = "," } {$8="71"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="115";  print}' 71.out.csv > 71.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="72"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="116";  print}' 72.out.csv > 72.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="73"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="117";  print}' 73.out.csv > 73.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="74"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="118";  print}' 74.out.csv > 74.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="75"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="119";  print}' 75.out.csv > 75.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="76"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="120";  print}' 76.out.csv > 76.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="77"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="121";  print}' 77.out.csv > 77.new.csv
# VHF-S, M
awk -F "," 'BEGIN { OFS = "," } {$8="78"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="122";  print}' 78.out.csv > 78.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="79"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="123";  print}' 79.out.csv > 79.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="80"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="124";  print}' 80.out.csv > 80.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="81"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="125";  print}' 81.out.csv > 81.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="82"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="126";  print}' 82.out.csv > 82.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="83"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="127";  print}' 83.out.csv > 83.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="84"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="128";  print}' 84.out.csv > 84.new.csv
# CD-X, F
awk -F "," 'BEGIN { OFS = "," } {$8="85"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_4"; $13="129";  print}' 85.out.csv > 85.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="86"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_4"; $13="130";  print}' 86.out.csv > 86.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="87"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_4"; $13="131";  print}' 87.out.csv > 87.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="88"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_4"; $13="132";  print}' 88.out.csv > 88.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="89"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_4"; $13="133";  print}' 89.out.csv > 89.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="90"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_4"; $13="134";  print}' 90.out.csv > 90.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="91"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_4"; $13="135";  print}' 91.out.csv > 91.new.csv
# CD-X, M
awk -F "," 'BEGIN { OFS = "," } {$8="92"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_4"; $13="136";  print}' 92.out.csv > 92.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="93"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_4"; $13="137";  print}' 93.out.csv > 93.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="94"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_4"; $13="138";  print}' 94.out.csv > 94.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="95"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_4"; $13="139";  print}' 95.out.csv > 95.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="96"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_4"; $13="140";  print}' 96.out.csv > 96.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="97"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_4"; $13="141";  print}' 97.out.csv > 97.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="98"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_4"; $13="142";  print}' 98.out.csv > 98.new.csv
# VHF-X, F
awk -F "," 'BEGIN { OFS = "," } {$8="99"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="143";  print}' 99.out.csv > 99.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="100"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="144";  print}' 100.out.csv > 100.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="101"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="145";  print}' 101.out.csv > 101.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="102"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="146";  print}' 102.out.csv > 102.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="103"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="147";  print}' 103.out.csv > 103.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="104"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="148";  print}' 104.out.csv > 104.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="105"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_4"; $13="149";  print}' 105.out.csv > 105.new.csv
# VHF-X, M
awk -F "," 'BEGIN { OFS = "," } {$8="106"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="150";  print}' 106.out.csv > 106.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="107"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="151";  print}' 107.out.csv > 107.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="108"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="152";  print}' 108.out.csv > 108.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="109"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="153";  print}' 109.out.csv > 109.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="110"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="154";  print}' 110.out.csv > 110.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="111"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="155";  print}' 111.out.csv > 111.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="112"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_4"; $13="156";  print}' 112.out.csv > 112.new.csv

# WEEK 8
# CD-S, F
awk -F "," 'BEGIN { OFS = "," } {$8="113"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_8"; $13="101";  print}' 113.out.csv > 113.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="114"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_8"; $13="102";  print}' 114.out.csv > 114.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="115"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_8"; $13="103";  print}' 115.out.csv > 115.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="116"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_8"; $13="104";  print}' 116.out.csv > 116.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="117"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_8"; $13="105";  print}' 117.out.csv > 117.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="118"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_8"; $13="106";  print}' 118.out.csv > 118.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="119"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_8"; $13="107";  print}' 119.out.csv > 119.new.csv
# CD-S, M
awk -F "," 'BEGIN { OFS = "," } {$8="120"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_8"; $13="108";  print}' 120.out.csv > 120.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="121"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_8"; $13="109";  print}' 121.out.csv > 121.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="122"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_8"; $13="110";  print}' 122.out.csv > 122.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="123"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_8"; $13="111";  print}' 123.out.csv > 123.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="124"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_8"; $13="112";  print}' 124.out.csv > 124.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="125"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_8"; $13="113";  print}' 125.out.csv > 125.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="126"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_8"; $13="114";  print}' 126.out.csv > 126.new.csv
# VHF-S, F
awk -F "," 'BEGIN { OFS = "," } {$8="127"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="115";  print}' 127.out.csv > 127.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="128"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="116";  print}' 128.out.csv > 128.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="129"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="117";  print}' 129.out.csv > 129.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="130"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="118";  print}' 130.out.csv > 130.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="131"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="119";  print}' 131.out.csv > 131.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="132"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="120";  print}' 132.out.csv > 132.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="133"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="121";  print}' 133.out.csv > 133.new.csv
# VHF-S, M
awk -F "," 'BEGIN { OFS = "," } {$8="134"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="122";  print}' 134.out.csv > 134.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="135"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="123";  print}' 135.out.csv > 135.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="136"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="124";  print}' 136.out.csv > 136.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="137"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="125";  print}' 137.out.csv > 137.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="138"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="126";  print}' 138.out.csv > 138.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="139"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="127";  print}' 139.out.csv > 139.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="140"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="128";  print}' 140.out.csv > 140.new.csv
# CD-X, F
awk -F "," 'BEGIN { OFS = "," } {$8="141"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_8"; $13="129";  print}' 141.out.csv > 141.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="142"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_8"; $13="130";  print}' 142.out.csv > 142.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="143"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_8"; $13="131";  print}' 143.out.csv > 143.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="144"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_8"; $13="132";  print}' 144.out.csv > 144.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="145"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_8"; $13="133";  print}' 145.out.csv > 145.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="146"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_8"; $13="134";  print}' 146.out.csv > 146.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="147"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_8"; $13="135";  print}' 147.out.csv > 147.new.csv
# CD-X, M
awk -F "," 'BEGIN { OFS = "," } {$8="148"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_8"; $13="136";  print}' 148.out.csv > 148.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="149"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_8"; $13="137";  print}' 149.out.csv > 149.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="150"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_8"; $13="138";  print}' 150.out.csv > 150.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="151"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_8"; $13="139";  print}' 151.out.csv > 151.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="152"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_8"; $13="140";  print}' 152.out.csv > 152.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="153"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_8"; $13="141";  print}' 153.out.csv > 153.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="154"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_8"; $13="142";  print}' 154.out.csv > 154.new.csv
# VHF-X, F
awk -F "," 'BEGIN { OFS = "," } {$8="155"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="143";  print}' 155.out.csv > 155.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="156"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="144";  print}' 156.out.csv > 156.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="157"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="145";  print}' 157.out.csv > 157.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="158"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="146";  print}' 158.out.csv > 158.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="159"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="147";  print}' 159.out.csv > 159.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="160"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="148";  print}' 160.out.csv > 160.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="161"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_8"; $13="149";  print}' 161.out.csv > 161.new.csv
# VHF-X, M
awk -F "," 'BEGIN { OFS = "," } {$8="162"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="150";  print}' 162.out.csv > 162.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="163"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="151";  print}' 163.out.csv > 163.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="164"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="152";  print}' 164.out.csv > 164.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="165"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="153";  print}' 165.out.csv > 165.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="166"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="154";  print}' 166.out.csv > 166.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="167"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="155";  print}' 167.out.csv > 167.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="168"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_8"; $13="156";  print}' 168.out.csv > 168.new.csv

# WEEK 12
# CD-S, F
awk -F "," 'BEGIN { OFS = "," } {$8="169"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_12"; $13="101";  print}' 169.out.csv > 169.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="170"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_12"; $13="102";  print}' 170.out.csv > 170.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="171"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_12"; $13="103";  print}' 171.out.csv > 171.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="172"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_12"; $13="104";  print}' 172.out.csv > 172.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="173"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_12"; $13="105";  print}' 173.out.csv > 173.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="174"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_12"; $13="106";  print}' 174.out.csv > 174.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="175"; $9="Sedentary"; $10="Control"; $11="Female"; $12="Week_12"; $13="107";  print}' 175.out.csv > 175.new.csv
# CD-S, M
awk -F "," 'BEGIN { OFS = "," } {$8="176"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_12"; $13="108";  print}' 176.out.csv > 176.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="177"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_12"; $13="109";  print}' 177.out.csv > 177.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="178"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_12"; $13="110";  print}' 178.out.csv > 178.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="179"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_12"; $13="111";  print}' 179.out.csv > 179.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="180"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_12"; $13="112";  print}' 180.out.csv > 180.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="181"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_12"; $13="113";  print}' 181.out.csv > 181.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="182"; $9="Sedentary"; $10="Control"; $11="Male"; $12="Week_12"; $13="114";  print}' 182.out.csv > 182.new.csv
# VHF-S, F
awk -F "," 'BEGIN { OFS = "," } {$8="183"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="115";  print}' 183.out.csv > 183.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="184"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="116";  print}' 184.out.csv > 184.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="185"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="117";  print}' 185.out.csv > 185.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="186"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="118";  print}' 186.out.csv > 186.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="187"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="119";  print}' 187.out.csv > 187.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="188"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="120";  print}' 188.out.csv > 188.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="189"; $9="Sedentary"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="121";  print}' 189.out.csv > 189.new.csv
# VHF-S, M
awk -F "," 'BEGIN { OFS = "," } {$8="190"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="122";  print}' 190.out.csv > 190.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="191"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="123";  print}' 191.out.csv > 191.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="192"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="124";  print}' 192.out.csv > 192.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="193"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="125";  print}' 193.out.csv > 193.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="194"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="126";  print}' 194.out.csv > 194.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="195"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="127";  print}' 195.out.csv > 195.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="196"; $9="Sedentary"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="128";  print}' 196.out.csv > 196.new.csv
# CD-X, F
awk -F "," 'BEGIN { OFS = "," } {$8="197"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_12"; $13="129";  print}' 197.out.csv > 197.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="198"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_12"; $13="130";  print}' 198.out.csv > 198.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="199"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_12"; $13="131";  print}' 199.out.csv > 199.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="200"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_12"; $13="132";  print}' 200.out.csv > 200.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="201"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_12"; $13="133";  print}' 201.out.csv > 201.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="202"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_12"; $13="134";  print}' 202.out.csv > 202.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="203"; $9="Exercise"; $10="Control"; $11="Female"; $12="Week_12"; $13="135";  print}' 203.out.csv > 203.new.csv
# CD-X, M
awk -F "," 'BEGIN { OFS = "," } {$8="204"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_12"; $13="136";  print}' 204.out.csv > 204.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="205"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_12"; $13="137";  print}' 205.out.csv > 205.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="206"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_12"; $13="138";  print}' 206.out.csv > 206.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="207"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_12"; $13="139";  print}' 207.out.csv > 207.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="208"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_12"; $13="140";  print}' 208.out.csv > 208.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="209"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_12"; $13="141";  print}' 209.out.csv > 209.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="210"; $9="Exercise"; $10="Control"; $11="Male"; $12="Week_12"; $13="142";  print}' 210.out.csv > 210.new.csv
# VHF-X, F
awk -F "," 'BEGIN { OFS = "," } {$8="211"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="143";  print}' 211.out.csv > 211.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="212"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="144";  print}' 212.out.csv > 212.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="213"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="145";  print}' 213.out.csv > 213.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="214"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="146";  print}' 214.out.csv > 214.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="215"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="147";  print}' 215.out.csv > 215.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="216"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="148";  print}' 216.out.csv > 216.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="217"; $9="Exercise"; $10="VeryHighFat"; $11="Female"; $12="Week_12"; $13="149";  print}' 217.out.csv > 217.new.csv
# VHF-X, M
awk -F "," 'BEGIN { OFS = "," } {$8="218"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="150";  print}' 218.out.csv > 218.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="219"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="151";  print}' 219.out.csv > 219.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="220"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="152";  print}' 220.out.csv > 220.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="221"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="153";  print}' 221.out.csv > 221.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="222"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="154";  print}' 222.out.csv > 222.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="223"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="155";  print}' 223.out.csv > 223.new.csv
awk -F "," 'BEGIN { OFS = "," } {$8="224"; $9="Exercise"; $10="VeryHighFat"; $11="Male"; $12="Week_12"; $13="156";  print}' 224.out.csv > 224.new.csv

# BASELINE SAMPLES
# Female
awk -F "," 'BEGIN { OFS = "," } {$8="225"; $9="Baseline"; $10="Baseline"; $11="Female"; $12="D_0"; $13="Female_baseline";  print}' 225.out.csv > 225.new.csv
# Male
awk -F "," 'BEGIN { OFS = "," } {$8="226"; $9="Baseline"; $10="Baseline"; $11="Male"; $12="D_0"; $13="Male_baseline";  print}' 226.out.csv > 226.new.csv

# Combine all the files using cat
cat *.new.csv > tot.csv
# Comfirm line numbers
wc -l tot.csv





