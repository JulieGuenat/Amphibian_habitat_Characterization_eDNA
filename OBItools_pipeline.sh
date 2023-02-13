# Processing amphibians sequences using OBItools

### 1. Cut sequences files
#Because we have a large amount of sequences (more 180 millions read), we need to split the mother file to run illuminapairedend command
#in parallel to reduce the computation time.
# To do so we first created a text file containing the line number at which we want to split the motherfile.
# We printed the number of lines contained in the mother file (730'689'392) and divided it by 100 (number of children file we wanted).
# Then, we created the text file using:
nano cut_delim.txt

# We changed the seperator, so the file would be line read. (change tab for ;) using:

cat cut_delim.txt | sed 's/\t/;/g' > cut_del.txt  # 's' means separate and 'g' means general

#Then Guillaume Lavanchy created the following loop to split the mother file:

!#/bin/bash/

# takes 3 arguments:
# 1: a file with three ";"-separated columns as input
# 2: the file with single reads (R1)
# 3: the file with paired reads (R2)

for i in $(cat $1)
do
START=$(echo $i | cut -f 2 -d ';')
NUM=$(echo $i | cut -f 1 -d ';')
tail -n +${START} $2 | head -n 7306894 > Batra_R1_${NUM}.fq
tail -n +${START} $3 | head -n 7306894 > Batra_R2_${NUM}.fq
done

#we end up with 100 files R1 and 100 files R2, named Batra_R*_x.fq. I then created a folder named R1_R2_Batr100 that contained all those files.

### 2. Runing the illuminapairedend command on bsub
# We were not able run the illu_pairedend_bsub_scripts (created by Dan), because Obitools is desactivated before the command.
# Thus we ran the illuminapairedend command via interactive jobs with bsub. Ran the following loop:
for i in $(ls Batra_R1_* | sed 's/Batra_R1_//g' | sed 's/.fq//g'); do bsub <<< "illuminapairedend --score-min=40 -r Batra_R2_${i}.fq Batra_R1_${i}.fq > batra_${i}.fq"; done

#The  output is the alignment of the single and paired sequences. Again I created a folder containing all those files named "Batra_illu_pairedend"

### 3. Concatenate all alligned files in one files
#Once the illuminapairedend command finished to run, we had to combine all files in one uniq file to continue the obitools analysis.
#To do so I used:
cat batr_* > batr.all.illu.fq

#To be able to see the distribution of the scores for the batr sequences, I used the following function to have the count for each score
bsub -n 4 -M 8000000 -J "obistat" 'obistat -c score batr_all_illu.fg > batr_all_illu_scores.txt'

#Then I have to modify the separator of the file to be able to further open it in R. Since it is space separated, but with differential amount f space, I have to cut the file by columns and then past the columns with ";" as delimiter.
awk '{print $1}' batr_all_illu_scores.txt > col1
awk '{print $2}' batr_all_illu_scores.txt > col2
awk '{print $3}' batr_all_illu_scores.txt > col3

paste -d";" col1 col2 > x
paste -d";" x col3 > batr_all_illu_scores_scsep.txt


### 4. Remove unaligned sequences record

nohup obigrep -p 'mode!="joined"' batr.illu.fq > batr_illu_ali.fq &
# The command above corresponds to the one for all sequences.
#'mode!="joined"' means that it will remove sequences that are only joined but not aligned.
#however, Obitools command does not support the nohup function, thus I ran this command with bsub.

#The following commands were run on the splited sequences files:
#I first tried to remove unaligned sequences using the parameter 'score>50'.
#for i in $(ls batra_* | sed 's/batra_//g' | sed 's/.fq//g'); do bsub <<< "obigrep -p 'score>50' batra_${i}.fq > batra_${i}_ali.fq"; done
#However, the threshold is set in an arbitrary manner, and I prefered to remove only sequences that were not alligned (that recieved the joined attribute during the illumina paired end processus)

#I performed obigrep by removing joined sequences and not by setting a threshold
for i in $(ls batra_*.fq | sed 's/batra_//g' | sed 's/.fq//g'); do bsub <<< "obigrep -p 'mode!=\"joined\"' batra_${i}.fq > batra_${i}_ali_rmjoined.fq"; done

#after this step I also concatenate all sequences together to be able to make a second histogramm with the scores distribution.
bsub -M 8000000 -J "obistat" 'obistat -c score batraALL_ali_rmjoined > batraALL_ali_rmjoined_countscore'
awk '{print $1}' batraALL_ali_rmjoined_countscore.txt > col1
awk '{print $2}' batraALL_ali_rmjoined_countscore.txt > col2
awk '{print $3}' batraALL_ali_rmjoined_countscore.txt > col3

paste -d";" col1 col2 > x
paste -d";" x col3 > batr_all_illu_scores_scsep.txt


### 5. Assign marker/sample combination
#for sequences selected with an alignment score above 50:
#for i in $(ls batra_*_ali.fq | sed 's/batra_//g' | sed 's/_ali.fq//g'); do bsub <<< "ngsfilter -t generic_ngsfilter_Batr01.txt -u unidentified_batr01.fq batra_${i}_ali.fq > batra_ali_assigned_${i}.fq"; done

#for sequences selected with the parameter 'mode!="joined"'
for i in $(ls batra_* | sed 's/batra_//g' | sed 's/_ali_rmjoined.fq//g'); do bsub <<< "ngsfilter -t generic_ngsfilter_Batr01.txt -u unidentified_Batra.fq batra_${i}_ali_rmjoined.fq > batra_ali_rmjoined_assigned_${i}.fq"; done

# We first needed to concatenated the sequences alligned, selected with "mode!=joined" and assigned to a sample name:
cat batra_ali_rmjoined_assigned_* > batrALL_ali_rmjoined_assigned.fq
# to be able to visualize the number of reads in function of the sample:
bsub <<< "obistat -c sample batrALL_ali_rmjoined_assigned.fq"

### 6. DEreplicate reads into uniq sequences
# Ran the command
bsub -n 4 -M 16000000 -q "long" -J "obiuniq" 'obiuniq -m sample batrALL_ali_rmjoined_assigned.fq > batrALL_alijoined_assigned_uniq.fq'


### 7. Remove singleton

#to have a look at the number of singleton, I used the following command:
obistat -c count batrALL_alijoined_assigned_uniq.fq | sort -nk1 | head -20

count             count     total
1                424419    424419
2                 48963     97926
3                 22089     66267
4                 13121     52484
5                  8912     44560
6                  6411     38466
7                  5013     35091
8                  3934     31472
9                  3401     30609
10                 2710     27100
11                 2383     26213
12                 2032     24384
13                 1826     23738
14                 1614     22596
15                 1381     20715
16                 1281     20496
17                 1140     19380
18                  977     17586
19                  973     18487


#I first chose to removed uniquely singleton, to be as conservative as possible:
obigrep -p "count>1" batrALL_alijoined_assigned_uniq.fq > batrALL_alijoined_assigned_uniq_C1.fq

# Following the OBITools tutorial (Wolf data): clean the sequences from PCR/sequencing errors
#bsub -n 4 -M 16000000 -q "long" -J "Obiclean" 'obiclean -s merged_sample -r 0.05 -H  batrALL_alijoined_assigned_uniq_C1.fq > batrALL_alijoined_assigned_uniq_C1_clean.fq'
#pas bien d'aggréger les séquences avant de les assigner!

### 8. Matching sequences to ref database
bsub -n 4 -M 16000000 -q "long" -J "Ecotag" 'ecotag -d ../Vert/embl_last_vert -R Batr_3m_final.fa -m 0.95 batrALL_alijoined_assigned_uniq_C1.fq > batrALL_alijoined_assigned_uniq_C1_tags.fq'

obigrep -p 'taxid!=1' batrALL_alijoined_assigned_uniq_C1_tags.fq | obigrep -p '(best_identity.values())[0]==1' | obistat -c species_name | sort -nk4

None                                  5  47810685
species_name                        count   total
Panurus biarmicus                     1      7644
Pseudacris maculata                   1      9608
Fringilla coelebs                     1     13241
Phalacrocorax carbo                   1     42430
Fulica atra                           1    130896
Abramis brama                         1    147056
Perca fluviatilis                     1    415960
Netta rufina                          1    486105
Esox lucius                           1    614072
Rutilus rutilus                       1    647803
Bufo bufo                             1    783929
Homo sapiens                         12   1362013
Pelophylax ridibundus                 1   1783757
Lepomis gibbosus                      1   1845215
Rana temporaria                       3   2481399
Scardinius erythrophthalmus           1   7645316
Rhodeus amarus                        1   8768630
Tinca tinca                           1  13992160
Rana arvalis                          1  32330034

#Retry with the database of Eduard!
obiconvert -t TAXO --ecopcrdb-output=taxonomy # ../
bsub -n 4 -M 16000000 -J "Ecotag" 'ecotag -t ../taxonomy -R ref_database/batr01.uniq.fastq -m 0.95 batrALL_alijoined_assigned_uniq_C1.fq > batrALL_alijoined_assigned_uniq_C1_tags_Eduard.fq'

obigrep -p 'taxid!=1' batrALL_alijoined_assigned_uniq_C1_tags_Eduard.fq | obigrep -p '(best_identity.values())[0]==1' | obistat -c species_name | sort -nk4


### Adding missing species in the ref_database
#1. create a header corresponding to : >merged_taxid={70019:1}; species_name=Pelophylax saharicus; family_name=Ranidae; taxid=70019; reverse_match=ACACCGCCCGTCACCCT
#2. copy past the amplicon in lower-case letter.

### Retry with my restreinte database:
bsub -n 4 -M 16000000 -q "long" -J "Ecotag" 'ecotag -d ../Vert/embl_last_vert -R ref_database/Batr_3m_final.fa -m 0.95 batrALL_alijoined_assigned_uniq_C1.fq > batrALL_alijoined_assigned_uniq_C1_tags.fq'

### Retry with Eduard databases
nohup ecotag -d ../TAXOEduard/obiTAXO_17.10.2018 -R ../ref_db_eduard/db_batr01_addsp.fasta -m 0.95 batrALL_alijoined_assigned_uniq_C1.fq > batrALL_alijoined_assigned_uniq_C1_tagsEduard_notrest.fq 2>stdoutEdbatr &

nohup ecotag -d ../TAXOEduard/obiTAXO_17.10.2018 -R ../ref_db_eduard/db_batr01.restr_addsp.fasta -m 0.95 batrALL_alijoined_assigned_uniq_C1.fq > batrALL_alijoined_assigned_uniq_C1inseALL_alijoined_assigned_uniq_C1_tagsEduard_restr.fq 2>stdoutEdbatrRestr95 &


"batrALL_alijoined_assigned_uniq_C1.fq.restr.tag" name of the final file file once tag are added

### Remove PCR/sequencing errors
obiclean -r 0.25 -H batrALL_alijoined_assigned_uniq_C1.fq.restr.tag > batrALL_alijoined_assigned_uniq_C1_restr.tag_clean.fq


####### LAST VERSION 2022
#assignement data base with full db added L. helveticus and P. bergeri

ecotag -d Vert19/embl_last_vert -R ref_database/db_batr01_addsp.fasta 6_Batr_clean_all.fq > 7_Batr_taxid_asigne.fq

obigrep -p 'taxid!=1' 7_Batr_taxid_asigne.fq | obigrep -p '(best_identity.values())[0]==1' | obistat -c species_name | sort -nk4

None                                 27  48081039
Pseudacris n. sp. ECM-2007            2        25
species_name                                                               count           total
Achromobacter sp. FZ97                1        42
Simocephalus sp. EIZ-2012             1        22
Variovorax sp. Tibet-S862             1       824
uncultured Fusobacterium sp.          1         2
uncultured Synechococcus sp.          1      3726
Micromys minutus                      1         2
Phacus ranula                         1         2
Synedra fragilaroides                 1         2
Acanthobrama persidis                 1         4
Phacus orbicularis                    1         4
Stauroneis constricta                 1         4
Asarcornis scutulata                  1         7
Fringilla montifringilla              1         7
Pseudacris brachyphona                1         8
Pelophylax kurtmuelleri               1         9
Japyx solifugus                       1        12
Pseudacris feriarum                   1        21
Pseudacris brimleyi                   1        23
Rhodeus sericeus                      1        23
Alburnus mossulensis                  1        24
Eunotia naegelii                      1        24
Acinetobacter sp.                     1        33
Anas crecca                           1        35
Ardea purpurea                        1        38
Bos taurus                            4        47
Podiceps cristatus                    1        64
Circus aeruginosus                    1        68
Pan troglodytes                       2        72
Fulica americana                      1        81
Arachis hypogaea                      1        98
Cervus elaphus                        1       162
Gorilla gorilla                       1       179
Capreolus capreolus                   1       230
Perca flavescens                      1       401
Synedra hyperborea                    1       412
Gymnoris dentata                      1       450
Anas platyrhynchos                    1       490
Pupa strigosa                         3       511
Perca schrenkii                       1       607
Neomys fodiens                        1       881
Gallus gallus                         2       887
Ulnaria acus                          1      1156
Squalius lepidus                      1      1213
Gasterosteus aculeatus                1      1416
Cygnus olor                           1      1633
Bufo verrucosissimus                  1      1807
Anas chathamica                       1      2523
Esox flaviae                          1      2747
Xenopus tropicalis                    1      4390
Lissotriton vulgaris                  1      6429
Panurus biarmicus                     1      7644
Alburnus alburnus                     1     10191
Hyla arborea                          1     11783
Fringilla coelebs                     1     13241
Alburnus tarichi                      1     14850
Sus scrofa                            5     16057
Gobio gobio                           1     18117
Rana pyrenaica                        1     25778
Pseudacris maculata                   3     27942
Pseudochondrostoma polylepis          1     33926
Phalacrocorax carbo                   1     42430
Cairina moschata                      1     48191
Fulica atra                           1    130896
Abramis brama                         1    147056
Acorus calamus                        1    157825
Perca fluviatilis                     1    415960
Netta rufina                          1    486105
Esox lucius                           2    614528
Rutilus rutilus                       2    647884
Homo sapiens                         12    694658
Bufo bufo                             2    784384
Pelophylax ridibundus                 1   1783757
Lepomis gibbosus                      2   1848076
Pelophylax bergeri                    1   1860967
Rana temporaria                       4   2485071
Scardinius erythrophthalmus           1   7645316
Rhodeus amarus                        1   8768630
Tinca tinca                           1  13992160
Rana arvalis                          1  32330034

obiclean -r 0.25 -H 7_Batr_taxid_asigne.fq > 8_Batr_taxid_dbassigned_PCRcleaned.fq

obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
  --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
  --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
  --delete-tag=obiclean_head --delete-tag=taxid_by_db --delete-tag=obiclean_headcount \
  --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=order_name \
   8_Batr_taxid_dbassigned_PCRcleaned.fq > \
  9_Batr_woannotate.fq

  obitab -o 9_Batr_woannotate.fq > 10_Batr_cleaned.tab
