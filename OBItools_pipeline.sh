  # Curated Version that can be published 2022

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

### 3. Remove unaligned sequences record
#I performed obigrep by removing joined sequences and not by setting a threshold
for i in $(ls batra_*.fq | sed 's/batra_//g' | sed 's/.fq//g'); do bsub <<< "obigrep -p 'mode!=\"joined\"' batra_${i}.fq > batra_${i}_ali_rmjoined.fq"; done

### 4. Assign marker/sample combination
#for sequences selected with the parameter 'mode!="joined"'
for i in $(ls batra_* | sed 's/batra_//g' | sed 's/_ali_rmjoined.fq//g'); do bsub <<< "ngsfilter -t generic_ngsfilter_Batr01.txt -u unidentified_Batra.fq batra_${i}_ali_rmjoined.fq > batra_ali_rmjoined_assigned_${i}.fq"; done

# We then needed to concatenated the sequences alligned, selected with "mode!=joined" and assigned to a sample name:
cat batra_ali_rmjoined_assigned_* > batrALL_ali_rmjoined_assigned.fq

### 6. Dereplicate reads into uniq sequences
# Ran the command
bsub -n 4 -M 16000000 -q "long" -J "obiuniq" 'obiuniq -m sample batrALL_ali_rmjoined_assigned.fq > batrALL_alijoined_assigned_uniq.fq'


### 7. assignement data base with full db added L. helveticus and P. bergeri

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
