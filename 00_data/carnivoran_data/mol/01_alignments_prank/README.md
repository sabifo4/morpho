# Generate gene alignments 

## Generate input file for PRANK 

Once we generated a file for each gene under each species directory - including a dummy one 
for when the gene was not available - we generated an input file for 
[`PRANK v.150803`](http://wasabiapp.org/software/prank/) so we could generate the gene alignments. 

First, we needed to generate a file that contained all the available gene sequences - this means 
that not all specimens were included in this file as not all genes were available for all the species.

```bash 
# The commands used in this snippet assume that you run this
# code inside `01_alignments_prank`

# Declare global array in which we define the name of all the 
# genes available
vars=("COX1_" "COX2_" "COX3_" "CYTB_" "NADH1_" "NADH2_" "NADH3_" "NADH4_" "NADH4L_" "NADH5_" "ATP6_" "ATP8_")

# Loop over the filtered sequences 12 times to check if, for each species,
# the gene is available

for j in `seq 0 11`
do

for i in ../00_collect_NCBI_data/NCBI_seqs_download_onelinefasta/*/*${vars[j]}*.fasta
do

# If this is not a dummy sequence,
# then please append it to the file that will be 
# used to generate the gene alignment
if [[ ! $i =~ "nothing" ]]
then
lines=$(wc -l $i)
num=$(echo $lines | sed 's/ ..*//')
if [[ $num == 1 ]]
then
printf "Something is wrong, you should not be here...\n"
else
cat $i >> ${vars[j]}mtgenes.fasta
printf "\n" >> ${vars[j]}mtgenes.fasta
fi
fi

done

done

# Create one directory per gene in which each 
# corresponding file will be moved. 
# PRANK will be used later to proceed with the alignments 
mkdir 00_alignments 

for i in ${vars[@]}
do 

name=$( echo $i | sed 's/\_//' )
mkdir -p 00_alignments/$name 
mv $i*mtgenes*fasta 00_alignments/$name 

done 
```

This generates one directory per gene with the input file for `PRANK` 
inside the `00_alignments` directory.

## Run `PRANK` 

We used the following code to run `PRANK` and generate the 
12 gene alignments:

```bash 
# The commands used in this snippet assume that you run this
# code inside `01_alignments_prank`

for i in 00_alignments/*/*fasta 
do

namedir=$( echo $i | sed 's/00\_alignments\///' | sed 's/\/..*//' )
file=$( echo $i | sed 's/..*\///' )
cd 00_alignments/$namedir
prank -d=$file +F -o=$namedir &
cd ../..

done
```

Then, we renamed the files output by `PRANK` and identified the length of the 
sequences so we could add the corresponding amount of gaps to the sequences for 
the rest of species for which the gene was not available:

```bash 
# The commands used in this snippet assume that you run this
# code inside `01_alignments_prank`

# Change the name of the output file to 
# $name_gene"_prankaln.fasta"
for i in 00_alignments/*/*best.fas
do
name=$( echo $i | sed 's/00\_Alignments\///' | sed 's/\/..*//' )
path=$(echo $i | sed 's/\/'$name'\.best\.fas//')
mv $i $path/$name"_prankaln.fasta"
done

# Now get each sequence in the alignment in one line:
for i in 00_alignments/*/*prankaln*
do
perl ../00_collect_NCBI_data/one_line_fasta.pl $i
done

# Count length of genes
# First grep the pattern and get only this line '/\>/{n;p;}'
# As we want only the first coincidence, we have to add " ; "
# to stop at first occurrence with "q": /\>/q"
for i in 00_alignments/*/*one_line*
do
name=$( echo $i | sed 's/00\_Alignments\///' | sed 's/\/..*//' )
path=$(echo $i | sed 's/\/'$name'..*fasta//')
file=$(sed -n '/\>/{n;p;} ; /\>/q' $i | sed 's/\n//')
num=$(echo -n $file | wc -m)
echo The sequence length for $i is $num
echo The sequence length for $i is $num >> log_seqlength_perspecies.txt 
## Print the dashes $num times
NAs=$(printf -- '-%.0s' `seq 1 $num`)
#echo $NAs
#printf "\n"
echo $NAs > $path/$name/$name"_NAs_to_add.txt"
done 
```

After that, we ran the next code to (i) rename the fasta headers (compatible PHYLIP format 
to be read by `MCMCtree`) and (ii) add the gap sequences of the species for which the 
genes were not available:

```bash 
# The commands used in this snippet assume that you run this
# code inside `01_alignments_prank`

for i in 00_alignments/*/*prankaln_one_line.fasta 
do

name=$( echo $i | sed 's/\.fasta/\_renamed\.fasta/' )
gene=$( echo $i | sed 's/00\_alignments\///' | sed 's/\/..*//' )
path=$(echo $i | sed 's/\/'$gene'..*fasta//')
cp $i $name

sed -i 's/..*Canis..*/\>Can\_lyc/' $name
tmp=`grep 'Can_lyc' $name`
if [[ ! $tmp ]]
then 
printf ">Can_lyc\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

sed -i 's/..*Cuon..*/\>Cuo\_alp/' $name
tmp=`grep 'Cuo_alp' $name`
if [[ ! $tmp ]]
then 
printf ">Cuo_alp\n" >> $name 
cat $path/$gene*NAs_to_add.txt >> $name
fi

sed -i 's/..*Speothos..*/\>Spe\_ven/' $name
tmp=`grep 'Spe_ven' $name`
if [[ ! $tmp ]]
then 
printf ">Spe_ven\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

sed -i 's/..*Vulpes..*/\>Vul\_vul/' $name
tmp=`grep 'Vul_vul' $name`
if [[ ! $tmp ]]
then 
printf ">Vul_vul\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

sed -i 's/..*Ailurus..*/\>Ail\_ful/' $name
tmp=`grep 'Ail_ful' $name`
if [[ ! $tmp ]]
then 
printf ">Ail_ful\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

sed -i 's/..*Nandinia..*/\>Nan\_bin/' $name
tmp=`grep 'Nan_bin' $name`
if [[ ! $tmp ]]
then 
printf ">Nan_bin\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

sed -i 's/..*Ursus..*/\>Urs\_ame/' $name
tmp=`grep 'Urs_ame' $name`
if [[ ! $tmp ]]
then 
printf ">Urs_ame\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

sed -i 's/..*Cerdocyon..*/\>Cer\_tho/' $name
tmp=`grep 'Cer_tho' $name`
if [[ ! $tmp ]]
then 
printf ">Cer_tho\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

sed -i 's/..*Otocyon..*/\>Oto\_meg/' $name
tmp=`grep 'Oto_meg' $name`
if [[ ! $tmp ]]
then 
printf ">Oto_meg\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

sed -i 's/..*Paradoxurus..*/\>Par\_her/' $name
tmp=`grep 'Par_her' $name`
if [[ ! $tmp ]]
then 
printf ">Par_her\n" >> $name 
cat $path/$gene/*NAs_to_add.txt >> $name
fi

done
```

The last step was to generate a partitioned alignment. We first partitioned each gene into two partitions: (i) first and second codon 
positions (1st+2nd CPs) and (ii) third codon positions (3rd CPs).
After that, we concatenated all 1st+2nd CPs of each gene in one partition and the 3rd CPs 
of each gene in a second partition. We used the [`fasta-phylip-partitions` pipeline](https://github.com/sabifo4/fasta-phylip-partitions) 
in order to generate the partitioned alignment by gene and codon positions: 

```bash 
# The commands used in this snippet assume that you run this
# code inside `01_alignments_prank`

# 1. Copy all renamed alignments into another directory 
#    in which all partitioned alignments will be generated 
cp 00_alignments/*/*renamed.fasta 00_alignments 
mkdir -p 01_partitioned_alignments/00_partitioning_step 
mv 00_alignments/*renamed.fasta 01_partitioned_alignments/00_partitioning_step  

# 2. Create the `species_names.txt` file needed by the 
#    `fasta-phylip-partitions` pipeline 
grep '>' 01_partitioned_alignments/00_partitioning_step/ATP6_prankaln_one_line_renamed.fasta | sed 's/>//g' > 01_partitioned_alignments/00_partitioning_step/species_names.txt
Run_tasks.sh  01_partitioned_alignments/00_partitioning_step/ mol_canids partY

# 3. Get concatenated 1st+2nd CPs partition and 3rd CPs partition in one 
#    unique file and copy the separated 1st+2nd CPs and 3rd CPs
mkdir -p 01_partitioned_alignments/01_partitioned_alignment 
cat 01_partitioned_alignments/00_partitioning_step/phylip_format/02_concatenated_alignments/part12/part12_mol_canids_concat.aln > 01_partitioned_alignments/01_partitioned_alignment/mit_12.3.aln
printf "\n\n" >> 01_partitioned_alignments/01_partitioned_alignment/mit_12.3.aln
cat 01_partitioned_alignments/00_partitioning_step/phylip_format/02_concatenated_alignments/part3/part3_mol_canids_concat.aln >> 01_partitioned_alignments/01_partitioned_alignment/mit_12.3.aln

cp 01_partitioned_alignments/00_partitioning_step/phylip_format/02_concatenated_alignments/part12/part12_mol_canids_concat.aln  01_partitioned_alignments/01_partitioned_alignment/mit_12.aln

cp 01_partitioned_alignments/00_partitioning_step/phylip_format/02_concatenated_alignments/part3/part3_mol_canids_concat.aln  01_partitioned_alignments/01_partitioned_alignment/mit_3.aln

```

The final alignments that will be used in `MCMCtree` for Bayesian divergence times estimation have been 
saved inside `01_partitioned_alignments/01_partitioned_alignment/`. They have also been 
saved in the [`alignments`](https://github.com/sabifo4/morpho/blob/master/alignments) directory.