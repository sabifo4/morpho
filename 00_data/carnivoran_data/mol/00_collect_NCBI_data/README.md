# Format NCBI sequences

## Get sequences in one line
In order to get the downloaded gene sequences for each specimen in one line, 
we ran the following loop using the perl script `one_line_fasta.pl`

```bash 
mkdir NCBI_seqs_download_onelinefasta

for i in NCBI_seqs_download/*/*fasta
do 
fname=$( echo $i | sed 's/..*\///' | sed 's/\.fasta//' )
dir=$( echo $i | sed 's/NCBI\_seqs\_download\///' | sed 's/\/..*//' )
mkdir -p NCBI_seqs_download_onelinefasta/$dir
perl one_line_fasta.pl $i 
mv NCBI_seqs_download/$dir/$fname"_one_line.fasta" NCBI_seqs_download_onelinefasta/$dir
done

tar -zcvf NCBI_seqs_download.tar.gz NCBI_seqs_download/
rm -r NCBI_seqs_download/
```

This generates a new folder called `NCBI_seqs_download_onelinefasta` with 
the sequences of each gene for each species written in one line. The directory 
`NCBI_seqs_download` is tarred and then deleted to save space.

## Generate dummy files for missing genes under species directories
Then, we checked which genes were available for each species. If the gene 
was not available, then a file including the key word "nothing" in its name 
was created, in which only a fasta header with the species name was 
included. The following code was used for this purpose: 

```bash 
# Move to directory with dirs with NCBI sequences
cd NCBI_seqs_download_onelinefasta

# Create global arrays to define genes
vars=("COX1" "COX2" "COX3" "CYTB" "NADH1" "NADH2" "NADH3" "NADH4" "NADH4L" "NADH5" "ATP6" "ATP8")
visited=()

# Loop over each directory to find which genes are missing in 
# each species 
for j in *
do

for i in $j/*fasta
do
dir=$( echo $i | sed 's/\/..*//' )
spname=$( echo $i | sed 's/..*\///' | sed 's/\_/\-/2' | sed 's/\-..*//' )
#echo $spname 
tag=$( echo $i | sed 's/..*\///' | sed 's/'${spname}'\_//' | sed 's/\_..*//' )
#echo $tag  
visited+=( $tag )
done

dif=(`echo ${visited[@]} ${vars[@]} | tr ' ' '\n' | sort | uniq -u`)
for missing in "${dif[@]}"
do
#echo $missing
name=$( echo $spname | sed 's/\_/\ /' )
printf "Gene "$missing" not available for species "$spname"\n"
printf "Gene "$missing" not available for species "$spname"\n" >> log_missing_genes_per_species.txt
printf ">"$missing"_blank "$name"\n" > $j/$spname"_"$missing"_nothing.fasta"
done

visited=()

done 

``` 

This creates a log file called `log_missing_genes_per_species.txt`, which 
contains a description of the genes not available in specific species.
