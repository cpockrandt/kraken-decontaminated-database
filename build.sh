#!/bin/bash -e

# YOU SHOULD ADJUST THIS
#########################################################
THREADS=16
PATH_TO_KMC="./kmc" # path to binary of kmer counter KMC3
#########################################################

CONTERMINATOR_FILE="./conterminator_refseq_and_genbank"
PATH_TO_CONTERMINATOR_HELPER="./build/bin/decontam_conterminator"
PATH_TO_CONTAMINATION_CLEANUP="./build/bin/decontam_kmer"

function check_command {
  if ! [ -x "$(command -v $1)" ]; then
    echo "Error: $1 is not installed." >&2
    exit 1
  fi
}

function krakenuniq_download_wrapper {
  if [[ "$1" =~ ^GCA ]]
  then
    krakenuniq-download --db genbank genbank/$2/assembly_accession=$1 --min-seq-len 50000
    krakenuniq-download --db genbank genbank/$2/assembly_accession=$1 --fna=rna_from_genomic || true
  else
    krakenuniq-download --db refseq refseq/$2/assembly_accession=$1 --min-seq-len 50000
    krakenuniq-download --db refseq refseq/$2/assembly_accession=$1 --fna=rna_from_genomic || true
  fi
}

check_command "parallel"
check_command "krakenuniq-download"
check_command "dustmasker"
check_command "wget"
check_command "gunzip"

echo "This script runs GNU parallel. Please check out 'parallel --citation'"

# The following commands download viruses, vector sequences, fungi, bacteria, archaea, selected eucaryotic pathogens and the human genome. Feel free to comment in or add anything you need
# NOTE: refseq and genbank data has to be downloaded into different folders (--db refseq and --db genbank)

echo "Downloading files ..."

krakenuniq-download --db refseq --threads ${THREADS} refseq/viral
krakenuniq-download --db refseq --threads ${THREADS} viral-neighbors
krakenuniq-download --db refseq --threads ${THREADS} contaminants
krakenuniq-download --db refseq --threads ${THREADS} refseq/bacteria refseq/archaea

krakenuniq-download --db refseq refseq/vertebrate_mammalian/Chromosome/species_taxid=9606
krakenuniq-download --db refseq refseq/vertebrate_mammalian/Chromosome/species_taxid=9606 --fna rna 

krakenuniq-download --db refseq refseq/fungi
krakenuniq-download --db refseq refseq/fungi --fna=rna_from_genomic # THIS WAS MISSING in florians script

krakenuniq-download --db refseq refseq/protozoa
krakenuniq-download --db refseq refseq/protozoa --fna=rna_from_genomic # THIS WAS MISSING in florians script

krakenuniq-download --db genbank genbank/fungi/Any/assembly_accession=GCA_002775015.1 --min-seq-len 20000 # Candida auris assembly by CDC
krakenuniq-download --db genbank genbank/fungi/Any/assembly_accession=GCA_002775015.1 --fna=rna_from_genomic || true # Candida auris assembly by CDC

krakenuniq-download --db genbank genbank/fungi/Any/assembly_accession=GCA_000497105.2 --min-seq-len 20000 # Candida glabrata assembly by CDC
krakenuniq-download --db genbank genbank/fungi/Any/assembly_accession=GCA_000497105.2 --fna=rna_from_genomic || true # Candida glabrata assembly by CDC

ASP_NIDULANS=GCA_000011425.1
BEAUVERIA_BASS=GCA_003337105.1

for AC in $ASP_NIDULANS $BEAUVERIA_BASS GCA_001599635.1 GCA_001599015.1 GCA_001600735.1 GCA_001598995.1 GCA_000696995.1 GCF_000151125.1 GCA_001511075.1 GCA_001445615.1 GCA_000812125.1 GCA_003336255.1 GCA_000151595.1 GCA_000982555.2 GCA_001077315.1 GCA_001046915.1 GCA_001430955.1 GCA_000697235.1 GCA_001566745.1 GCA_000697335.1 GCA_000786445.1 GCA_000835815.1 GCA_000149475.3 GCA_000855695.1 GCA_000835755.1 GCA_000697315.1 GCA_000697015.1 GCA_000743335.1 GCA_002572885.1 GCA_001402995.1 GCA_000723665.1 GCA_000697395.1 GCA_001275765.2 GCA_002551515.1 GCA_001264725.1 GCA_001264965.1 GCA_001264885.1 GCF_000150975.2 GCA_000697295.1 GCA_000587855.1 GCA_000697255.1 GCF_000226095.1 GCA_000611715.1 GCA_002080375.1 GCA_001599375.1 GCA_000697605.1 GCA_003055205.1 GCA_000697055.1 GCA_000817185.1 GCA_000812075.1 GCA_000820605.1 GCA_001630435.1 GCA_000696955.1 GCA_003025155.1 GCA_000151175.1 GCA_000616785.1 GCA_000616865.1 GCA_000151455.1 GCA_001651435.1 GCA_001752605.1 GCA_001752585.1 GCA_001752625.1 GCA_001752645.1 GCA_002798055.1 ; do 
  krakenuniq_download_wrapper $AC "fungi/Any"
done

ALTERNARIA_ALTERNATA=GCF_001642055.1
ASP_CLAVATUS=GCF_000002715.2
ASP_FLAVUS=GCF_000006275.2
ASP_FUMIGATUS=GCF_000002655.1
ASP_NIGER=GCF_000002855.3
ASP_ORYZ=GCF_000184455.2
ASP_TEREUS=GCF_000149615.1
CAND_ALBICANS=GCF_000182965.3
CAND_DUBL=GCF_000026945.1
CAND_ORTHO=GCF_000315875.1
CAND_TROP=GCF_000006335.3
CHAETOM_GLOB=GCF_000143365.1
CLAD_BANT=GCF_000835475.1
CLAD_CARR=GCF_000365165.1
COCC_IMM=GCF_000149335.2
COCC_POSA=GCF_000151335.2
CRYPTO_GATTII=GCF_000185945.1
CRYPTO_NEO_VAR_GRUB=GCF_000149245.1
CRYPTO_NEO_VAR_NEO=GCF_000091045.1
CRYPTO_NEO_VAR_NEO2=GCF_000149385.1
Enterocytozoon_bieneusi=GCA_000209485.1
Encephalitozoon_cuniculi=GCA_000091225.2
Encephalitozoon_hellem=GCF_000277815.2
Encephalitozoon_intestinalis=GCF_000146465.1

for AC in $ALTERNARIA_ALTERNATA $ASP_CLAVATUS $ASP_FLAVUS $ASP_FUMIGATUS $ASP_NIGER $ASP_ORYZ $ASP_TEREUS $CAND_ALBICANS $CAND_DUBL $CAND_ORTHO $CAND_TROP $CHAETOM_GLOB $CLAD_BANT $CLAD_CARR $COCC_IMM $COCC_POSA $CRYPTO_GATTII $CRYPTO_NEO_VAR_GRUB $CRYPTO_NEO_VAR_NEO $CRYPTO_NEO_VAR_NEO2 GCF_000230625.1 GCF_000836115.1 GCF_000835505.1 GCF_000835455.1 GCF_900079805.1 GCF_000149955.1 GCF_000149555.1 GCF_000149585.1 GCF_001417885.1 GCF_000149685.1 GCF_000181695.1 GCF_001278385.1 GCF_000349305.1 GCF_000149425.1 GCF_000151145.1 GCF_000150735.1 GCF_000150705.2 GCF_001477535.1 GCF_000968595.1 GCF_000835555.1 GCF_000439145.1 GCF_002708625.1 GCF_000146045.2 GCF_000732125.1 GCF_000961545.1 GCF_000001985.1 GCF_000151425.1 GCF_000151505.1 GCF_000293215.1 GCF_000836295.1 GCF_000002525.2 $Enterocytozoon_bieneusi $Encephalitozoon_cuniculi $Encephalitozoon_hellem $Encephalitozoon_intestinalis; do
  krakenuniq_download_wrapper $AC "fungi/Any"
done

Ehistolytica_2016=GCF_000208925.1
Entamoeba_dispar=GCF_000209125.1
Giardia_intestinalis=GCF_000002435.1
Cryptosporidium_parvum=GCF_000165345.1
Cyclospora_cayetanensis=GCF_002999335.1
Plasmodium_falciparum=GCF_000002765.4
Plasmodium_falciparum2=GCA_900632045.1
Plasmodium_vivax=GCF_000002415.2
Plasmodium_ovale=GCA_900090025.2
Plasmodium_ovale=GCA_900090035.2
Plasmodium_malariae=GCF_900090045.1
Plasmodium_knowlesi=GCF_000006355.1
Babesia_divergens=GCA_001077455.2
Babesia_microti=GCF_000691945.2
Trypanosoma_brucei_gambiense=GCF_000210295.1
Trypanosoma_brucei_brucei=GCF_000002445.2
Trypanosoma_cruzi=GCF_000209065.1
Leishmania_major_strain_Friedlin=GCF_000002725.2
Leishmania_infantum=GCA_003020905.1
Leishmania_donovani=GCA_003719575.1
Leishmania_panamensis=GCF_000755165.1
Leishmania_braziliensis=GCF_000002845.2
Leishmania_arabica=GCA_000410695.2
Sarcocystis_neurona=GCA_000875885.1
Toxoplasma_gondii_ME49=GCF_000006565.2
Naegleria_gruberi=GCF_000004985.1
Naegleria_lovaniensis=GCA_003324165.1
Naegleria_fowleri=GCA_000499105.1
Acanthamoeba_castellanii_Neff=GCF_000313135.1
Acanthamoeba_polyphaga=GCA_001567625.1
Acanthamoeba_lenticulata=GCA_002179805.1
Acanthamoeba_astronyxis=GCA_000826245.1
Acanthamoeba_culbertsoni=GCA_000826265.1
Acanthamoeba_divionensis=GCA_000826405.1
Acanthamoeba_lugdunensis=GCA_000826425.1
Entamoeba_moshkovskii=GCA_002914575.1
Trichomonas_vaginalis=GCA_002891335.1

for AC in $Ehistolytica_2016 $Entamoeba_dispar $Giardia_intestinalis $Cryptosporidium_parvum $Cyclospora_cayetanensis $Plasmodium_falciparum $Plasmodium_falciparum2 $Plasmodium_vivax $Plasmodium_ovale $Plasmodium_ovale $Plasmodium_malariae $Plasmodium_knowlesi $Babesia_divergens $Babesia_microti $Trypanosoma_brucei_gambiense $Trypanosoma_brucei_brucei $Trypanosoma_cruzi $Leishmania_major_strain_Friedlin $Leishmania_infantum $Leishmania_donovani $Leishmania_panamensis $Leishmania_braziliensis $Leishmania_arabica $Sarcocystis_neurona $Toxoplasma_gondii_ME49 $Naegleria_gruberi $Naegleria_lovaniensis $Naegleria_fowleri $Acanthamoeba_castellanii_Neff $Acanthamoeba_polyphaga $Acanthamoeba_lenticulata $Acanthamoeba_astronyxis $Acanthamoeba_culbertsoni $Acanthamoeba_divionensis $Acanthamoeba_lugdunensis $Entamoeba_moshkovskii $Trichomonas_vaginalis; do
  krakenuniq_download_wrapper $AC "protozoa/Any"
done

## Metazoa
Enterobius_vermicularis=GCA_900576705.1
Necator_americanus=GCF_000507365.1
Trichuris_trichiura=GCA_000613005.1
Strongyloides_sterncoralis=GCF_001040885.1
Brugia_malayi=GCA_000002995.5
Brugia_timori=GCA_900618025.1
Loa_loa=GCF_000183805.2
Onchocerca_volvulus=GCA_000499405.2
Toxocara_canis=GCA_000803305.1
Trichinella_spiralis=GCA_000181795.3
Taenia_saginata=GCA_001693075.2
Taenia_solium=GCA_001870725.1
Hymenolepis_nana=GCA_900617975.1
Echinococcus_granulosus=GCF_000524195.1
Echinococcus_multilocularis=GCA_000469725.3
Schistosoma_mansoni=GCF_000237925.1
Schistosoma_japonicum=GCA_000151775.1
Schistosoma_haematobium=GCF_000699445.1
Clonorchis_sinensis=GCA_003604175.1
Opisthorchis_viverrini=GCF_000715545.1
Opisthorchis_felineus=GCA_004794785.1
Fasciola_hepatica=GCA_900302435.1
Anisakis_simplex=GCA_900617985.1

for AC in $Enterobius_vermicularis $Necator_americanus $Trichuris_trichiura $Strongyloides_sterncoralis $Brugia_malayi $Brugia_timori $Loa_loa $Onchocerca_volvulus $Toxocara_canis $Trichinella_spiralis $Taenia_saginata $Taenia_solium $Hymenolepis_nana $Echinococcus_granulosus $Echinococcus_multilocularis $Schistosoma_mansoni $Schistosoma_japonicum $Schistosoma_haematobium $Clonorchis_sinensis $Opisthorchis_viverrini $Opisthorchis_felineus $Fasciola_hepatica $Anisakis_simplex; do
  krakenuniq_download_wrapper $AC "invertebrate/Any"
done

### Merge refseq and genbank folders into new folder "database"
echo "Merging refseq and genbank folders into a single one"
rsync -a genbank/library refseq && rm -rf genbank && mv refseq database

### Remove low complexity with dustmasker
echo "Running duskmasker on all files"
find ./database/library \( -name "*.fna" -o -name "*.fasta" \) | parallel --will-cite -j ${THREADS} "/ccb/sw/bin/dustmasker -in {} -infmt fasta -level 20 -outfmt fasta | sed '/^>/! s/[agct]/N/g' > {}.dusted && mv {}.dusted {}"

### Removing contamination with Conterminator by Martin Steinegger (in a simplified way: when there is an unspaced alignment of a certain length between two species of different kingdoms)

echo "Removing contamination with conterminator (by Martin Steinegger)"

# download contamination information
wget -O - -o /dev/null ftp://ftp.ccb.jhu.edu/pub/data/conterminator/genbank.gz | gunzip > ${CONTERMINATOR_FILE}
wget -O - -o /dev/null ftp://ftp.ccb.jhu.edu/pub/data/conterminator/refseq.gz | gunzip >> ${CONTERMINATOR_FILE}

# extract sequence ids from dna data
echo "id\tdb\tdomain\tassembly-level\tseqlen" > dna_ids
for d in $(ls -d database/*/*/*/); do
  #echo -e "$d"
  find $d -type f -name '*.fna' -print | xargs awk -F ' ' '($0 ~ /^>/ && !(FILENAME ~ /rna_from/)) { if (NR > 1) { print len } len = 0; n = split(FILENAME, a, "/"); printf "%s\t%s\t%s\t%s\t%s\t", substr($1, 2), a[1], a[3], a[4], FILENAME } (!($0 ~ /^>/) && !(FILENAME ~ /rna_from/)) { len += length($0) } END { print len }' >> dna_ids
done

# extract sequence ids from rna data (will miss some tRNAs)
echo "id\tdb\tdomain\tassembly-level\tseqlen" > rna_ids
for d in $(ls -d database/*/*/*/); do
  #echo -e "$d"
  if ls $d*rna_from_genomic.fna > /dev/null 2>&1; then # NOTE: nullglob can not be set!
    awk -F ' ' '($0 ~ /^>/) { if (NR > 1) print len; len = 0; n = split(FILENAME, a, "/"); rnapos = index($0, "_mrna_"); id=substr($0, 6, rnapos - 6); printf "%s\t%s\t%s\t%s\t%s\t", id, a[1], a[3], a[4], FILENAME } !($0 ~ /^>/) { len += length($0) } END { print len }' $d*rna_from_genomic.fna >> rna_ids
  fi

  if ls $d*rna.fna > /dev/null 2>&1; then
    awk -F ' ' '($0 ~ /^>/) { if (NR > 1) print len; len = 0; n = split(FILENAME, a, "/"); rnapos = index($0, "_mrna_"); id=substr($0, 6, rnapos - 6); printf "%s\t%s\t%s\t%s\t%s\t", id, a[1], a[3], a[4], FILENAME } !($0 ~ /^>/) { len += length($0) } END { print len }' $d*rna.fna >> rna_ids
  fi
done

# extract sequence ids from dna viral neighbors
echo "id\tseqlen" > dna_viral_neighbors_ids
find datebase/library/viral/Neighbors -type f -name '*.fasta' -print | xargs awk -F ' ' '($0 ~ /^>/) { if (NR > 1) { print len } len = 0; printf "%s\t%s\t", substr($1, 2), FILENAME } !($0 ~ /^>/) { len += length($0) } END{ print len }' >> dna_viral_neighbors_ids

# compute overlaps with conterminator data
awk -F$'\t' 'FNR==NR{ f[$2]=$0; next } { if ($1 in f) { print $0 "\t" f[$1]} }' ${CONTERMINATOR_FILE} dna_ids > dna_ids_overlap
awk -F$'\t' 'FNR==NR{ f[$2]=$0; next } { if ($1 in f) { print $0 "\t" f[$1]} }' ${CONTERMINATOR_FILE} rna_ids > rna_ids_overlap
awk -F$'\t' 'FNR==NR{ f[$2]=$0; next } { if ($1 in f) { print $0 "\t" f[$1]} }' ${CONTERMINATOR_FILE} dna_viral_neighbors_ids > dna_viral_neighbors_ids_overlap

# remove contamination
while read p; do
        # triplets=$(awk -v path="${p}" -F$'\t' '($5 == path) { printf "%s,%s,%s ", $1, $11, $12 }' dna_ids_overlap)
        triplets=$(awk -v path="${p}" -F$'\t' '($5 == path) { printf "%s %s %s\n", $1, $11, $12 }' dna_ids_overlap)
        p="${p}"; # TODO: if it does not work, prepend: ${PATH_TO_THIS_DIR} - also fix it in the next two blocks
        #echo "$p $triplets";
        ${PATH_TO_CONTERMINATOR_HELPER} -F "${p}" -C <(echo $triplets)
        mv "${p}.tmp.fa" "${p}"
done < <(cut -f5 dna_ids_overlap | sort -u)

while read p; do
        #triplets=$(awk -v path="${p}" -F$'\t' '($5 == path) { printf "%s,%s,%s ", $1, $11, $12 }' rna_ids_overlap)
        triplets=$(awk -v path="${p}" -F$'\t' '($5 == path) { printf "%s %s %s\n", $1, $11, $12 }' rna_ids_overlap)
        p="${p}";
        #echo "$p $triplets";
        ${PATH_TO_CONTERMINATOR_HELPER} -F "${p}" -C <(echo $triplets) -R
        mv "${p}.tmp.fa" "${p}"
done < <(cut -f5 rna_ids_overlap | sort -u)

while read p; do
        #triplets=$(awk -v path="${p}" -F$'\t' '($5 == path) { printf "%s,%s,%s ", $1, $11, $12 }' dna_viral_neighbors_ids_overlap)
        triplets=$(awk -v path="${p}" -F$'\t' '($5 == path) { printf "%s %s %s\n", $1, $11, $12 }' dna_viral_neighbors_ids_overlap)
        p="${p}";
        #echo "$p $triplets";
        ${PATH_TO_CONTERMINATOR_HELPER} -F "${p}" -C <(echo $triplets) -R
        mv "${p}.tmp.fa" "${p}"
done < <(cut -f5 rna_ids_overlap | sort -u)

# Removing human contamination (building a kmer-database on the human genome and masking everything that matches it. removes sequences and genomes if they are too contaminated with human DNA)

echo "Removing human contamination with kmers"

find ./database/library/vertebrate_mammalian -iname '*homo*sapiens*.fna' -exec readlink -f {} \; > kmc_human_files # put all paths to human dna/rna into a files

${PATH_TO_KMC} -k31 -ci1 -cs255 -cx3294967295 -t${THREADS} -fm @kmc_human_files kmc_human_db /tmp # index all the human dna/rna (all 31-mers)

# have to remove the human genome temporarily, otherwise it will be replaced with Ns entirely

for d in $(find ./database/library/vertebrate_mammalian/Chromosome -iname '*homo*sapiens*.fna'); do
  mv $d ./
done

# mask human dna in all genomes with Ns

${PATH_TO_CONTAMINATION_CLEANUP} -l ./contamination-cleanup.log  -t ${THREADS} --remove-seq-threshold 0.1 --remove-file-threshold  0.1 ./database/library ./kmc_human_db

# move human files back

for d in $(find . -maxdepth 1 -iname '*homo*sapiens*.fna'); do
  mv $d ./database/library/vertebrate_mammalian/Chromosome/
done

# Building kraken1 / kraken1uniq database

echo "Finally build kraken database"

krakenuniq-build --build --db database --threads ${THREADS}
