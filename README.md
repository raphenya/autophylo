# The Automatic Phylogentic Tree Builder

The **Auto**matic **Phylo**gentic (AutoPhylo) is used to generate phylogentic trees automatically by sampling the selected database (e.g., NCBI nr) and performs all tasks associated with traditional phylogentic tree building which includes trimming, dropping overly similar sequences, and generating an maximum likelihood (ML) tree using RAxML.

The user provides one protein sequence or multiple sequences and blast databases are used to sample sequences
from phyla selected by the user e.g., The following can be used for phyla found in the human gut i.e., *Actinobacteriota*, *Bacteroidota*, *Desulfobacterota*, *Firmicutes*, *Proteobacteria*, *Synergistota*, *Verrucomicrobiota*, *Fusobacteria*.

# Installation

The tool requires Python >= 3.11 and conda >= 4.12.0. The latest release can be installed directly from pip or this repository.

```
pip install autophylo
```

Or

Create a conda environment using the `environment.yml` file which installs all the dependencies (listed below).

```
conda env create -f environment.yml
```

# Usage

```
conda activate autophylo
```

# Install autophylo using tarball

Install the `autophylo` application within the created `autophylo` conda environment using a tarball.

```
python3 -m pip install /path/to/autophylo-1.0.0.tar.gz
```

# Dependencies

The following are required dependencies (listed below):

- NCBI BLAST 2.15.0
- BLAST databases (version 5)
- FastTree 2.1.11
- MUSCLE 5.1.0
- RAxML 8.2.13
- GBLOCKS 0.91b
- USEARCH 12.0_beta
- seqkit 2.8.2
- Trimal 1.5.0
- biopython 1.84
- joblib 1.4.2



# Download pre-formatted blast databases

(https://ftp.ncbi.nlm.nih.gov/blast/db/v5/README)

- The pre-formatted databases offer the following advantages:
    * Pre-formatting removes the need to run makeblastdb;
    * Species-level taxonomy ids are included for each database entry;
    * Databases are broken into smaller-sized volumes and are therefore easier
      to download;
    * Sequences in FASTA format can be generated from the pre-formatted databases
      by using the blastdbcmd utility;
    * A convenient script (update_blastdb.pl) is available in the blast+ package
      to download the pre-formatted databases.

# download nr

```
update_blastdb.pl --source ncbi --decompress --blastdb_version 5  --verbose 2 --num_threads 30 nr > log.nr 2>&1
```

# Update `config` file 

Obtain path to dependencies programs using the `$CONDA_PREFIX` variable.

After activating the `autophylo` run the following command to get the path and use it to update the `config` file

```
(autophylo) echo $CONDA_PREFIX
```

# Example `config` file

```
(autophylo) amos@Amogelangs-MacBook-Pro autophylo % echo $CONDA_PREFIX
/Users/amos/miniconda3/envs/autophylo
```

# Updated `config` file

NOTE: The databases can be place anywhere in the filesystem and in this example they are in `/Users/amos/datalake`.

```
[DEFAULT]
ServerAliveInterval = 45
Compression = yes
CompressionLevel = 9
ForwardX11 = yes

[ALIGNMENT]
MUSCLE=/Users/amos/miniconda3/envs/autophylo/bin/muscle
TRIMAL=/Users/amos/miniconda3/envs/autophylo/bin/trimal

[TREE]
FastTree=/Users/amos/miniconda3/envs/autophylo/bin/FastTree
RAxML_PTHREADS=/Users/amos/miniconda3/envs/autophylo/bin/raxmlHPC-PTHREADS-AVX
RAxML_HYBRID=/Users/amos/miniconda3/envs/autophylo/bin/raxmlHPC-HYBRID-AVX

[CLUSTERING]
GBLOCKS=/Users/amos/miniconda3/envs/autophylo/bin/Gblocks
USEARCH=/Users/amos/miniconda3/envs/autophylo/bin/usearch

[DATABASES]
BLASTDB=/Users/amos/datalake/BLASTDB/NR/nr
TAXIDS=/Users/amos/datalake/BLASTDB/NR/taxids
```