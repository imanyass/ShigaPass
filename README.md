# ShigaPass 

ShigaPass is a new *in silico* tool used to predict *Shigella* serotypes and differentiate between *Shigella*, EIEC (Enteroinvasive *E. coli*), and non *Shigella*/EIEC using assembled whole genomes.

## Dependecies
Blast+ version 2.12.0
## Installation
**1.** Clone this repository with the following command line:
```sh
git clone https://github.com/imanyass/ShigaPass.git
```
**2.** Give the execute permission to the file ShigaPass.sh:
```
chmod +x ShigaPass.sh
```
**3.** Execute ShigaPass  with the following command line model:
```sh
./ShigaPass.sh  [options]
```
## Usage 
Run ShigaPass with -h option to read the following documentation:
````
usage : ShigaPass_v1.2.sh -l <your_list> -o <output_directory> -p <databases_pathway>

-l        List file contains the path of FASTA files (mandatory)
-o        Output directory (mandatory)
-p        Path to databases directory (mandatory)
-u        Update the databases (Optional)
-k        Keep intermediate files (Optional)
-h        Display this help and exit

Example: ShigaPass.sh -l list_of_fasta.ls -o ShigaPass_Results -p ShigaPass/ShigaPass_DataBases -u
Please note that the -u option should be used when running the script for the first time
````

# Output Files
* In the output directory, two files will be written:
  1. ShigaPass_summary_date.csv: tab-delimited file with one row per genome inclinding the sample name, type of *rfb*, number of *rfb* hits, MLST ST, type of *fliC*, CRISPR spacers, the presence of *ipaH*, number of *ipaH* hits, the predicted serotype and *S. flexneri* subserotype.
  2. ShigaPass_Flex_summary_date.csv: tab-delimited file detailing the phage and plasmid-encoded O-antigen modification (POAC) genes detected for the predicted *S. flexneri* genomes
 * In case -k option is used, a directory will be created for every assembled genome and will contain several files
 
| Extension       | Description |
| :------------- |:-------------| 
|blastout.txt  |Blast results in tabular format|
|allrecords.txt     |Blast hits that passed the selected thresholds| 
|records.txt     |The best blast hit that passed the selected thresholds| 
|hits.txt     |Name and number of hits that passed the selected thresholds (only for k-mers databases: *rfb*, *ipaH* and POAC genes)| 
