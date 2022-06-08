# ShigaPass 

ShigaPass is a new *in silico* tool used to predict *Shigella* serotypes and differentiate between *Shigella*, EIEC (Enteroinvasive *E. coli*), and non *Shigella*/EIEC using assembled whole genomes.

## Dependecies
Blast+ version 2.12.0
## Installation
**1.** Clone this repository with the following command line:
```sh
git clone https://github.com/imanyass/ShigaPass.git
```
**2.** B. Give the execute permission to the file ShigaPass_v1.0.sh:
```
chmod +x ShigaPass_v1.0.sh
```
**3.** Execute ShigaPass  with the following command line model:
```sh
./ShigaPass_v1.0.sh  [options]
```
## Usage 
Run ShigaPass with -h option to read the following documentation:
````
usage : bash ShigaPass_v1.0.sh -l <your_list> -o <output_directory> -p <databases_pathway>

-l        List file contains the path of FASTA files (mandatory)
-o        Output directory (mandatory)
-p        Path to databases directory (mandatory)
-u        Update the databases (Optional)
-k        Keep intermediate files (Optional)
-h        Display this help and exit

Example: bash ShigaPass_v1.0.sh -l list_of_fasta.ls -o ShigaPass_Results -p ShigaPass/DATABASES -u
Please note that the -u option should be used when running the script for the first time"
````
