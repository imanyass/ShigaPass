#!/bin/bash


#module load blast+/2.12.0

usage() {
        echo
        echo "###### This script is used for determining Shigella serotypes  #####"
        echo "usage : ShigaPass_v1.2.sh -l <your_list> -o <output_directory> -p <databases_pathway>"
        echo
        echo "options :"
        echo "-l        List file contains the path of FASTA files (mandatory)"
        echo "-o        Output directory (mandatory)"
        echo "-p        Path to databases directory (mandatory)"
        echo "-u        Update the databases (Optional)"
	echo "-k	Keep intermediate files (Optional)"
        echo "-h        Display this help and exit"
        echo "Example: ShigaPass_v1.2.sh -l list_of_fasta.ls -o ShigaPass_Results -p ShigaPass/ShigaPass_DataBases -u"
        echo "Please note that the -u option should be used when running the script for the first time"
}

MKDB=0 
KEEP=0

while getopts "l:o:p:huk" option; do
        case "${option}" in
                l) LIST=${OPTARG};;
                o) OUTDIR=${OPTARG};;
                p) DBPATHWAY=${OPTARG};;
                u) MKDB=1;; #To update databases
                h) # display usage
                        usage
                        exit 0;;
		k) KEEP=1;; #To keep intermediate files
                \?) # incorrect option
                        echo "Error: Invalid option" 
                        usage
                        exit 1;;
        esac
done

echo $LIST
echo $OUTDIR

# Create the output directory if not present
if [ ! -d ${OUTDIR} ]
then
        mkdir ${OUTDIR}
fi



BLAST_OPT="-num_threads 4 -num_alignments 10000 -outfmt 6 -word_size 11 -dust no"
date="`date '+%d_%m_%Y__%H_%M_%S'`"
echo "Name;rfb;rfb_hits;MLST;fliC;CRISPR;ipaH;ipaH_hitss;Predicted_Serotype;Predicted_FlexSerotype" > ${OUTDIR}/ShigaPass_summary_${date}.csv


BLAST_awk () {
blastn -db ${DBPATHWAY}/$1 -query ${f} -out ${OUTDIR}/${NAMEDIR}/$2_blastout.txt -num_threads 4 -num_alignments 10000 -outfmt 6 -word_size 11 -dust no
awk -v ID="$3" -v COV="$4" -F "\t|_" '{if ($6>=ID && ($7/$5)*100>=COV) print $0}' ${OUTDIR}/${NAMEDIR}/$2_blastout.txt > ${OUTDIR}/${NAMEDIR}/$2_allrecords.txt
}

Hits_awk () {
awk -F "\t|_" 'OFS=";"{a[$2]++;} END{for(i in a) print i,a[i]}' ${OUTDIR}/${NAMEDIR}/$1_allrecords.txt |sort -k 2 -t ";" -nr -o ${OUTDIR}/${NAMEDIR}/$1_hits.txt
}



if [ $MKDB = 1 ]
then
	for database in ${DBPATHWAY}/*.fasta; do makeblastdb -dbtype nucl -in ${database}; done 
fi


for y in $(cat ${LIST})
do
        echo ${y}
        RFB=""
        MLST=""
        FLIC=""
        CRISPR=""
        gtrI=""
	gtrIC=""
        gtrII=""
        gtrX=""
        gtrIV=""
        gtrV=""
        oac=""
        oac1b=""
        optII=""
        FLEXSEROTYPE=""
        SEROTYPE=""
        ipaH=""
        ipah_hits=""
        confidence=""

    	NAMEDIR=$(basename $y .fasta)
        if [ ! -d ${OUTDIR}/${NAMEDIR} ]
        then
                mkdir ${OUTDIR}/${NAMEDIR}
	else 
                rm -r ${OUTDIR}/${NAMEDIR}/*
        fi

	sed 's/_/~/g' ${y} > ${OUTDIR}/${NAMEDIR}/${NAMEDIR}_parsed.fasta # Removing "_" from Fasta's name
	f="${OUTDIR}/${NAMEDIR}/${NAMEDIR}_parsed.fasta"

	FastaFile=${y##*/}
        FastaName=${FastaFile%%.*}

	echo "### ipaH checkpoint ###"
        BLAST_awk ipaH_150-mers.fasta ipaH 98 95
	Hits_awk ipaH
        ipaH=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/ipaH_hits.txt)
        if [[ ! -z "$ipaH" ]]
        then
                ipaH="ipaH+"
                ipaH_hits=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/ipaH_hits.txt  | head -n 1 | cut -f 2 -d ";") 
        else
                ipaH="ipaH-"
                ipaH_hits="0"
        fi
        echo $ipaH
        echo $ipaH_hits
	
	if [[  "$ipaH" == "ipaH-" ]]
	then
		SEROTYPE="Not Shigella/EIEC"
		RFB="ND"
		MLST="ND"
		FLIC="ND"
		CRISPR="ND"
	else

		echo "### Determining rfb ###" 
		BLAST_awk ABC-serotypes_150-mers.fasta rfb 98 95 
		Hits_awk rfb
		RFB=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/rfb_hits.txt | head -n 1 | cut -f 1 -d ";") 
		hit=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/rfb_hits.txt  | head -n 1 | cut -f 2 -d ";") 
		echo $RFB

		#Array declaration
		RFBs=("A2" "A3a" "A3b" "C10")
		FILEs=("AprovBEDP02-5104_150-mers.fasta" "A16_rfbU_150-mers.fasta" "A16_rfbU_150-mers.fasta" "taurine_SB6.fasta")
		NewRFBs=("AprovBEDP02-5104" "A16" "A16" "C6")

		for index in ${!RFBs[@]}
		do 
			if [[  "$RFB" == "${RFBs[$index]}" ]]
			then
				BLAST_awk ${FILEs[$index]} additionalrfb 98 95
				Hits_awk additionalrfb
				RFB=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/additionalrfb_hits.txt | head -n 1 | cut -f 1 -d ";")
				if [[ ! -z "$RFB" ]]
				then
					RFB=${NewRFBs[$index]}
					echo "rfb has changed to" ${NewRFBs[$index]}
					if [[ "$RFB" != "C6" ]]
                                        then
                                                hit=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/additionalrfb_hits.txt | head -n 1 | cut -f 2 -d ";")
                                        else
                                                hit=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/rfb_hits.txt  | head -n 1 | cut -f 2 -d ";")
                                        fi
				else
					RFB=${RFBs[$index]}
					echo "rfb has remained" ${RFBs[$index]}
				fi
				break
			fi
		done
	
		if [[  "$RFB" == "A3b" ]] #A3b is the common part found in both A3 and A16
		then
			RFB=$(sort -k 1 -t ";" ${OUTDIR}/${NAMEDIR}/rfb_hits.txt |head -n 1 | cut -f 1 -d ";")
			if [[ "$RFB" == "A3a" ]]
			then
				RFB="A3"
				echo "rfb has changed to A3; hits detected are unique for A3"
			else
				RFB="A3/A16"
				echo "hits detected are commun with A3 and Aprov97-10607"
			fi
		elif [[  "$RFB" == "A3a" ]]
		then
			RFB="A3"
                 	echo "rfb has changed to A3; hits detected are unique for A3"
		elif [[  "$RFB" == "C1" ]]  # search for galF gene which is normally present in SB1
		then
			BLAST_awk galF_SB1.fasta additionalrfb 98 95
			Hits_awk additionalrfb
			RFB=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/additionalrfb_hits.txt |head -n 1 | cut -f 1 -d ";")
			if [[ ! -z "$RFB" ]]
			then
				RFB="C1"
				echo "rfb has remained C1"
			else
				RFB="C20"
				echo "rfb has changed to C20"
			fi
	
		elif [[  "$RFB" == "B1-5" ]]
		then
			echo "### Determining phage and plasmid encoded O-antigen modification genes ###"
			BLAST_awk POAC-genes_150-mers.fasta POAC 98 95
			Hits_awk POAC
			sed 's/gtrX/32/g' ${OUTDIR}/${NAMEDIR}/POAC_hits.txt| sed 's/gtrII/4/g' |sed 's/gtrIC/2/g' | sed 's/gtrIV/8/g' | sed 's/gtrV/16/g' | sed 's/gtrI/1/g' | sed 's/oac1b/128/g' | sed 's/oac/64/g' |sed 's/optII/256/g' |\
			cut -f 1 -d ";" | awk '{total += $1} END{print "score="total + 0}' >${OUTDIR}/${NAMEDIR}/score.txt
			score=$(cut -f 2 -d "="  ${OUTDIR}/${NAMEDIR}/score.txt)
			echo $score
			SCOREs=("1" "129" "3" "131" "4" "36" "96" "64" "128" "8" "264" "72" "328" "16" "80" "48" "112" "32" "288" "0" "256")
			FLEX=("1a" "1b" "1c(7a)" "7b" "2a" "2b" "3a" "3b" "3b atypical (oac1b)" "4a" "4av" "4b" "4bv" "5a" "5a" "5b" "5b" "X" "Xv" "Y" "Yv")
			FLEXSEROTYPE="Unknown"
			for i in ${!SCOREs[@]}
			do
				if [[  "$score" == "${SCOREs[$i]}" ]]
				then
					FLEXSEROTYPE=${FLEX[$i]}
					break 
				fi
			done
		echo $FLEXSEROTYPE
		phages=$(sort -k 1 -t ";" ${OUTDIR}/${NAMEDIR}/POAC_hits.txt |cut -f 1 -d ";" |awk 'BEGIN { ORS = ";" } { print }' )
		echo "$FastaName;$phages;$FLEXSEROTYPE" | sed 's/;;/;/g' | tee -a ${OUTDIR}/ShigaPass_Flex_summary_${date}.csv

		elif [[  -z "$RFB"  || "$RFB" == "" ]] # if no rfb hit is detected, search for the presence of SS rfb
		then
			BLAST_awk D_150-mers.fasta additionalrfb 100 100
			Hits_awk additionalrfb
			hit=$(sort -k 2 -t ";" -n -r ${OUTDIR}/${NAMEDIR}/additionalrfb_hits.txt  | head -n 1 | cut -f 2 -d ";") 
			if [[ "$hit" -ge 30 ]]
			then
				RFB="D"
				echo "rfb has changed to D"
			else
				RFB="none"  # if there is no match at this point, we assume that there is no satisfying match
				echo $RFB
			fi
		else
			echo "NO additional blast is needed"
		fi

		echo "### Determining ST for the 7 genes  ###"
	
		for g in adk fumC gyrB icd mdh purA recA
		do 
			declare GENEST="ST_${g}"
			blastn -db ${DBPATHWAY}/${g}_len.fasta -query $f ${BLAST_OPT} -out ${OUTDIR}/${NAMEDIR}/${g}_blastout.txt
			GENEST=$(awk -F "\t|:" '{if ($5==100 && $6==$4) print $2}' ${OUTDIR}/${NAMEDIR}/${g}_blastout.txt | cut -f 2 -d "-") 
			#GENEST=$(blastn -db ${DBPATHWAY}/${g}_len.fasta -query $f ${BLAST_OPT} | awk -F '\t|:' '{if ($5==100 && $6==$4) print $2}' | cut -f 2 -d "-" )
			[[ -z $GENEST ]] && GENEST='ND'
			echo "${g}:$GENEST" | tee -a ${OUTDIR}/${NAMEDIR}/mlst_alleles.txt
		done

		echo "### Infering MLST ###"
		MOFILE="${OUTDIR}/${NAMEDIR}/mlst_alleles.txt"
		echo "$(grep adk $MOFILE | cut -f 2 -d ":" ) $(grep fumC $MOFILE |cut -f 2 -d ":" ) $(grep gyrB $MOFILE |cut -f 2 -d ":" ) $(grep icd $MOFILE |cut -f 2 -d ":" ) $(grep mdh $MOFILE |cut -f 2 -d : ) $(grep purA $MOFILE |cut -f 2 -d : ) $(grep recA $MOFILE |cut -f 2 -d : )"
		awk -v AWK_adk="$(grep adk $MOFILE | cut -f 2 -d ":" )" \
		-v AWK_fumC="$(grep fumC $MOFILE |cut -f 2 -d ":" )" \
		-v AWK_gyrB="$(grep gyrB $MOFILE |cut -f 2 -d ":" )" \
		-v AWK_icd="$(grep icd $MOFILE |cut -f 2 -d ":" )" \
		-v AWK_mdh="$(grep mdh $MOFILE |cut -f 2 -d ":" )" \
		-v AWK_purA="$(grep purA $MOFILE |cut -f 2 -d ":" )" \
		-v AWK_recA="$(grep recA $MOFILE |cut -f 2 -d ":" )" \
		'{if ($2==AWK_adk && $3==AWK_fumC && $4==AWK_gyrB && $5==AWK_icd && $6==AWK_mdh && $7==AWK_purA && $8==AWK_recA) print "ST"$1}' ${DBPATHWAY}/ST_profiles.txt > ${OUTDIR}/${NAMEDIR}/mlst_ST.txt # if a line of ST databank matches every ST, we print it
		MLST=$(cut -f 2 ${OUTDIR}/${NAMEDIR}/mlst_ST.txt | cut -f 2 -d ":")
		[[ ! -z "$MLST" ]] || MLST="none" # if there is no match at this point, we assume that there is definetely no match
		echo $MLST

		
		echo "### Determining fliC ###"
                BLAST_awk fliC_Shigella_v1.fasta flic 98 95
                sort -k 12 -n -r  ${OUTDIR}/${NAMEDIR}/flic_allrecords.txt|head -n 1 > ${OUTDIR}/${NAMEDIR}/flic_records.txt
                FLIC=$(cut -f 2 ${OUTDIR}/${NAMEDIR}/flic_records.txt | cut -f 2 -d "_" )
                [[ ! -z "$FLIC" ]] || awk -F '\t|_' '{if ($6>=98 && ($7/$5)*100>=45) print $0}' ${OUTDIR}/${NAMEDIR}/flic_blastout.txt |\
                (sort -k 3,3 -k 4,4 -k 12,12  -n -r | head -n 1) > ${OUTDIR}/${NAMEDIR}/flic_records.txt # if IS is suspected, we lower the %id 
                FLIC=$(cut -f 2 ${OUTDIR}/${NAMEDIR}/flic_records.txt | cut -f 2 -d "_" )
                [[ ! -z "$FLIC" ]] || FLIC="none" # if there is no match at this point, we assume that there is no satisfying match
                echo $FLIC

		echo "### Determining CRISPR-type ###"

		BLAST_awk CRISPR_spacers.fasta crispr 100 100
	        awk -F '\t|_' '{if ($12<$13) {print | "sort -nk7"} else if ($12>$13) {print | "sort -nrk7"}}' ${OUTDIR}/${NAMEDIR}/crispr_allrecords.txt |\
        	awk -F "_" '{a[$1]=a[$1]","$2} END {for (i in a) print a[i]}' | sed 's/,//' |sort -r |paste -sd ',' > ${OUTDIR}/${NAMEDIR}/crispr_records.txt
	        CRISPR=($(cat ${OUTDIR}/${NAMEDIR}/crispr_records.txt ))
        	[[ ! -z "$CRISPR" ]] || CRISPR="none" # if there is no match at this point, we assume that there is no satisfying match
	        echo $CRISPR
	

		echo "### Combining data into a serotype ###"
		
		SEROTYPE=$(awk -v AWK_mlst=${MLST} -v AWK_flic=${FLIC} -v AWK_crispr=${CRISPR} -v AWK_rfb=${RFB} -F ";"\
        	'{if ( $1==AWK_mlst && $2==AWK_flic && $3==AWK_crispr && $4==AWK_rfb ) print $5}' ${DBPATHWAY}/ShigaPass_Serotype_Profiles_v3.csv) # Comparing the obtained profile with our profiles dataset to infer the serotype
        	if [[ ! -z "$SEROTYPE" ]]
        	then
                        echo $SEROTYPE
                        echo "Profile matching 100%"
                else
                # if the serotype is not detected, compare the obtained profile by taking always the rfb and 2/3 of the rest of databases
                	SEROTYPE=$(awk -v AWK_mlst=${MLST} -v AWK_flic=${FLIC} -v AWK_crispr=${CRISPR} -v AWK_rfb=${RFB} -F ";"\
                	'{if ($1==AWK_mlst && $2==AWK_flic  && $4==AWK_rfb || $1==AWK_mlst && $3==AWK_crispr  && $4==AWK_rfb ||$2==AWK_flic && $3==AWK_crispr && $4==AWK_rfb ) print $5}' ${DBPATHWAY}/ShigaPass_Serotype_Profiles_v3.csv |cut -f 1 -d " " |head -n 1)
                	if [[ ! -z "$SEROTYPE" ]]
                	then
                                echo $SEROTYPE
                                echo "Profile matching 75%"
			else
			# if the serotype is not detected, search for a known Shigella MLST
				SEROTYPE=$(awk -v AWK_mlst=${MLST} -v AWK_flic=${FLIC} -v AWK_crispr=${CRISPR} -v AWK_rfb=${RFB} -F ";" \
				'{if ($1==AWK_mlst && $2==AWK_flic || $1==AWK_mlst && $3==AWK_crispr || $1==AWK_mlst ) print "unknown"}' ${DBPATHWAY}/ShigaPass_Serotype_Profiles_v3.csv |cut -f 1 -d " " |head -n 1)
				if [[ ! -z "$SEROTYPE" ]] && [[  "$rfb" != "none" ]] 
				then
					SEROTYPE="Shigella spp."
					echo $SEROTYPE
					echo "No profile matching with rfb, more probably contamination"
				elif [[ ! -z "$SEROTYPE" ]] && [[  "$rfb" == "none" ]]
				then
					SEROTYPE="Shigella spp."
					echo $SEROTYPE
					echo "No profile matching with rfb, more probabaly bad sequence quality"
				else 
					SEROTYPE="EIEC"
					echo $SEROTYPE
					FLEXSEROTYPE=""
					echo "No profile matching and ipaH+, More probabaly EIEC"
				fi
			fi	
		fi
	fi
	

echo "$FastaName;$RFB;$hit;$MLST;$FLIC;$CRISPR;$ipaH;$ipaH_hits;$SEROTYPE;$FLEXSEROTYPE" | tee -a ${OUTDIR}/ShigaPass_summary_${date}.csv
rm ${f}

if [ $KEEP = 0 ]  # Delete intermediate files
then
        rm -r ${OUTDIR}/${NAMEDIR}
fi


done 
 

