import os
import re
from glob import glob

configfile: "snakemake/config.yaml"
report: "code/report/workflow.rst"

db_version = "FEB_2022"

rule all:
    input:
        expand("Universal_Microbiomics_Alignment_Database/UNIPROT_INFO_{version}.txt.gz",version = db_version),
        expand("outputs/TAXONOMY_DB_{version}.txt", version = db_version),
        expand("outputs/all_compound_info_{version}.txt", version = db_version),
        expand("outputs/all_reaction_info_{version}.txt", version = db_version),
        expand("Universal_Microbiomics_Alignment_Database/UNIPROT_INFO_{version}.txt.gz", version = db_version),
        expand("Fix_RNACentral_Taxonomy/rnacentral_clean_{version}.fasta.gz", version = db_version),
        expand("outputs/all_reaction_info_{version}.txt", version = db_version)
        

rule download_taxdump:
    output: "INPUTS/new_taxdump.tar.gz"
    shell:
        """
        cd INPUTS
        wget -O new_taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
        """

rule build_tax_db:
    # Build the taxonomy database
    input:
        taxdump = "INPUTS/new_taxdump.tar.gz",
        img_genomes = "INPUTS/manual_inputs/tax_db/All_IMG_Genomes.txt",
        ictv = "INPUTS/manual_inputs/tax_db/ICTV.txt"
    output:
        tax_db_tsv = "outputs/TAXONOMY_DB_{version}.txt",
        tax_db_cyto = "outputs/TAXONOMY_DB_{version}.cyto"
    log: "logs/tax_db_{version}.log"
    benchmark: "benchmark/tax_db_{version}.txt"
    resources: cpus = 1, mem_mb = 16000, time_min = 7200
    shell:
        """
        # Cleanup old database
        rm -rf Universal-Taxonomy-Database
        
        printf "\n\n######\nStarted at $(date)\n######\n\n" | tee {log}
        start=`date +%s` | tee -a {log}

        # Download GitHub Repository
        printf "\n\n######\nCloning GitHub Repository\n######\n\n" | tee -a {log}
        git clone https://github.com/TealFurnholm/Universal-Taxonomy-Database.git

        # Copy files obtained manually to working directory
        cp manual_inputs/tax_db/* Universal-Taxonomy-Database/

        # Enter working directory
        cd $PWD/Universal-Taxonomy-Database

        # Run the database construction scripts
        printf "\n\n######\nRunning perl script 1 of 3\n######\n\n" | tee -a {log}
        perl Create_Taxonomy_Database_1of3.pl | tee -a {log}

        printf "\n\n######\nRunning perl script 2 of 3\n######\n\n" | tee -a {log}
        perl Create_Taxonomy_Database_2of3.pl | tee -a {log}

        printf "\n\n######\nRunning perl script 3 of 3\n######\n\n" | tee -a {log}
        perl Create_Taxonomy_Database_3of3.pl | tee -a {log}

        cp TAXONOMY_DB_*.txt ../outputs/
        cp TAXONOMY_DB_*.cyto ../outputs/

        end=`date +%s` | tee -a {log}
        printf "\n\n######\nEnded at $(date)\nRuntime was `expr $end - $start` seconds\n######\n\n" | tee -a {log}
        """


# rule git_clone_compounds_DB:
#     output: 
#         db_dir = directory("Universal_Biological_Compounds_Database")
#     shell:
#         """
#         printf "\n\n######\nCloning GitHub Repository\n######\n\n" | tee -a $PROJ_ROOT/{log}
#         git clone https://github.com/TealFurnholm/Universal_Biological_Compounds_Database.git | tee -a $PROJ_ROOT/{log}
#         """


rule prepare_compounds_db_inputs:
    input: "Universal_Biological_Compounds_Database"
    output:
        biocyc_file_list = "Universal_Biological_Compounds_Database/biocyc_files.txt",
        biocyc_dir = directory("Universal_Biological_Compounds_Database/BIOCYC_NF")
    log: "logs/prepare_compounds_db_inputs.log"
    benchmark: "benchmarks/prepare_compounds_db_inputs.txt"
    shell:
        """
        touch {log} #make log file
        PROJ_ROOT=$(pwd)

        # Set variables here
        BIOCYC_FLATS_URL="http://brg-files.ai.sri.com/subscription/dist/flatfiles-52983746/index.html"
        BIOCYC_USER="biocyc-flatfiles"
        BIOCYC_PASS="data-20541"

        # Download latest code for Universal Compounds database
        printf "\n\n######\nStarted at $(date)\n######\n\n" | tee $PROJ_ROOT/{log}
        start=`date +%s` | tee -a $PROJ_ROOT/{log}
        
        # Enter directory containing Universal Compounds Database
        cd $PWD/Universal_Biological_Compounds_Database

        printf "\n\n######\nGetting function names\n######\n\n" | tee -a $PROJ_ROOT/{log}
        perl Get_Function_Names.pl | tee -a $PROJ_ROOT/{log}

        printf "\n\n######\nDownloading flat file index\n######\n\n" | tee -a $PROJ_ROOT/{log}
        wget -O biocyc_web.txt $BIOCYC_FLATS_URL | tee -a $PROJ_ROOT/{log}

        printf "\n\n######\nParsing file file index\n######\n\n" | tee -a $PROJ_ROOT/{log}
        grep -oP "http.*\.tar\.gz" biocyc_web.txt > biocyc_files.txt

        printf "\n\n######\nGetting BioCyc with Get_BioCyc.pl\n######\n\n" | tee -a $PROJ_ROOT/{log}
        perl Get_BioCyc.pl -pwd=$BIOCYC_PASS -usr=$BIOCYC_USER -in=biocyc_files.txt | tee -a $PROJ_ROOT/{log}
        """

rule download_pathbank_metabolites:
    output: "INPUTS/pathbank_all_metabolites.csv.zip"
    shell:
        """
        cd INPUTS
        wget -N https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip
        """

rule download_hmdb_metabolites:
    output: "INPUTS/hmdb_metabolites.zip"
    shell:
        """
        cd INPUTS
        wget -N https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip
        """

rule download_pubchem_chebi:
    output: "INPUTS/PubChem_substance_text_chebi_summary.csv"
    shell:
        """
        cd INPUTS
        wget -O PubChem_substance_text_chebi_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22chebi%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_chebi%22}'
        """

rule download_pubchem_biocyc:
    output: "INPUTS/PubChem_substance_text_biocyc_summary.csv"
    shell:
        """
        cd INPUTS
        wget -O PubChem_substance_text_biocyc_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22biocyc%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_biocyc%22}'
        """

rule download_pubchem_hmdb:
    output: "INPUTS/PubChem_substance_text_hmdb_summary.csv"
    shell:
        """
        cd INPUTS
        wget -O PubChem_substance_text_hmdb_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22hmdb%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_hmdb%22}'
        """

rule download_pubchem_kegg:
    output: "INPUTS/PubChem_substance_text_kegg_summary.csv"
    shell:
        """
        cd INPUTS
        wget -O PubChem_substance_text_kegg_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22kegg%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_kegg%22}'
        """

rule download_chebi_ontology:
    output: "INPUTS/chebi.owl"
    shell:
        """
        cd INPUTS
        wget -N https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl
        """

rule build_compounds_db:
    # This script builds the compounds database
    # This script includes login information for BioCyc and should not be uploaded to a public directory while including such information
    input:
        rules.download_pathbank_metabolites.output,
        rules.download_hmdb_metabolites.output,
        rules.download_pubchem_chebi.output,
        rules.download_pubchem_biocyc.output,
        rules.download_pubchem_hmdb.output,
        rules.download_pubchem_kegg.output,
        rules.download_chebi_ontology.output,
        tax_db = rules.build_tax_db.output.tax_db_tsv
    output:
        compounds_db = "outputs/all_compound_info_{version}.txt",
        function_names = "outputs/Function_Names_{version}.txt"
    log: "logs/compounds_db_{version}.log"
    benchmark: "benchmark/compounds_db_{version}.txt"
    resources: cpus=1
    shell:
        """
        # Change to directory containing Universal Compounds Database
        cd $PWD/Universal_Biological_Compounds_Database

        perl Create_Compounds_Database.pl | tee {log}

        ln all_*.txt ../outputs/
        ln Function_Names_*.txt ../outputs/
        """

rule build_reactions_db:
    input: 
        tax_db = rules.build_tax_db.output.tax_db_tsv,
        compounds_db = rules.build_compounds_db.output.compounds_db
    output: 
        rxn_db = "outputs/all_reaction_info_{version}.txt"
    log: "logs/reactions_db_{version}.log"
    benchmark: "benchmark/reactions_db_{version}.txt"
    resources: cpus=1
    shell:
        """
        cd $PWD/Universal_Biological_Compounds_Database

        perl Create_Reactions_Database.pl | tee {log}
        cp all_*.txt ../outputs/
        """

rule get_alignemnt_db_git:
    output: directory("Universal_Microbiomics_Alignment_Database")
    shell:
        """
        git clone https://github.com/TealFurnholm/Universal_Microbiomics_Alignment_Database.git
        """

rule get_uniprot:
    output: "Universal_Microbiomics_Alignment_Database/uniprot-all.tab.gz"
    shell:
        """
        wget -O {output} 'https://www.uniprot.org/uniprot/?query=*&format=tab&force=true&columns=id,protein%20names,length,lineage-id,lineage(GENUS),lineage(SPECIES),organism,feature(SIGNAL),feature(TRANSMEMBRANE),database(TCDB),database(eggNOG),database(Pfam),database(TIGRFAMs),go-id,database(InterPro),ec,database(BioCyc),feature(DNA%20BINDING),feature(METAL%20BINDING),comment(SUBCELLULAR%20LOCATION),database(KEGG),rhea-id&compress=yes'
        """

rule get_uniref:
    output:
        uniref100 = "Universal_Microbiomics_Alignment_Database/uniref100.fasta.gz"
    shell:
        """
        cd Universal_Microbiomics_Alignment_Database

        wget -N https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
        """

rule get_uniparc:
    output:
        uniparc_all = "Universal_Microbiomics_Alignment_Database/uniparc_all.xml.gz"
    shell:
        """
        cd Universal_Microbiomics_Alignment_Database

        wget -N https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_all.xml.gz
        """

rule get_uniprot_mapping:
    output:
        uniprot_mapping = "Universal_Microbiomics_Alignment_Database/idmapping.dat.gz",
        uniprot_to_uniref = "Universal_Microbiomics_Alignment_Database/uniprot_to_uniref.txt.gz"
    shell:
        """
        cd Universal_Microbiomics_Alignment_Database

        wget -N https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
        zcat idmapping.dat.gz | grep -P "(UniRef100|UniRef90)" | gzip > uniprot_to_uniref.txt.gz
        """

rule get_tcdb:
    output:
        TCDB_get_script = "Universal_Microbiomics_Alignment_Database/getSubstrates.py",
        tcdb = "Universal_Microbiomics_Alignment_Database/tcdb.faa"
    shell:
        """
        cd Universal_Microbiomics_Alignment_Database

        wget --no-check-certificate -N https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py
        wget -O tcdb.faa http://www.tcdb.org/public/tcdb
        """

rule annotate_TCDB:
    input:
        tcdb = rules.get_tcdb.output.tcdb,
        uniref100 = rules.get_uniref.output.uniref100
    output:
        tcdb_diamond_db = "Universal_Microbiomics_Alignment_Database/tcdb.dmnd",
        UR100vsTCDB = "Universal_Microbiomics_Alignment_Database/UR100vsTCDB.m8"
    conda: "config/conda_yaml/main.yaml"
    resources: cpus = 16
    shell:
        """
        diamond makedb --in {input.tcdb} -d Universal_Microbiomics_Alignment_Database/tcdb
        diamond blastp -d {output.tcdb_diamond_db} -q {input.uniref100} -o {output.UR100vsTCDB} --query-cover 50 --top 0.5 --threads {resources.cpus} --strand both -f 6 qseqid qlen sseqid slen qstart qend sstart send evalue pident mismatch qcovhsp scovhsp
        """

rule build_alignment_db:
    input:
        rules.get_uniref.output.uniref100,
        rules.get_uniparc.output.uniparc_all,
        rules.get_uniprot_mapping.output.uniprot_mapping,
        rules.get_uniprot_mapping.output.uniprot_to_uniref,
        rules.get_tcdb.output.TCDB_get_script,
        rules.get_tcdb.output.tcdb,
        rules.build_tax_db.output.tax_db_cyto,
        uniprot = rules.get_uniprot.output,
        tcdb_diamond_db = rules.annotate_TCDB.output.tcdb_diamond_db,
        UR100vsTCDB = rules.annotate_TCDB.output.UR100vsTCDB,
        taxonomy_db = rules.build_tax_db.output.tax_db_tsv,
        compounds_db = rules.build_compounds_db.output.compounds_db,
        reactions_db = rules.build_reactions_db.output.rxn_db,
        function_names = rules.build_compounds_db.output.function_names
    output: 
        uniprot_info = "Universal_Microbiomics_Alignment_Database/UNIPROT_INFO_{version}.txt.gz"
    conda: "snakemake/conda_yaml/main.yaml"
    log: "Universal_Microbiomics_Alignment_Database/create_database_{version}.log"
    benchmark: "benchmark/build_alignment_db_{version}.txt"
    shell:
        """
        ln -f {input.taxonomy_db} Universal_Microbiomics_Alignment_Database/
        ln -f outputs/all*.txt Universal_Microbiomics_Alignment_Database/
        ln -f {input.function_names} Universal_Microbiomics_Alignment_Database/

        cd Universal_Microbiomics_Alignment_Database

        #perl CurateUniProt.pl 1>create_database_{wildcards.version}.log 2>&1

        echo "\n\n***** CurateUniProt.pl is complete *****\n\n" 1>create_database_{wildcards.version}.log 2>&1

        cp UNIPROT_INFO_{wildcards.version}int.txt UNIPROT_INFO_{wildcards.version}.txt
        gzip UNIPROT_INFO_{wildcards.version}.txt

        perl CurateUniRef.pl 1>>create_database_{wildcards.version}.log 2>&1

        echo "\n\n***** CurateUniRef.pl is complete *****\n\n" 1>create_database_{wildcards.version}.log 2>&1
        """


rule download_RNAcentral:
    output:
        id_mapping = "Fix_RNACentral_Taxonomy/id_mapping.tsv.gz",
        species_ids = "Fix_RNACentral_Taxonomy/rnacentral_species_specific_ids.fasta.gz"
    shell:
        """
        cd Fix_RNACentral_Taxonomy
        wget http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz
        wget http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_species_specific_ids.fasta.gz
        """

rule build_ncRNA_db:
    input:
        script = "Fix_RNACentral_Taxonomy/FIX_RNACENTAL_HEADERS.pl",
        rnacentral_mapping = rules.download_RNAcentral.output.id_mapping,
        rnacentral_speices_ids = rules.download_RNAcentral.output.species_ids,
        tax_db = rules.build_tax_db.output.tax_db_tsv
    output: "Fix_RNACentral_Taxonomy/rnacentral_clean_{version}.fasta.gz"
    conda: "snakemake/conda_yaml/main.yaml"
    log:  "Fix_RNACentral_Taxonomy/build_log_{version}.log"
    benchmark:  "benchmark/build_ncRNA_db_{version}.txt"
    shell:
        """
        ln -f {input.tax_db} Fix_RNACentral_Taxonomy/
        cd Fix_RNACentral_Taxonomy

        perl $(basename {input.script}) #1> $(basename {log}) 2>&1

        ln -f rnacentral_clean.fasta.gz rnacentral_clean_{wildcards.version}.fasta.gz
        """
