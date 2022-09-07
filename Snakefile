import os
import re
from glob import glob

configfile: "snakemake/config.yaml"
report: "code/report/workflow.rst"

db_version = "JUL_2022"

rule all:
    input:
        expand("outputs/TAXONOMY_DB_{version}.txt",version = "JUN_2022"),
        expand("Universal_Microbiomics_Alignment_Database/UNIPROT_INFO_{version}.txt.gz",version = "JUN_2022"),
        "BIOCYC_NF",
        ".get_function_names",
        ".biocyc_monomers_converted",
        "KEGG_GENES_RXN.txt",
        ".convert_kegg_genes",
        "BIOCYC_RXN_DB.txt",
        "BIOCYC_CPD_DB.txt",
        "KEGG_RXN_DB.txt",
        "KEGG_CPD_DB.txt",
        "PATHBANK_CPD_DB.txt",
        "HMDB_CPD_DB.txt",
        "PUBCHEM_CPD_DB.txt",
        "CHEBI_CPD_DB.txt",
        "RHEA_RXN_DB.txt",
        "MERGED_CPD_DB.txt",
        "MERGED_RXN_DB.txt"

rule make_rulegraph:
    output:
        "rulegraph.pdf",
        "rulegraph.png"
    shell:
        """
        snakemake all --rulegraph --dry-run | dot -Tpdf > rulegraph.pdf
        snakemake all --rulegraph --dry-run | dot -Tpng > rulegraph.png
        """

rule download_taxdump:
    output: "INPUTS/taxonomy/new_taxdump.tar.gz"
    shell:
        """
        mkdir -p INPUTS/taxonomy
        cd INPUTS/taxonomy
        wget -O new_taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
        """

rule build_tax_db:
    # Build the taxonomy database
    input:
        taxdump = "INPUTS/taxonomy/new_taxdump.tar.gz",
        img_genomes = "INPUTS/manual_inputs/tax_db/All_IMG_Genomes.txt",
        ictv = "INPUTS/manual_inputs/tax_db/ICTV.txt"
    output:
        tax_db_tsv = "outputs/TAXONOMY_DB_{version}.txt",
        tax_db_cyto = "outputs/TAXONOMY_DB_{version}.cyto"
    log: "snakemake/logs/tax_db/tax_db_{version}.log"
    benchmark: "snakemake/benchmarks/tax_db/tax_db_{version}.txt"
    resources: cpus = 1, mem_mb = 16000, time_min = 7200
    shell:
        """
        PROJ_ROOT=$PWD
        
        # Cleanup old database
        rm -rf Universal-Taxonomy-Database
        
        printf "\n\n######\nStarted at $(date)\n######\n\n" | tee $PROJ_ROOT/{log}
        start=`date +%s` | tee -a $PROJ_ROOT/{log}

        # Download GitHub Repository
        printf "\n\n######\nCloning GitHub Repository\n######\n\n" | tee -a $PROJ_ROOT/{log}
        git clone https://github.com/TealFurnholm/Universal-Taxonomy-Database.git | tee -a $PROJ_ROOT/{log}

        # Copy files obtained manually to working directory
        ln INPUTS/manual_inputs/tax_db/* Universal-Taxonomy-Database/

        # Enter working directory
        cd $PWD/Universal-Taxonomy-Database

        # Run the database construction scripts
        printf "\n\n######\nRunning perl script 1 of 3\n######\n\n" | tee -a $PROJ_ROOT/{log}
        perl Create_Taxonomy_Database_1of3.pl | tee -a $PROJ_ROOT/{log}

        printf "\n\n######\nRunning perl script 2 of 3\n######\n\n" | tee -a $PROJ_ROOT/{log}
        perl Create_Taxonomy_Database_2of3.pl | tee -a $PROJ_ROOT/{log}

        printf "\n\n######\nRunning perl script 3 of 3\n######\n\n" | tee -a $PROJ_ROOT/{log}
        perl Create_Taxonomy_Database_3of3.pl | tee -a $PROJ_ROOT/{log}

        cp TAXONOMY_DB_*.txt ../outputs/
        cp TAXONOMY_DB_*.cyto ../outputs/

        end=`date +%s` | tee -a $PROJ_ROOT/{log}
        printf "\n\n######\nEnded at $(date)\nRuntime was `expr $end - $start` seconds\n######\n\n" | tee -a $PROJ_ROOT/{log}
        """


rule get_biocyc:
    input: 
        script = "Get_BioCyc.pl"
    output:
        #directory("BIOCYC_NF")
        touch("get_biocyc.done")
    params:
        BIOCYC_FLATS_URL="http://brg-files.ai.sri.com/subscription/dist/flatfiles-52983746/index.html",
        BIOCYC_USER="biocyc-flatfiles",
        BIOCYC_PASS="data-20541"
    log: "snakemake/logs/biocyc/get_biocyc.log"
    benchmark: "snakemake/benchmarks/biocyc/get_biocyc.log"
    shell:
        """
        PROJ_ROOT=$PWD
        #mkdir -p {output}
        #cd {output}

        printf "\n\n######\nDownloading flat file index\n######\n\n" | tee -a $PROJ_ROOT/{log}
        wget -O biocyc_web.txt {params.BIOCYC_FLATS_URL} | tee -a $PROJ_ROOT/{log}

        printf "\n\n######\nParsing file file index\n######\n\n" | tee -a $PROJ_ROOT/{log}
        grep -oP "http.*\.tar\.gz" biocyc_web.txt > biocyc_files.txt

        printf "\n\n######\nGetting BioCyc with Get_BioCyc.pl\n######\n\n" | tee -a $PROJ_ROOT/{log}
        perl $PROJ_ROOT/{input.script} -pwd={params.BIOCYC_PASS} -usr={params.BIOCYC_USER} -in=biocyc_files.txt | tee -a $PROJ_ROOT/{log}
        """

rule get_function_names:
    input: 
        script = "Get_Function_Names.pl"
    output:
        #directory("INPUTS/function_names")
        "Function_Names.txt",
        touch(".get_function_names")
    log: "snakemake/logs/get_function_names/get_function_names.log"
    benchmark: "snakemake/benchmarks/get_function_names/get_function_names.log"
    shell:
        """
        PROJ_ROOT=$PWD
        #mkdir -p {output}
        #cd {output}

        perl $PROJ_ROOT/{input.script} | tee $PROJ_ROOT/{log}

        printf "\n\n######\nCompleted.\n######\n\n" | tee -a $PROJ_ROOT/{log}
        """
    
rule convert_BioCyc_monomers:
    input: 
        "get_biocyc.done",
        script = "Convert_BioCyc_Monomers.pl",
        #biocyc_dir = "BIOCYC_NF"
        #biocyc_dir = directory("/geomicro/data2/kiledal/old_UMRAD/Universal_Biological_Compounds_Database/BIOCYC_NF") #for testing before new download is done
    params:
        biocyc_dir = "BIOCYC_NF"
    output:
        touch(".biocyc_monomers_converted"),
        "BIOCYC_NF/BIOCYC_MONO_RXNS.txt"
    log: "snakemake/logs/convert_BioCyc_monomers/convert_BioCyc_monomers.log"
    benchmark: "snakemake/benchmarks/convert_BioCyc_monomers/convert_BioCyc_monomers.log"
    shell:
        """
        PROJ_ROOT=$PWD
        cd {params.biocyc_dir}

        perl $PROJ_ROOT/{input.script} | tee $PROJ_ROOT/{log}
        """

rule convert_kegg_genes:
    input: 
        script = "Convert_Kegg_Genes.pl"
    output: 
        #directory("INPUTS/kegg")
        "KEGG_GENES_RXN.txt",
        touch(".convert_kegg_genes")
    log: "snakemake/logs/convert_kegg_genes/convert_kegg_genes.log"
    benchmark: "snakemake/benchmarks/convert_kegg_genes/convert_kegg_genes.log"
    shell:
        """
        PROJ_ROOT=$PWD
        #mkdir -p {output}
        #cd {output}

        perl $PROJ_ROOT/{input.script} | tee $PROJ_ROOT/{log}
        """

rule create_biocyc_cpd_rxn_db:
    input: 
        script = "Create_CPD-RXN-DB_BIOCYC.pl",
    output: 
        "BIOCYC_RXN_DB.txt",
        "BIOCYC_CPD_DB.txt"
    conda: "snakemake/conda_yaml/main.yaml"
    log: "snakemake/logs/create_biocyc_cpd_rxn_db/create_biocyc_cpd_rxn_db.log"
    benchmark: "snakemake/benchmarks/create_biocyc_cpd_rxn_db/create_biocyc_cpd_rxn_db.log"
    shell:
        """
        PROJ_ROOT=$PWD

        perl $PROJ_ROOT/{input.script} | tee $PROJ_ROOT/{log}
        """

rule create_kegg_cpd_rxn_db:
    input: 
        script = "Create_CPD-RXN-DB_KEGG.pl",
        kegg_genes = "KEGG_GENES_RXN.txt"
    output:
        "KEGG_RXN_DB.txt",
        "KEGG_CPD_DB.txt"
    conda: "snakemake/conda_yaml/main.yaml"
    log: "snakemake/logs/create_kegg_cpd_rxn_db/create_kegg_cpd_rxn_db.log"
    benchmark: "snakemake/benchmarks/create_kegg_cpd_rxn_db/create_kegg_cpd_rxn_db.log"
    shell:
        """
        PROJ_ROOT=$PWD

        perl $PROJ_ROOT/{input.script} | tee $PROJ_ROOT/{log}
        """

rule create_other_cpd_rxn_db:
    input: 
        script = "Create_CPD-RXN-DB_OTHERDB.pl"
    output:
        "PATHBANK_CPD_DB.txt",
        "HMDB_CPD_DB.txt",
        "PUBCHEM_CPD_DB.txt"
    conda: "snakemake/conda_yaml/main.yaml"
    log: "snakemake/logs/create_other_cpd_rxn_db/create_other_cpd_rxn_db.log"
    benchmark: "snakemake/benchmarks/create_other_cpd_rxn_db/create_other_cpd_rxn_db.log"
    shell:
        """
        PROJ_ROOT=$PWD

        perl $PROJ_ROOT/{input.script} | tee $PROJ_ROOT/{log}
        """

rule create_rhea_cpd_rxn_db:
    input: 
        script = "Create_CPD-RXN-DB_RHEA.pl",
    output:
        "CHEBI_CPD_DB.txt",
        "RHEA_RXN_DB.txt"
    conda: "snakemake/conda_yaml/main.yaml"
    log: "snakemake/logs/create_rhea_cpd_rxn_db/create_rhea_cpd_rxn_db.log"
    benchmark: "snakemake/benchmarks/create_rhea_cpd_rxn_db/create_rhea_cpd_rxn_db.log"
    shell:
        """
        PROJ_ROOT=$PWD

        perl $PROJ_ROOT/{input.script} | tee $PROJ_ROOT/{log}
        """

rule merge_dbs:
    input: 
        "PATHBANK_CPD_DB.txt",
        "PUBCHEM_CPD_DB.txt",
        "HMDB_CPD_DB.txt",
        "KEGG_CPD_DB.txt",
        "BIOCYC_CPD_DB.txt",
        "CHEBI_CPD_DB.txt",
        script = "Merge_CPD_DB.pl"
    output:
        "MERGED_CPD_DB.txt",
        "MERGED_RXN_DB.txt"
    conda: "snakemake/conda_yaml/main.yaml"
    log: "snakemake/logs/merge_dbs/merge_dbs.log"
    benchmark: "snakemake/benchmarks/merge_dbs/merge_dbs.log"
    shell:
        """
        PROJ_ROOT=$PWD

        perl $PROJ_ROOT/{input.script} | tee $PROJ_ROOT/{log}

        cat *RXN_DB.txt > MERGED_RXN_DB.txt
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
    conda: "snakemake/conda_yaml/main.yaml"
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
        compounds_db = "MERGED_CPD_DB.txt",
        reactions_db = "MERGED_RXN_DB.txt"
        #function_names = rules.build_compounds_db.output.function_names
    output: 
        #uniprot_info = "Universal_Microbiomics_Alignment_Database/UNIPROT_INFO_{version}.txt.gz"
        uniprot_out = "Universal_Microbiomics_Alignment_Database/OUT_UNIPROTtest.txt"
    conda: "snakemake/conda_yaml/main.yaml"
    log: "Universal_Microbiomics_Alignment_Database/create_database.log"
    benchmark: "benchmark/build_alignment_db.txt"
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
        perl CurateUniRef.pl 1>>create_database.log 2>&1
        echo "\n\n***** CurateUniRef.pl is complete *****\n\n" 1>create_database.log 2>&1
        """

rule build_alignment_db2:
    input:
        script = "latest_Create_Alignment_DB.pl",
        inidm = "Universal_Microbiomics_Alignment_Database/idmapping.dat.gz",
        inup = "Universal_Microbiomics_Alignment_Database/uniprot-all.tab.gz",
        inpar = "Universal_Microbiomics_Alignment_Database/uniparc_all.xml.gz",
        inkegn = "KEGG_GENES_RXN.txt",
        inbmon = "BIOCYC_NF/BIOCYC_MONO_RXNS.txt",
        inrhrx = "RHEA_RXN_DB.txt",
        inkgrx = "KEGG_RXN_DB.txt",
        inbcrx = "BIOCYC_RXN_DB.txt",
        intcdb = "Universal_Microbiomics_Alignment_Database/UR100vsTCDB.m8",
        intrch = "Universal_Microbiomics_Alignment_Database/getSubstrates.py",
        infn = "Function_Names.txt",
        inurfa = "Universal_Microbiomics_Alignment_Database/uniref100.fasta.gz"
    output: 
        #uniprot_info = "Universal_Microbiomics_Alignment_Database/UNIPROT_INFO_{version}.txt.gz"
        #uniprot_out = "OUT_UNIPROTtest.txt"
        build_done = ".done_build_alignment_db"
    conda: "snakemake/conda_yaml/main.yaml"
    log: "create_alignment_database_20220803.log"
    benchmark: "benchmark/build_alignment_db.txt"
    resources: partition = "largemem", cpus = 1, mem_mb = 1200000
    shell:
        """
        ln -f {input.inidm} ./
        ln -f {input.inup} ./
        ln -f {input.inpar} ./
        ln -f {input.inbmon} ./
        ln -f {input.intcdb} ./
        ln -f {input.intrch} ./
        ln -f {input.inurfa} ./

        printf "*** Input files linked ***\n\n" 1>{log} 2>&1
       
        perl {input.script}  1>>{log} 2>&1
        printf "\n\n***** {input.script} is complete *****\n\n" 1>>{log} 2>&1
        """


rule download_RNAcentral:
    output:
        id_mapping = "id_mapping.tsv.gz",
        species_ids = "rnacentral_species_specific_ids.fasta.gz"
    shell:
        """
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