use strict;
use warnings;
use Sort::Naturally 'nsort';


#die;
##############################################################
#OPEN FILES - CHECK
##############################################################
my $inidm  = 'idmapping.dat.gz';		-e $inidm  || die "1 no $inidm  - please place into this folder and restart\n";
my $inup  	= 'uniprot-all.tab.gz';		-e $inup   || die "2 no $inup   - please place into this folder and restart\n";
my $inpar	= 'uniparc_all.csv.gz';		-e $inpar  || die "3 no $inpar   - please place into this folder and restart\n";
my $urfa	= 'uniref100.fasta.gz';         -e $urfa   || die "4 no $urfa   - please place into this folder and restart\n";
my $inkegn	= 'KEGG_GENES_RXN.txt';		-e $inkegn || die "5 no $inkegn - please place into this folder and restart\n";
my $inbmon	= 'BIOCYC_MONO_RXNS.txt';	-e $inbmon || die "6 no $inbmon - please place into this folder and restart\n";
my $inrhrx	= 'RHEA_RXN_DB.txt';		-e $inrhrx || die "7 no $inrhrx - please place into this folder and restart\n";
my $inkgrx = 'KEGG_RXN_DB.txt';		-e $inkgrx || die "8 no $inkgrx - please place into this folder and restart\n";
my $inbcrx	= 'BIOCYC_RXN_DB.txt';		-e $inbcrx || die "9 no $inbcrx - please place into this folder and restart\n";
my $intcdb	= 'UR100vsTCDB.m8';		-e $intcdb || die "10 no $intcdb - please place into this folder and restart\n";
my $intrch	= 'getSubstrates.py';		-e $intrch || die "11 no $intrch - please place into this folder and restart\n";
my $infn	= 'Function_Names.txt';		-e $infn   || die "12 no $infn - please place into this folder and restart\n";


#output
open(OUTPINT, ">", "OUT_UNIPROTtest.txt")||die;
open(OUTREF, ">", "OUT_UNIREFtest.txt")||die;


##############################################################
#			BEGIN LOADING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



my (%UR100_INFO, %UR100_LEN, %UR100_NAME);


##############################################################
##### 		Fix/Get UniRef #s and add basic info	######
##############################################################
#GET PROT LENGTH AND ODD AA DISTs
my $on = 0; my $time = localtime;
print "GET PROT LEN AND ODD AA STATS $time\n";
open(INURFA, "unpigz -c $urfa |") || die "failed to open $urfa";
my $ur100; my @seq = (); my ($len, $tid, $name);
while(<INURFA>){
    $_ = uc;
    if (/^>/) {
        # finish current record
        $len = length(join("",@seq));
        if($len < 10) {
            # first line or seq too short, don't store any of this
        } else {
            $UR100_INFO{$ur100}{7}{$tid}++;
            $UR100_LEN{$ur100}=$len;
            $UR100_NAME{$ur100}=$name;
            if(progress($on)){$time=localtime; print "on $on time $time ur100 $ur100 len $len tid $tid name $UR100_NAME{$ur100}\n";} $on++;
            #if($on>10000000){last;}#!!!!
        }

        # start new record
        m/(UNIREF100\_\S+)\s+(.*?)\s+N\=\d.*TAXID\=(\d+|N\/A)/ or die "failed parsing fasta header: $_";
        $ur100=$1;
        # $name = CleanNames($2);
        $name = $2;
        $tid = $3;
        @seq = ();
    } else {
        # collect sequence
        chomp;
        push(@seq, $_);
    }
}
close INURFA;
# finish final record (if any)
$len = length(join("",@seq));
if($ur100 and $len >= 10) {
    $UR100_INFO{$ur100}{7}{$tid}++;
    $UR100_LEN{$ur100}=$len;
    $UR100_NAME{$ur100}=$name;
    $on++;
}
my $num_info_recs = keys %UR100_INFO;
my $num_len_recs = keys %UR100_LEN;
my $num_name_recs = keys %UR100_NAME;
print "DONE $on lines #INFO=$num_info_recs #LEN=$num_len_recs #NAME=$num_name_recs\n\n";

#LOAD IDMAP
my (%KGEN_KRXN, %UR100_UR90, %UPID_UR100, %UR90_INFO);
$on=0; $time=localtime; my $seq_missing=0; my $print_progress = 0; my $no_ur100 = 0;
print "INPUT UNIREF TO UNIPROT MAPPING $time\n";
my $upid=''; $ur100=''; my $ur90=''; my @KEGG=(); $tid='';
my ($totupid, $totur1);
open(INMAP,  "unpigz -c $inidm |") || die "failed opening $inidm";
while(<INMAP>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;

	(my $prot, my $type, my $id)=split("\t",$_); 
	if($upid eq ''){$upid = $prot; $ur100=''; $ur90=''; @KEGG=(); $tid='';}
	if($prot eq $upid){ 
		   if($type eq "UNIREF100"){ 	$ur100=$id; }
        	elsif($type eq "UNIREF90"){  	$ur90=$id; }
		elsif($type eq "NCBI_TAXID"){	$tid=$id; }
		elsif($type eq "KEGG"){foreach my $k (keys %{$KGEN_KRXN{$id}}){ push(@KEGG,$k);}}
		else{}
	}
	else{	#reached the next protein, set current data and reset
		$on++;
                if(progress($on)){ $print_progress=1};
		if($ur100 !~/UNIREF100/){$upid = ''; $no_ur100++; next;}
                if (!exists($UR100_LEN{$ur100})) {$seq_missing++;}
		if($ur90 =~ /UNIREF90/){$UR100_UR90{$ur100}=$ur90;}
		$UPID_UR100{$upid}=$ur100;
		foreach my $k (@KEGG){ if($k=~/^R\d+$/){ $UR100_INFO{$ur100}{17}{$k}++; $UR90_INFO{$ur90}{17}{$k}++; }}
		if($tid=~/^\d+$/){$UR100_INFO{$ur100}{7}{$tid}++;}
                if($print_progress){
                    $print_progress = 0;
                    $time = localtime;
                    $totupid = keys %UPID_UR100;
                    $totur1  = keys %UR100_UR90;
                    $time = localtime;
                    print "on $on time $time type $type upid $upid prot $prot ur100 $ur100 ur90 $ur90 id $id upid_ur100 $totupid UR100_UR90 $totur1\n";
                }
		$upid='';
	}
	#if($on>1000000){last;}#!!!!
}
close INMAP;
$num_info_recs = keys %UR100_INFO;
my $num_upid2ur100 = keys %UPID_UR100;
print "DONE $on UPIDs; but no sequence: $seq_missing, entries w/o UR100: $no_ur100; #INFO=$num_info_recs #UPID2UR100=$num_upid2ur100\n\n";



##############################################################
##############################################################
#GET TCDB SUBSTRATES
# - right now used to pick best TCDB
# - previous used as left/right cpds - compile from RXN_DB chebi
$time=localtime; $on=0;
print "INPUT TRANSPORTER SUBSTRATES $time\n";
open(INTRCH, $intrch) || die "failed to open $intrch";
my %TCDB_CPDS;
my (@CHEBS, $chebi,);
while(<INTRCH>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\n\r]+//;
        (my $tcdb, my $chebs)=split("\t", $_);
        if($tcdb!~/TCDB/){$tcdb="TCDB:".$tcdb;}
        @CHEBS=();
        @CHEBS = ( $chebs =~ /(CHEBI.\d+)/g );
        $chebi=join(";",@CHEBS);
        $TCDB_CPDS{$tcdb}=$chebi;
	if(progress($on)){$time=localtime; print "on $on $time tcdb $tcdb chebi $chebi\n";} $on++;
}
close INTRCH;
print "DONE $on lines\n\n";

#INPUT TCDB ALIGNMENTS
$time=localtime; $on=0;
print "INPUT UNIREF100 vs TCDB MATCHES $time\n";
open(INTCDB, $intcdb) || die "failed to open $intcdb";
my (@stuff, $tcdb, $pid, $cov, $sco, %UR100_TCDB_SCO, %HASCPD);
my %UR100_TCDB;
while(<INTCDB>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_,-1);
        $ur100 = $stuff[0];
        if($stuff[2]=~/([^\|]+$)/){$tcdb = "TCDB:".$1;}
        else{$tcdb='';}
        $pid  = $stuff[9];	if($pid < 80){next;}
        $cov  = $stuff[11];	if($cov < 80){next;}
        $sco  = $pid*$cov;

        #get best hit, weight to tcdbID w/cpds
        if(!exists($UR100_TCDB{$ur100})){ 	#create ur100 - tcdb
		$UR100_TCDB{$ur100}=$tcdb; 	
		$UR100_TCDB_SCO{$ur100}=$sco; 	
		if($TCDB_CPDS{$tcdb}=~/CHEBI/){$HASCPD{$ur100}=1;}
		else{$HASCPD{$ur100}=0;}
	}
	else{ 	#already have a ur100-tcdb
		if($sco > $UR100_TCDB_SCO{$ur100}){ 								#is a better scoring tcdb match
			if($TCDB_CPDS{$tcdb}=~/CHEBI/){ 							#tcdb has cpd  
				$HASCPD{$ur100}=1; $UR100_TCDB{$ur100}=$tcdb; $UR100_TCDB_SCO{$ur100}=$sco;} 	#current has better score and a chebi 
			elsif($TCDB_CPDS{$tcdb}!~/CHEBI/ && $HASCPD{$ur100}==0){
				$UR100_TCDB{$ur100}=$tcdb; $UR100_TCDB_SCO{$ur100}=$sco;} 			#no chebi but better score
			else{}
	}	}

	if(progress($on)){$time=localtime; print "on $on time $time tcdb $tcdb ur100 $ur100 sco $sco hascpd $HASCPD{$ur100}\n";} $on++;
        #  if($on>10000){last;}
}
close INTCDB;
print "DONE $on lines\n\n";
undef(%HASCPD);
undef(%UR100_TCDB_SCO);
undef(%TCDB_CPDS);
##############################################################
##############################################################




##############################################################
##############################################################
#LOAD KEGG GENES
$time=localtime; $on=0;
print "INPUT KEGG GENES $time\n";
open(INKEGN, $inkegn) || die "failed to open $inkegn";
while(<INKEGN>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $kgene, my $rxn)=split("\t",$_);
	$KGEN_KRXN{$kgene}{$rxn}++;
        if(progress($on)) {$time=localtime;
                print "on $on time $time kgene $kgene rxn $rxn cnt $KGEN_KRXN{$kgene}{$rxn}\n";
        } $on++;
	#if($on>100000){last;}#!!!!
}
close INKEGN;
print "DONE $on lines\n\n";

#LOAD BIOCYC MONOMERS
$time=localtime; $on=0;
print "INPUT MONOMERS $time\n";
open(INBMON, $inbmon) || die "failed to open $inbmon";
my %MONO_BRXN;
while(<INBMON>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $mono, my $rxn)=split("\t",$_);
	$MONO_BRXN{$mono}{$rxn}++;
        if(progress($on)) {$time=localtime;
                print "on $on time $time mono $mono rxn $rxn cnt $MONO_BRXN{$mono}{$rxn}\n";
        } $on++;
	#if($on>100000){last;}#!!!!
}
close INBMON;
print "DONE $on lines\n\n";

#GET FUNCTION NAMES
$time=localtime; $on=0;
print "INPUT $infn $time\n";
open(INFN, $infn) || die "failed to open $infn";
my %FUNC_NAMES;
while(<INFN>){
	if($_ !~ /\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	(my $id, my $name)=split("\t",$_,-1);
	$id=~s/\s//g;
	$FUNC_NAMES{$id}=$name;
        if(progress($on)) {$time=localtime;
		print "on $on time $time id $id name $name\n";
        } $on++;
}
close(INFN);
print "DONE $on lines\n\n";
##############################################################
##############################################################


#INPUT UNIPARC FUNCTIONS
$on=0; my $picked=0; my $skipped=0; $time=localtime; $print_progress=0;
print "INPUT UNIPARC $time\n";
open(INPARTAB,  "unpigz -c $inpar |") || die "failed to open $inpar";
my (@PFAM, @TIGR, @IPR);
my ($pfam, $tigr, $ipr);
while(<INPARTAB>){
        chomp;
        ($upid, $pfam, $tigr, $ipr) = split "\t";

        $on ++;
	$ur100="UNIREF100_".$upid;
	if (! exists($UR100_LEN{$ur100})) {$skipped++; next;} #no ur100, skip

        @PFAM = split(';', $pfam);
        @TIGR = split(';', $tigr);
        @IPR = split(';', $ipr);

	$ur90=$UR100_UR90{$ur100};
	if($ur90!~/UNIREF90/) {$ur90="UNIREF90_".$upid;}
	if(!exists($UR90_INFO{$ur90})) {$ur90='';}

        #output functions
        foreach my $id (@PFAM){ $UR100_INFO{$ur100}{12}{$id}++; if($ur90=~/\w/){$UR90_INFO{$ur90}{12}{$id}++;}}
        foreach my $id (@TIGR){ $UR100_INFO{$ur100}{13}{$id}++; if($ur90=~/\w/){$UR90_INFO{$ur90}{13}{$id}++;}}
        foreach my $id (@IPR ){ $UR100_INFO{$ur100}{15}{$id}++; if($ur90=~/\w/){$UR90_INFO{$ur90}{15}{$id}++;}}

        $picked++;
        if (progress($picked)) {
            $time = localtime;
            print "park $on time $time picked=$picked upid $upid ur100 $ur100 ur90 $ur90 pfam @PFAM tigr @TIGR ipr @IPR\n";
        }
}
close INPARTAB;
$num_info_recs = keys %UR100_INFO;
print "DONE $on records total; picked: $picked; skipped (no sequence): $skipped; #NAME=$num_info_recs\n\n";
##############################################################
##############################################################




##############################################################
##############################################################
#LOAD UPID AND EC -> RXNS
#RHEA
$on=0; $time=localtime;
print "INPUT RHEA RXN $time\n";
open(INRHRX, $inrhrx) || die "failed to open $inrhrx";
my (%EC_RRXN, %UPID_RRXN);
my ($rxn, $ec);
while(<INRHRX>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff	=split("\t",$_);
	$rxn	=$stuff[0];
	if($rxn !~ /^\d+$/){next;}
	$upid	=$stuff[24];
	$ec	=$stuff[25];
	if($ec=~/[\d\.]+/){ $EC_RRXN{$ec}{$rxn}++; }
	if($upid=~/\w/){    $UPID_RRXN{$upid}{$rxn}++; }
	if(progress($on)) { print "rrxn $rxn ec $ec upid $upid\n"; } $on++;
}
close INRHRX;
print "DONE $on lines\n\n";
#KEGG
$on=0; $time=localtime;
print "INPUT KEGG RXN $time\n";
open(INKGRX, $inkgrx) || die "failed to open $inkgrx";
my (%EC_KRXN, %UPID_KRXN);
while(<INKGRX>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff	=split("\t",$_);
	$rxn	=$stuff[0];
	$upid	=$stuff[24];
	$ec	=$stuff[25];
	if($ec=~/[\d\.]+/){ $EC_KRXN{$ec}{$rxn}++; }
	if($upid=~/\w/){    $UPID_KRXN{$upid}{$rxn}++; }
	if(progress($on)) { print "krxn $rxn ec $ec upid $upid\n"; } $on++;
}
close INKGRX;
print "DONE $on lines\n\n";
#BIOCYC
$on=0; $time=localtime;
print "INPUT BIOCYC RXN $time\n";
open(INBCRX, $inbcrx) || die "failed to open $inbcrx";
my (%EC_BRXN, %UPID_BRXN);
while(<INBCRX>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff	=split("\t",$_);
	$rxn	=$stuff[0];
	$upid	=$stuff[24];
	$ec	=$stuff[25];
	if($ec=~/[\d\.]+/){ $EC_BRXN{$ec}{$rxn}++; }
	if($upid=~/\w/){    $UPID_BRXN{$upid}{$rxn}++; }
	if(progress($on)) { print "bioc $rxn ec $ec upid $upid\n"; } $on++;
}
close INBCRX;
print "DONE $on lines\n\n";
##############################################################
##############################################################



##############################################################
##### 		INPUT / OUTPUT UNIPROT DOWNLOAD		######
##############################################################
$on=0; $time=localtime; $skipped = 0;
print "INPUT UNIPROT $time\n";
print OUTPINT "UP-ID\tUR100\tUR90\tName\tLength\t";
print OUTPINT "SigPep\tTMS\tDNA\tTaxonId\tMetal\tLoc\t";
print OUTPINT "TCDB\tCOG\tPfam\tTigr\tGene_Ont\tInterPro\tECs\t";
print OUTPINT "kegg\trhea\tbiocyc\n";
open(INUP, "unpigz -c $inup |") || die "failed to open $inup";
my ($plen, @NAMES, @GN, @KNS, $n, $sig, $tms, @TMS, $tcp, @TCDBS, %seen, $cog, @COG);
my ($pfa, $tig, $gos, @GOS, $ecs, @ECS, @MON, $dna, $cln, $met, @METS, @GME, $loc);
my (@KEGS, @RHEA, %BRXN, %RRXN, %KRXN, @BIO, @RHE, @KEG, $kegg, $rhea, $bioc, @FIN);
my (@IDS, @GIDS, $out);
while(<INUP>){
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_);

	#screen bad entries
	if($stuff[0] =~ /^(ENTRY|PID)/i){	$skipped++; next;}	#skip column headers
	if($stuff[0] !~ /\w+/){			$skipped++; next;}	#no PID
	if($stuff[2] < 10 ){			$skipped++; next;}	#prot <10 a.a.
	if($stuff[3] !~ /^\d+$/){		$skipped++; next;}	#no taxid

	#SET BASIC VARIABLES
	$upid   = $stuff[0];
	if(!exists($UPID_UR100{$upid})){	$skipped++; next;}
	$ur100	=$UPID_UR100{$upid};
        # $ur90 	=$UPID_UR90{$upid};  # FIXME: what is UPID_UR90 for?
	$name	=$stuff[1];
	$plen	=$stuff[2];
	if($plen!~/\d/){$plen=$UR100_LEN{$ur100};}
	$tid	=$stuff[3];

	# FIX NAMES
	@NAMES=split('\(',$stuff[1]);
	@GN=(); @KNS=();
	#clean and get names 4+ characters
	foreach my $i (@NAMES){ $i = $i; if($i!~/[A-Z]{4,}/){next;} push(@GN,$i);}  # TODO: what about the name cleaning?
	while($GN[0]=~/\w/){ #get rid of duplicated names after clean
		$n=shift(@GN);
		if(grep $n, @KNS || grep $n, @GN){}
		else{ push(@KNS,$n); }
        }
	@KNS=nsort(@KNS); 
	$name=join(";",@KNS);
	if($UR100_NAME{$ur100}=~/[A-Z]{4,}/ && $name !~ /[A-Z]{4,}/){ $name=$UR100_NAME{$ur100};}

	#COMPILE/CLEAN UNIPROT ANNOTATIONS
	$sig=''; 		if($stuff[7]=~/(\d+)\.\.(\d+)/){						$sig = "SIGNAL:".$1."..".$2;}

	$tms=''; @TMS=(); 	if($stuff[8]=~/TRANSMEM|TMS/){@TMS = ($stuff[8]=~/(\d+\.\.\d+)/g); 		$tms  = join(";",@TMS); $tms="TMS:".$tms;}

	$tcp=''; @TCDBS=();	@TCDBS = ($stuff[9]=~/([A-Z\d\.]+)/g); #UNIPROT ANNOTATED TCDBS
		if($UR100_TCDB{$ur100}=~/[A-Z\d\.]+/){push(@TCDBS,$UR100_TCDB{$ur100});} #ALIGNMENT ANNOTATED TCDBS
		for my $i (0..$#TCDBS){ if( $TCDBS[$i] !~ /^TCDB/ && $TCDBS[$i]=~/\w/){ $TCDBS[$i]="TCDB:".$TCDBS[$i]; }} 
		%seen=(); @TCDBS = grep{ !$seen{$_}++ } @TCDBS;							$tcp = join(";", @TCDBS);

	$cog=''; @COG=();	if($stuff[10]=~/[CK]OG\d+/){@COG  = ($stuff[10]=~/\b([CK]OG\d{4})\b/g);		$cog = join(";",@COG);}

	$pfa=''; @PFAM=();	if($stuff[11]=~/(PF\d+)/){  @PFAM = ($stuff[11]=~/(PF\d+)/g); 			$pfa = join(";",@PFAM);}

	$tig=''; @TIGR=();	if($stuff[12]=~/(TIGR\d+)/){@TIGR = ($stuff[12]=~/(TIGR\d+)/g); 		$tig = join(";",@TIGR);}

	$gos=''; @GOS=();	if($stuff[13]=~/GO.\d+/){   @GOS  = ($stuff[13]=~/(GO.\d+)/g); 			$gos = join(";",@GOS);}

	$ipr=''; @IPR=();	if($stuff[14]=~/(IPR\d+)/){ @IPR  = ($stuff[14]=~/(IPR\d+)/g);			$ipr = join(";",@IPR);}

	$ecs=''; @ECS=();	if($stuff[15]=~/[\d\.]+/){  @ECS  = ($stuff[15]=~/(\d+\.\d+\.\d+\.\d+)/g);	$ecs = join(";",@ECS);}

		 @MON=();  	if($stuff[16]=~/MONOMER/){  @MON  = ($stuff[16] =~ /([\-\w]*MONOMER[\-\w]*)/g); } 

	$dna=''; 		if($stuff[17]=~/(\d+\.\.\d+).*?NOTE\=\"([^\"]+)/){ $cln=$2; 	$dna = "DNA:".$1."|".$cln;}

	$met=''; @METS=();	if($stuff[18]=~/NOTE/){	    @METS = ($stuff[18]=~/NOTE\W+([\w\s]+\w).*?[\"\;]+/g );
					@GME=(); for my $i (0..$#METS){ $METS[$i] = $METS[$i];
						if($METS[$i]=~/\w/){push(@GME, $METS[$i]);}}			$met = join(";", @GME);}
						
        $loc = mangle_locations($stuff[19]);
		 @KEGS=(); 	if($stuff[20]=~/\w+/){	    @KEGS = ($stuff[20]=~/([^\;]+)/g);}

		 @RHEA=(); 	if($stuff[21]=~/\w+/){	    @RHEA = ($stuff[21]=~/(\d+)/g);}


	## GET REACTIONS - setting sort JIC want to switch to top hit - may be neccessary if excessive/general function annotations
	#just toss a "last;" in the if statement

	%BRXN=(); %RRXN=(); %KRXN=();
	# MATCH RHEA/KEGG/BIOCYC ECs data with uniprot IDs to get the related rxn

	# changed sort approach !!!test 
	foreach my $rxn (sort{$b<=>$a} keys %{$UPID_BRXN{$upid}}){ 	if($rxn=~/\w/){$BRXN{$rxn}++; }}
	foreach my $rxn (sort{$b<=>$a} keys %{$UPID_RRXN{$upid}}){ 	if($rxn=~/\w/){$RRXN{$rxn}++; }}
	foreach my $rxn (sort{$b<=>$a} keys %{$UPID_KRXN{$upid}}){ 	if($rxn=~/\w/){$KRXN{$rxn}++; }}

	foreach my $ec 	(@ECS){ #match RHEA/KEGG/BIOCYC ECs data with uniprot ECs to get the related rxn 
		foreach my $rxn (sort{$b<=>$a} keys %{$EC_BRXN{$ec}}){  if($rxn=~/\w/){$BRXN{$rxn}++; }}
		foreach my $rxn (sort{$b<=>$a} keys %{$EC_RRXN{$ec}}){  if($rxn=~/\w/){$RRXN{$rxn}++; }}
		foreach my $rxn (sort{$b<=>$a} keys %{$EC_KRXN{$ec}}){	if($rxn=~/\w/){$KRXN{$rxn}++; }}}

	foreach my $mon (@MON){#biocyc monomers
                foreach my $rxn (sort{$b<=>$a} keys %{$MONO_BRXN{$mon}}){if($rxn=~/\w/){$BRXN{$rxn}++; }}}
	foreach my $kg (@KEGS){#kegg genes
		foreach my $rxn (sort{$b<=>$a} keys %{$KGEN_KRXN{$kg}}){if($rxn=~/\w/){$KRXN{$rxn}++; }}}
	foreach my $rxn (@RHEA){ 					if($rxn=~/\w/){$RRXN{$rxn}++; }}

	@BIO=(); @RHE=(); @KEG=();
	foreach my $rxn (sort{$BRXN{$b}<=>$BRXN{$a}} keys %BRXN){ push(@BIO,$rxn); }
	foreach my $rxn (sort{$RRXN{$b}<=>$RRXN{$a}} keys %RRXN){ push(@RHE,$rxn); }
	foreach my $rxn (sort{$KRXN{$b}<=>$KRXN{$a}} keys %KRXN){ push(@KEG,$rxn); }
	$kegg=join(";",@KEG);
	$rhea=join(";",@RHE);
	$bioc=join(";",@BIO);

	#COMBINE AND CLEAN PROT INFO
	@FIN=();
	$FIN[0]=$ur100;		$UR100_INFO{$ur100}{0}++;  # TODO: need this?
	$FIN[1]=$ur90;		$UR100_INFO{$ur100}{1}=$ur90;
	$FIN[2]=$name;		$UR100_INFO{$ur100}{2}=$UR100_NAME{$ur100};
	$FIN[3]=$plen;		$UR100_INFO{$ur100}{3}=$UR100_LEN{$ur100};

	#coordinates or ambig
	$FIN[4]=$sig;
	$FIN[5]=$tms;
	$FIN[6]=$dna;		for my $i (4..6){ $UR100_INFO{$ur100}{$i}=$FIN[$i]; }
	$FIN[7]=$tid;		if($tid=~/^\d+$/){$UR100_INFO{$ur100}{7}{$tid}++;}  # TODO: overwrite existing??
	
	#basic Function IDs
	$FIN[8]=$met;
	$FIN[9]=$loc;
	$FIN[10]=$tcp;
	$FIN[11]=$cog;
	$FIN[12]=$pfa;
	$FIN[13]=$tig;
	$FIN[14]=$gos;
	$FIN[15]=$ipr;
	$FIN[16]=$ecs;

	#reactions
	$FIN[17]=$kegg;
	$FIN[18]=$rhea;
	$FIN[19]=$bioc;

	#get ur100/ur90 prot-func counts
	for my $i (8..19){
		@IDS=();
		@IDS=split(";", $FIN[$i]);
		%seen=(); @GIDS=();
		@IDS = grep{ !$seen{$_}++ } @IDS;
		foreach my $id (@IDS){ if($id!~/\w/){next;} push(@GIDS,$id); $UR100_INFO{$ur100}{$i}{$id}++; $UR90_INFO{$ur90}{$i}{$id}++;}
		@IDS = nsort(@GIDS);
		$FIN[$i]=join(";", @IDS);
	}
	
	#output cleaned uniprot
	$out=join("\t",@FIN);
	print OUTPINT "$upid\t$out\n";
	if(progress($on)) {$time=localtime; print "on $on time $time out upids, skipped $skipped\n"; } $on++;
}
close INUP;
undef(%EC_RRXN);
undef(%UPID_RRXN);
undef(%EC_KRXN);
undef(%UPID_KRXN);
undef(%EC_BRXN);
undef(%UPID_BRXN);
undef(%MONO_BRXN);
undef(%KGEN_KRXN);
undef(%UPID_UR100);
undef(%UPID_UR90);
undef(%UR100_TCDB);
close(OUTPINT);
$num_info_recs = keys %UR100_INFO;
print "DONE on $on uniprots, skipped: $skipped; #INFO=$num_info_recs\n\n";



#COMPILE AND OUTPUT UNIREF
$on=0; $time=localtime; my $nolenskip = 0;
print "OUTPUT UNIREF100 $time\n";
#go thru all UniRef100s
print OUTREF "UR100\tUR90\tName\tLength\t";
print OUTREF "SigPep\tTMS\tDNA\tTaxonId\tMetal\tLoc\t";
print OUTREF "TCDB\tCOG\tPfam\tTigr\tGene_Ont\tInterPro\tECs\t";
print OUTREF "kegg\trhea\tbiocyc\n";
my (@OUT, $nc, %extra_names, $ids);
foreach my $ur100 (keys %UR100_LEN){
	@OUT=();
	$OUT[0]=$ur100;
	$ur90=$UR100_UR90{$ur100};
	$OUT[1]=$ur90;
	$OUT[2]=$UR100_NAME{$ur100};
	$OUT[3]=$UR100_LEN{$ur100};

	$nc=0;
	for my $i (4..6){ $OUT[$i] = $UR100_INFO{$ur100}{$i};}
        %extra_names = ();
	for my $i (7..19){
		@IDS=();
		foreach my $id (keys %{$UR100_INFO{$ur100}{$i}}){ if($id=~/\w/){push(@IDS,$id);}}
                if ($ur90) {
		    if($IDS[0]!~/\w/){ foreach my $id (keys %{$UR90_INFO{$ur90}{$i}}){ if($id=~/\w/){push(@IDS,$id);}}}
                }
                %seen=(); @IDS = grep{ !$seen{$_}++ } @IDS;

		foreach my $id (@IDS){	
			if($OUT[2]=~/(HYPOTHETICAL|UNCHARACTERIZED|UNDESCRIBED|UNKNOWN|PUTATIVE|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)/){
				if($FUNC_NAMES{$id}!~/(HYPOTHETICAL|UNCHARACTERIZED|UNDESCRIBED|UNKNOWN|PUTATIVE|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)/){
					if($FUNC_NAMES{$id}=~/\w/ && $nc<4 ){ $extra_names{$FUNC_NAMES{$id}}++; $nc++; }  # TODO: split func names by semi-colon?
		}	}	}
		
		@IDS = nsort(@IDS);
		$ids=join(";",@IDS);
		if($ids !~/\w/){$ids='';}
		$OUT[$i]=$ids;
	}
        if (keys %extra_names) {
            $OUT[2]=join(';', sort(keys %extra_names)).";".$OUT[2];
        }
	$out=join("\t", @OUT);
	print OUTREF "$out\n";
	if(progress($on)) {$time=localtime; print "on $on time $time ur100 $ur100 name $OUT[2]\n";} $on++;
}
$time = localtime;
print "DONE $on lines -- $time -- skipped $nolenskip\n\n";
undef(%UR100_LEN);
undef(%UR100_NAME);
undef(%UR100_INFO);



### SUBROUTINES ###
sub CleanNames{
	my @GREEKS = ("α","β","γ","δ","ε","ζ","η","θ","ι","κ","λ","μ","ν","ξ","ο","π","ρ","ς","σ","τ","υ","φ","χ","ψ","ω");
	my @GREEKL = ("ALPHA","BETA","GAMMA","DELTA","EPSILON","ZETA","ETA","THETA","IOTA","KAPPA","LAMBDA","MU","NU","XI","OMICRON","PI","RHO","SIGMA","SIGMA","TAU","UPSILON","PHI","CHI","PSI","OMEGA");
        my $nameX = $_[0];
        #remove junk punctuation/standardize
        my $sta=0; my $end=1;
        while($end ne $sta){
                $sta=$nameX;
                #swap greek symbols for text
                for my $g (0..$#GREEKL){ #fix pathbank and other greek symbols
                        if($nameX =~/($GREEKS[$g])/){
                                $nameX =~ s/$GREEKS[$g]/$GREEKL[$g]/g;
                }       }
                $nameX =~ s/\%2B(\d*)/$1\+/g;   #fix html +/- code (HMDB db)
                $nameX =~ s/\%2D(\d*)/$1\-/g;   #fix html +/- code (HMDB db)
                $nameX =~ s/(ARROW|STEREO|RIGHT|LEFT|\-)*\&/\&/g; #fix html +/- code (rhea)
                $nameX =~ s/\&\w+\;\/*//g; #fix html +/- code (rhea)
                $nameX =~ s/\s+/_/g;
                $nameX =~ s/[^\w\-\+]+/_/g;
                $nameX =~ s/\_\+|\+\_/\+/g;
                $nameX =~ s/\_\-|\-\_/\-/g;
                $nameX =~ s/\-+/\-/g;
                $nameX =~ s/\++/\+/g;
                $nameX =~ s/\++\-+|\-+\++/\+/g;
                $nameX =~ s/\_+/\_/g;
                $nameX =~ s/(^[\_\W]+|[\_\W]+$)//g;

                #clear out junk descriptors
                $nameX =~ s/^(LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)//g;
                $nameX =~ s/^(HYPOTHETICAL|UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)//g;
                $nameX =~ s/[\b\_](LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)[\b\_]/\_/g;
                $nameX =~ s/[\b\_](HYPOTHETICAL|UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)[\b\_]/\_/g;
                $nameX =~ s/[\b\_](LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)$//g;
                $nameX =~ s/[\b\_](HYPOTHETICAL|UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)$//g;
                $end=$nameX;
        }
        return($nameX);
}


sub fix_taxonomy{
	my $namey = $_[0];
	
	#fix the species/strain name
	$namey =~ s/([A-Z]+)\s(PROTEOBACTER)(IA|IUM)/$1$2$3/;
	$namey =~ s/\bPROPIONIBACTERIUM/CUTIBACTERIUM/g;
	$namey =~ s/\bLEPIDOSAURIA/SAURIA/g;
	$namey =~ s/ENDOSYMBIONT.OF\s+/ENDOSYMBIONT-/;
	$namey =~ s/COMPOSITE.GENOME.*//;
	$namey =~ s/MARINE.GROUP.(\w+)/MARINE-GROUP-$1/;
	$namey =~ s/\s+METAGENOME//;
	$namey =~ s/OOMYCETES/OOMYCOTA/;
	$namey =~ s/LILIOPSIDA/MAGNOLIOPSIDA/;
	$namey =~ s/^NR[^A-Z]//;	
	$namey =~ s/.*INCERTAE.SEDIS.*//;
	$namey =~ s/\_(PHYLUM|CLASS|ORDER|FAMILY|GENUS)[\b\_].*/\_/;
	$namey =~ s/ENRICHMENT.CULTURE.CLONES*|ENRICHMENT.CULTURES*//;
	$namey =~ s/\_(SENSU\_LATO|AFF|GEN|CF)\_/\_/g;
	$namey =~ s/^(SENSU\_LATO|AFF|GEN|CF)\_//g;
	$namey =~ s/\b(SENSU\_LATO|AFF|GEN|CF)\_//g;

	#remove ambiguous junk
        $namey =~ s/^(LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)//g;
        $namey =~ s/^(UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)//g;
        $namey =~ s/[\b\_](LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)[\b\_]/\_/g;
        $namey =~ s/[\b\_](UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)[\b\_]/\_/g;
        $namey =~ s/[\b\_](LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)$//g;
        $namey =~ s/[\b\_](UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)$//g;
	
	return($namey);
}

sub mangle_locations {
    # the uniprot subcellular location field
    #
    # This function picks up words and spaces following the marker "SUBCELLULAR
    # LOCATION: ", of which there can be multiple per record.  It misses when
    # the marker is immediately followed by some "[...]" bracketed stuff
    # (TODO.)  Cases where the marker is followed by only a "Note=..." are
    # skipped.  There is quite a bit more going on in the field.
    #
    my @found = m/SUBCELLULAR LOCATION: ([\w ]+\w)/g;
    my %locs;
    foreach my $i (@found) {
        if ($i =~ m/^Note$/) { next; }  # marker was immediately followed by "Note=..."
        $locs{$i} = 1;
    }
    return join(';', sort(keys %locs));
}

sub progress {
    my $on = shift;
    return 1 if (grep {$_ == $on} (1, 10, 100, 1000, 10000, 100000, 1000000));
    return 1 if ($on % 10000000 == 0);
    return 0;
}
