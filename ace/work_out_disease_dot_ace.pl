#!/usr/bin/perl -w

# work_out_disease_pep_text.pl
# Purpose: Generate reasonable plaintext from DiseasePep orthologies.

# One set of input lines -- homology/orthology file, default "diseasepep_output":

# A2M	ZK337.1A	PZP	
# AAAS	K04G11.4	WDR5	
# AASS	R02D3.1	AASS	1
# ABAT	K04D7.3	ABAT	1

# Second set of input lines -- sequence plus locus, default "protein_sequence_locus_WS81.txt":

# "WP:CE00015"	"B0464.1"	"drs-1"
# "WP:CE00044"	"C02F5.8"	"tsp-1"
# "WP:CE00072"	"C14B9.1"	"hsp-12.2"
# "WP:CE00077"	"C14B9.6a"	"gei-8"

# Another set of input lines -- informative FASTA headers, default "DiseasePep.headers":

# >MPO (LocusID: 4353) OMIM:606989 MYELOPEROXIDASE [HOMO SAPIENS] gi|2160397
# >ALDOB (LocusID: 229) OMIM:229600 ALDOLASE B [HOMO SAPIENS] gi|2160383
# >CP (LocusID: 1356) OMIM:117700 CERULOPLASMIN [HOMO SAPIENS] gi|1620909
# >PHKA2 (LocusID: 5256) OMIM:306000 PHOSPHORYLASE KINASE ALPHA SUBUNIT [HOMO SAPIENS] gi|1217901

print "Input homology/orthology file (default \"diseasepep_output\"): ";
$orthology_input = <STDIN>;
chomp($orthology_input);
unless ($orthology_input =~ /\s+/) 
{
    $orthology_input = "diseasepep_output";
}

print "Input protein-sequence-locus file (default \"protein_sequence_locus_WS85.txt\"): ";
$locus_input = <STDIN>;
chomp($locus_input);
unless ($locus_input =~ /\s+/)
{
    $locus_input = "protein_sequence_locus_WS85.txt";
}

print "Input FASTA with informative headers file (default \"DiseasePep.headers\"): ";
$headers_input = <STDIN>;
chomp($headers_input);
unless ($headers_input =~ /\s+/)
{
    $headers_input = "DiseasePep.headers";
}

$output = $orthology_input . ".pre-ace-annots";
$seq_list = $orthology_input . ".disease_ortho_seqlist";

open (ORTHO_INPUT, "$orthology_input") || die;

open (OUTPUT, ">$output") || die;
open (SEQLIST, ">$seq_list") || die;

while (<ORTHO_INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    if ($input_line =~ /^(.+)\s(.+)\s.+\s1$/) 
    {
        $human_locus = $1;
        $worm_seq = $2;
        $gene_name = $worm_seq;

        if ($worm_seq =~ /(\S+)([A-Z]+)$/) 
        {
            $gene_name = $1;
            $worm_seq_suffix = $2;
            $worm_seq_suffix =~ tr/A-Z/a-z/;
            $worm_seq = $gene_name . $worm_seq_suffix;
        }

        open (LOCUS_INPUT, "$locus_input") || die;
        while (<LOCUS_INPUT>) 
        {
            $locus_input_line = $_;
            chomp ($locus_input_line);
            if ($locus_input_line =~ /^\"WP.+\"\s+\"$worm_seq\"\s+\"(.+)\"/) 	# "WP:CE00077"  "C14B9.6a"      "gei-8"
            {
                $gene_name = $1;
            }
        }
        close LOCUS_INPUT;

        open (HEADERS_INPUT, "$headers_input") || die;
        while (<HEADERS_INPUT>) 
        {
            $headers_input_line = $_;
            chomp ($headers_input_line);
            if ($headers_input_line =~ />$human_locus\s+.+LocusID.+(OMIM.\d+)\s+(.+)\s+.HOMO SAPIENS.\s+gi.\d+/)    
#                                     # >MPO (LocusID: 4353) OMIM:606989 MYELOPEROXIDASE [HOMO SAPIENS] gi|2160397
            {
                $omim_entry      = $1;
                $human_gene_name = $2;
            }
        }
        close HEADERS_INPUT;

        print OUTPUT "Locus : \"$gene_name\"\n";
        print OUTPUT "Provisional_description \"$gene_name\" Person_evidence \"Schwarz EM\"\n";
        print OUTPUT "Provisional_description \"$gene_name\" Person_evidence \"Chan J\"\n";
        print OUTPUT "\n";
        print OUTPUT "LongText : \"$gene_name\"\n";
        print OUTPUT "\n";
        print OUTPUT "The $gene_name gene is orthologous to the human gene ";
        print OUTPUT "$human_gene_name ($human_locus; $omim_entry), ";
        print OUTPUT "which when mutated leads to $human_locus disease.\n";
        print OUTPUT "\n";
        print SEQLIST "$worm_seq\n";
    }
}

close ORTHO_INPUT;

close OUTPUT;
close SEQLIST;
