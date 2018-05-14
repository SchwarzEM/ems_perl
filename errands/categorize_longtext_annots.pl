#!/usr/bin/perl -w

# categorize_longtext_annots.pl
# by Erich Schwarz, 6/12/02.

# Purpose: classify annotations.

print "LongText .ace file to be categorized?: ";
chomp($input = <STDIN>);
$rough_output = $input . ".rough_categories";
$output = $input . ".categories";
$summary = $input . ".summary";

$ref_testfile = $input . ".ref_testfile";

$rough_genelist = $input . ".rough_genelist";
$genelist = $input . ".genelist";

open (INPUT, "$input") || die "Couldn't open input file. $!\n";
open (ROUGH_OUTPUT, ">$rough_output") || die "Couldn't open rough output file. $!\n";
open (ROUGH_GENELIST, ">$rough_genelist") || die "Couldn't open rough gene list. $!\n";

open (REF_TESTFILE, ">$ref_testfile") || die "Couldn't open test reference list. $!\n";

# Disease (cgc4103) reference:    Paper_evidence "[cgc4103]"
# Disease (cgc4637) reference:    Paper_evidence "[cgc4637]"

$disease_filter_1 = "cgc4103";

$disease_filter_2 = "cgc4637";

# DnaJ references:                PMID_evidence "10318904"
#                                 PMID_evidence "11230128"

$dnaj_filter = "pmid10318904";

# Prion references:               PMID_evidence "10611975"
#                                 PMID_evidence "11050225" <-- key reference
#                                 PMID_evidence "11447696"
#                                 PMID_evidence "11685242"
#                                 PMID_evidence "11782551"
#                                 PMID_evidence "11832240"

$prion_filter = "pmid11050225";

# Apoptosis refs:                 PMID_evidence "11181990"

$apoptosis_filter = "pmid11181990";

# MYST acetyltransferase refs:    Paper_evidence "[cgc4587]"

$myst_filter = "cgc4587";

# Pred. mito. prots. ref:         PMID_evidence "11035803"

$mito_pred_filter = "pmid11035803";

# F-box ref:                      PMID_evidence "11178263"

$f_box_filter = "pmid11178263";

$gene_name = "";
@paper_references      = "";
$entry_name            = "";

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);

    if ($input_line =~ /LongTextEnd/) 
    {
        @paper_references      = "";
        $entry_name            = "";
        $gene_name             = "";
    }

    elsif ($input_line =~ /^Locus : \"(\S+)\"/) 
    {
        $gene_name = $1;
    }

    elsif ($input_line =~ /^Sequence : \"(\S+)\"/) 
    {
        $entry_name = $1;
        if ($entry_name =~ /(.+\..*\d+)\D/) 
        {
            $gene_name = $1;
        }
        else 
        {
            $gene_name = $entry_name;
        }
    }

    elsif ($input_line =~ /Paper_evidence \".(\S+).\"/)
    {
        $indiv_reference = $1;
        unless ($indiv_reference =~ /^cgc\d{1,}$/) 
        {
            print "Misparsed reference for $gene_name: $indiv_reference\n";
            die;
        }
        chomp ($indiv_reference);
        print REF_TESTFILE "$gene_name\t$indiv_reference\n";

        push (@paper_references, $indiv_reference);
    }

    elsif ($input_line =~ /PMID_evidence \"(\d+)\"/)
    {
        $indiv_reference = "pmid" . "$1";
        unless ($indiv_reference =~ /^pmid\d{5,}$/)
        {
            print "Misparsed reference for $gene_name: $indiv_reference\n";
            die;
        }
        chomp ($indiv_reference);
        print REF_TESTFILE "$gene_name\t$indiv_reference\n";

        push (@paper_references, $indiv_reference);
    }

    elsif ($input_line =~ /^LongText :/) 
    {
        @sorted_paper_refs = sort (@paper_references);
        $ordered_refs = join ":", @sorted_paper_refs;

        if ($ordered_refs =~ /$disease_filter_1/) 
        {
            print ROUGH_OUTPUT "$gene_name\tDisease_homology\n";
            print ROUGH_GENELIST "$gene_name\n";
        }

        if ($ordered_refs =~ /$disease_filter_2/) 
        {
            print ROUGH_OUTPUT "$gene_name\tDisease_homology\n";
            print ROUGH_GENELIST "$gene_name\n";
        }

        if ($ordered_refs =~ /$dnaj_filter/)
        {
            print ROUGH_OUTPUT "$gene_name\tDnaJ_domain\n";
            print ROUGH_GENELIST "$gene_name\n";
        }

        if ($ordered_refs =~ /$prion_filter/)
        {
            print ROUGH_OUTPUT "$gene_name\tPrion_like\n";
            print ROUGH_GENELIST "$gene_name\n";
        }

        if ($ordered_refs =~ /$apoptosis_filter/)
        {
            print ROUGH_OUTPUT "$gene_name\tApoptosis_homolog\n";
            print ROUGH_GENELIST "$gene_name\n";
        }

        if ($ordered_refs =~ /$myst_filter/)
        {
            print ROUGH_OUTPUT "$gene_name\tMYST_acetyltransferase\n";
            print ROUGH_GENELIST "$gene_name\n";
        }

        if ($ordered_refs =~ /$mito_pred_filter/)
        {
            print ROUGH_OUTPUT "$gene_name\tPred_mitoprot\n";
            print ROUGH_GENELIST "$gene_name\n";
        }

        if ($ordered_refs =~ /$f_box_filter/)
        {
            print ROUGH_OUTPUT "$gene_name\tF-box_protein\n";
            print ROUGH_GENELIST "$gene_name\n";
        }

# Delete all references associated with stereotypical annotations.

        $ordered_refs =~ s/cgc4103//g;
        $ordered_refs =~ s/cgc4637//g;
        $ordered_refs =~ s/pmid10318904//g;
        $ordered_refs =~ s/pmid11230128//g;
        $ordered_refs =~ s/pmid10611975//g;
        $ordered_refs =~ s/pmid11050225//g;
        $ordered_refs =~ s/pmid11447696//g;
        $ordered_refs =~ s/pmid11685242//g;
        $ordered_refs =~ s/pmid11782551//g;
        $ordered_refs =~ s/pmid11832240//g;
        $ordered_refs =~ s/pmid11181990//g;
        $ordered_refs =~ s/cgc4587//g;
        $ordered_refs =~ s/pmid11035803//g;
        $ordered_refs =~ s/pmid11178263//g;

        $ordered_refs =~ s/\:{2,}/\:/g;

# If any references are left, then, by definition, this was manually annotated.

        if (($ordered_refs =~ /cgc\d+/) || ($ordered_refs =~ /pmid\d+/)) 
        {
            print ROUGH_OUTPUT "$gene_name\tManual_annotation\n";
            print ROUGH_GENELIST "$gene_name\n";
        }
    }
}

close INPUT;
close ROUGH_OUTPUT;

system "sort $rough_output > $output";
system "rm $rough_output";

system "sort $rough_genelist > $genelist";
system "rm $rough_genelist";

# $disease_filter_1 = "cgc4103";
# $disease_filter_2 = "cgc4637";
# $dnaj_filter = "pmid10318904.+pmid11230128";
# $prion_filter = "pmid10611975.*pmid11050225.*pmid11447696.*pmid11685242.*pmid11782551.*pmid11832240";
# $apoptosis_filter = "pmid11181990";
# $myst_filter = "cgc4587";
# $mito_pred_filter = "pmid11035803";

open (SUMMARY, ">$summary") || die "Couldn't open summary. $!\n";

$rough_gene_count = `wc -l $genelist`;
if ($rough_gene_count =~ /^\s+(\d+)\s/) 
{
    $total_gene_count = $1;
}
print SUMMARY "$total_gene_count genes were annotated.\n";
print SUMMARY "\n";

# Disease_homology
# DnaJ_domain
# Prion_like
# Apoptosis_homolog
# MYST_acetyltransferase
# Pred_mitoprot
# F-box_protein

system "grep Manual_annotation $output > manual_count_scratchfile";
$count = `wc -l manual_count_scratchfile`;
system "rm manual_count_scratchfile";
$count =~ s/\w+_count_scratchfile//g;
$count =~ s/\s+//g;
print SUMMARY "$count\tgenes were manually annotated.\n";

system "grep Disease_homology $output > disease_count_scratchfile";
$count = `wc -l disease_count_scratchfile`;
system "rm disease_count_scratchfile";
$count =~ s/\w+_count_scratchfile//g;
$count =~ s/\s+//g;
print SUMMARY "$count\tgenes have disease homologies noted in cgc4103 or cgc4637.\n";

system "grep DnaJ_domain $output > dnaj_count_scratchfile";
$count = `wc -l dnaj_count_scratchfile`;
system "rm dnaj_count_scratchfile";
$count =~ s/\w+_count_scratchfile//g;
$count =~ s/\s+//g;
print SUMMARY "$count\tgenes encode DnaJ-domain proteins.\n";

system "grep Prion_like $output > prion_count_scratchfile";
$count = `wc -l prion_count_scratchfile`;
system "rm prion_count_scratchfile";
$count =~ s/\w+_count_scratchfile//g;
$count =~ s/\s+//g;
print SUMMARY "$count\tgenes encode prion-like proteins.\n";

system "grep Apoptosis_homolog $output > apoptosis_count_scratchfile";
$count = `wc -l apoptosis_count_scratchfile`;
system "rm apoptosis_count_scratchfile";
$count =~ s/\w+_count_scratchfile//g;
$count =~ s/\s+//g;
print SUMMARY "$count\tgenes encode apoptosis-related proteins.\n";

system "grep MYST_acetyltransferase $output > myst_count_scratchfile";
$count = `grep MYST_acetyltransferase $output | wc -l`;
system "rm myst_count_scratchfile";
$count =~ s/\w+_count_scratchfile//g;
$count =~ s/\s+//g;
print SUMMARY "$count\tgenes encode MYST acetyltransferase proteins.\n";

system "grep Pred_mitoprot $output > mito_count_scratchfile";
$count = `wc -l mito_count_scratchfile`;
system "rm mito_count_scratchfile";
$count =~ s/\w+_count_scratchfile//g;
$count =~ s/\s+//g;
print SUMMARY "$count\tgenes encode predicted mitochondrial proteins.\n";

system "grep F-box_protein $output > f_box_count_scratchfile";
$count = `grep F-box_protein $output | wc -l`;
system "rm f_box_count_scratchfile";
$count =~ s/\w+_count_scratchfile//g;
$count =~ s/\s+//g;
print SUMMARY "$count\tgenes encode F-box proteins.\n";
print SUMMARY "\n";

close SUMMARY;
