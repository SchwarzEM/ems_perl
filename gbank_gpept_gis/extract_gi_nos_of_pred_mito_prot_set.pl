#!/usr/bin/perl -w

# extract_gi_nos_of_pred_mito_prot_set.pl

# Wonky one-time script to extract gi numbers from 
# exported text of:
#     

# Sample of text to extract:
# 
# <http://www3.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=3877397&form=6&db=p&Dopt=g%5C>
# predicted using Genefinder; Similarity to Pea glycine cleavage system H
# protein (SW:GCSH_PEA); cDNA EST EMBL:D68820 comes from this gene; cDNA
# EST yk333a3.3 comes from this gene; cDNA EST yk333a3.5 comes from this gene
# P 3875395
# <http://www3.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=3875395&form=6&db=p&Dopt=g%5C>
# Similarity to B.subtilis CLPX protein, an ATP-dependent CLP protease
# (TR:E221213); cDNA EST EMBL:T01094 comes from this gene; cDNA EST
# EMBL:D64644 comes from this gene; cDNA EST EMBL:D67734 comes from this
# gene; cDNA EST yk504h6.3 comes from this gene;>
# K 3876393
# <http://www3.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=3876393&form=6&db=p&Dopt=g%5C>
# Similarity to Human 2-oxoisovalerate dehydrogenase (SW:ODBB_HUMAN)

print "File to extract? ";
$input_file = <STDIN>;
chomp ($input_file);
$pred_output_file = $input_file . ".theo-pred-mitoprot-gi-nos";
$known_output_file = $input_file . ".known-mitoprot-gi-nos";

open (INPUT_FILE, "$input_file") || die "Can't open input file $input_file. $! \n";
open (PRED_OUTPUT_FILE, ">$pred_output_file") || die "Can't open predicted output file $pred_output_file. $! \n";
open (KNOWN_OUTPUT_FILE, ">$known_output_file") || die "Can't open known output file $known_output_file. $! \n";

while (<INPUT_FILE>) 
{
    if ($_ =~ /^K (\d{4,12})\D+/) 
    {
        print KNOWN_OUTPUT_FILE "$1\n";
    }
    elsif ($_ =~ /^P (\d{4,12})\D+/)
    {
        print PRED_OUTPUT_FILE "$1\n";
    }
}
close INPUT_FILE;
close KNOWN_OUTPUT_FILE;
close PRED_OUTPUT_FILE;
