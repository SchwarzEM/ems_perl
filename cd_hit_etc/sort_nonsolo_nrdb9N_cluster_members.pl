#!/usr/bin/perl -w

# sort_nonsolo_nrdb9N_cluster_members.pl
# Erich Schwarz, 5/15/02.

# extract sheep versus goats from a marked-up nrdbN.pl (e.g. nrdb98.pl)
#     non-singleton "clusters" file.

# input:
# 3       141     441     *15638629  no_gene_name  acyl transferase-like protein  gi|15638629    Ustilago maydis
# 3       296     340     14456153  no_gene_name  putative acyl transferase  gi|14456153    Ustilago maydis


print "input file: ";
$input = <STDIN>;
chomp($input);
$output1 = $input . ".desired_gi_nos";
$output2 = $input . ".unwanted_gi_nos";

open (INPUT, "$input") || die;
open (OUTPUT1, ">$output1") || die;
open (OUTPUT2, ">$output2") || die;

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);

    if ($input_line =~ /\d+\s+\d+\s+\d+\s+\*(\d+)\s+/)
    {
        print OUTPUT1 "$1\n";
    }
    elsif ($input_line =~ /\d+\s+\d+\s+\d+\s+(\d+)\s+/)
    {
        print OUTPUT2 "$1\n";
    }
}

close INPUT;
close OUTPUT1;
close OUTPUT2;
