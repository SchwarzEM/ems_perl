#!/usr/bin/perl -w

# extract_GOified_loci.pl
# Erich Schwarz, 1/28/02

# goal: quick hack to get locus counts from slices of gene_association.wb

print "Name of file that has gene_association.wb lines from which to count loci? ";
$input = <STDIN>;
chomp ($input);
$rough_output = $input . ".rough_output";
$sorted_output = $input . ".sorted_output";

open (INPUT, "$input") || die "Can't open $input. $!\n";
open (ROUGH_OUTPUT, ">$rough_output") || die "Can't open $rough_output. $!\n";

while (<INPUT>) 
{
    $input_line = $_;
    if ($input_line =~ /^WB\t\S+\t(\S+)\t/) 
    { 
        print ROUGH_OUTPUT "$1\n";
    }
}

close INPUT;
close ROUGH_OUTPUT;

system ("sort $rough_output | uniq - > $sorted_output");
system ("rm $rough_output");
