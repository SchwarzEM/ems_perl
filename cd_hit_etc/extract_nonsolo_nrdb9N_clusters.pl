#!/usr/bin/perl -w 

# extract_nonsolo_nrdb9N_clusters.pl
# Erich Schwarz, 5/15/02

# Purpose: extract the subset of sequences from the "clusters"
#    output of nrdb90.pl (or nrdb98.pl, etc.) that in fact have
#    more than one member.

print "\n";
print "  The purpose of this script is to extract the subset of groups\n";
print "  from the \"clusters\" output of a nrdb9N.pl script (e.g., nrdb98.pl)\n";
print "  which, in fact, have more than one member per \'cluster\'.\n";
print "\n";
print "File to extract? (default: \"clusters\"): ";

$input = <STDIN>;
chomp ($input);
unless ($input =~ /\S+/) 
{
     $input = "clusters";
}

$output = $input . ".filtered_nrdb_cluster";

$input_line = "";
$print_line = "no";

open (INPUT_FILE, "$input")    || die "Can't open $input file. $!\n";
open (OUTPUT_FILE, ">$output") || die "Can't open $output file. $!\n";

while (<INPUT_FILE>) 
{
    $input_line = $_;
    chomp($input_line);
    if ($input_line =~ /Cluster \d+ with 1 members/)
    {
        $print_line = "no";
    }
    elsif ($input_line =~ /Cluster \d+ with \d+ members/) 
    {
        $print_line = "yes";
        print OUTPUT_FILE "\n";
        print OUTPUT_FILE "$input_line\n";
    }
    elsif ($print_line eq "yes") 
    {
        print OUTPUT_FILE "$input_line\n";
    }
}

close INPUT_FILE;
close OUTPUT_FILE;
