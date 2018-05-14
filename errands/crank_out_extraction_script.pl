#!/usr/bin/perl -w
# crank_out_extraction_script.pl
# Erich Schwarz <emsch@its.caltech.edu>, 11/17/02.

# Purpose: work with psi-blast outputs that have had past hits filtered out to detect more hits.
 
# Make a script from lines like this:
# -rw-rw-r--    1 schwarz  schwarz    107099 Nov  6 00:52 14574446_aa_2-66.psi-blast-x30.nr_1e-04.align-6

# Use user-supplied e-value and h-value so that varying sensitivies can be generated.

print "\n";
print "  This will generate a script for extracting all the\n";
print "  NUMBER_aa_STUFF_align-6 files in the current directory\n"; 
print "  with the script convg_blastpgp_m-6_to_germane_hits.pl.\n";
print "  Because convg_blastpgp_m-6_to_germane_hits.pl needs an\n";
print "  e-value and h-value supplied by the user for future\n";
print "  psi-blast searches, so does this script.\n";
print "\n";

print "User-supplied e-value will be?:  ";
$e_value=<STDIN>;
chomp ($e_value);

print "User-supplied h-value will be?:  ";
$h_value=<STDIN>;
chomp ($h_value);

$date="";
system 'date +"%s" > foobar.date.tmp';
open DATE, "foobar.date.tmp";
$date=<DATE>;
chomp($date);
system 'rm foobar.date.tmp';

$input_list_name=$date."_input_list";
print "Precursor to extraction script is $input_list_name\n";

system "ls -l *align-6 > $input_list_name";

$output_script=$date.".output_bulk_germane_script";

open INPUT_LIST, "$input_list_name";
open OUTPUT_SCRIPT, ">$output_script";

print "Extraction script is $output_script\n";

while (<INPUT_LIST>) 
{
    if ($_ =~ /.+\s+(\d+_aa_.+align-6)/) 
    {
        print OUTPUT_SCRIPT 'convg_blastpgp_m-6_to_germane_hits.pl ';
        print OUTPUT_SCRIPT "$1";
        print OUTPUT_SCRIPT '.filtered ';
        print OUTPUT_SCRIPT "$e_value ";
        print OUTPUT_SCRIPT "$h_value ";
        print OUTPUT_SCRIPT ";\n";
    }
}

close INPUT_LIST;
close OUTPUT_SCRIPT;

system "rm $input_list_name";

system "chmod +x $output_script";
