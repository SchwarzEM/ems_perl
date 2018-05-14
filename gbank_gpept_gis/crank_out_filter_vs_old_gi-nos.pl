#!/usr/bin/perl -w
# crank_out_filter_vs_old_gi-nos.pl
# Erich Schwarz, 11/17/02.

# Make a script from lines like this:
# -rw-rw-r--    1 schwarz  schwarz    107099 Nov  6 00:52 14574446_aa_2-66.psi-blast-x30.nr_1e-04.align-6

print "\n";
print "  This will generate a script for filtering out pre-detected \n";
print "  gi numbers from psi-blast searches.\n";
print "\n";

print "User-supplied list of prior gi numbers will be?  ";
chomp($prior_list = <STDIN>);

$date="";
system 'date +"%s" > foobar.date.tmp';
open DATE, "foobar.date.tmp";
$date=<DATE>;
chomp($date);
system 'rm foobar.date.tmp';

$input_list_name=$date."_input_list";
print "Precursor to extraction script is $input_list_name\n";

system "ls -l *align-6 > $input_list_name";

$output_script=$date.".output_filter_script";

open INPUT_LIST, "$input_list_name";
open OUTPUT_SCRIPT, ">$output_script";

print "Extraction script is $output_script\n";

while (<INPUT_LIST>) 
{
    if ($_ =~ /.+\s+(\d+_aa_.+align-6)/) 
    {
        $input_psi_blast_file = $1;
        print OUTPUT_SCRIPT "grep -v -f $prior_list ";
        print OUTPUT_SCRIPT "$input_psi_blast_file ";
        print OUTPUT_SCRIPT '> ';
        print OUTPUT_SCRIPT "$input_psi_blast_file";
        print OUTPUT_SCRIPT '.filtered';
        print OUTPUT_SCRIPT ";\n";
    }
}

close INPUT_LIST;
close OUTPUT_SCRIPT;

system "rm $input_list_name";
system "chmod +x $output_script";

