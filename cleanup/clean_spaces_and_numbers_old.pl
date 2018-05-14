#!/usr/bin/perl -w

# clean_spaces_and_numbers.pl
# Erich Schwarz, 2/21/01
# To clean spaces and numbers out of a text file,
#   e.g., to clean a DNA or protein sequence file
#   from GenBank/etc. into FASTA-able format.

print "What file do you want to clean? ";
$infile_name = <STDIN>;
chomp($infile_name);

# Safety and monitoring:

print "You want to clean $input_filename, right? ";
$answer = <STDIN>;
chomp $answer;
if ($answer ne "yes") 
{
    print "Anything but a yes means a no.";
    die
}
print "OK, chief.\n";

# Setting up cleanup:

$outfile_name = $infile_name . ".output";
open (INFILE, $infile_name);
open (OUTFILE, ">$outfile_name");

# One more check:

print "The cleaned version of $infile_name will be $outfile_name.\n";

while (<INFILE>) 
    {
        $line_to_output = $_;
        $line_to_output =~ s/\d//g;
        $line_to_output =~ s/ //g;
        print OUTFILE $line_to_output;
    }
