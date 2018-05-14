#!/usr/bin/perl -w

# clean_up_pre_FASTA.pl
# Erich Schwarz, 6/28/02

# Purpose: clean up various FASTA precursor files.

print "What file do you want to clean? ";
$infile_name = <STDIN>;
chomp($infile_name);

# Safety and monitoring:

print "You want to clean $infile_name, right? ";
$answer = <STDIN>;
chomp $answer;
if ($answer ne "yes") 
{
    print "Anything but a yes means a no.\n";
    die;
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
    if ($line_to_output =~ /^>/) 
    {
        print OUTFILE $line_to_output;
    }
    else 
    {
        $line_to_output =~ s/\d//g;
        $line_to_output =~ s/ //g;
        $line_to_output =~ s/\-//g;
        unless ($line_to_output =~ /^\n$/) 
        {
            print OUTFILE $line_to_output;
        }
    }
}
