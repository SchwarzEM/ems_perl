#!/usr/bin/perl -w

# clean_up_PomPep.pl
# Erich Schwarz <emsch@its.caltech.edu>, 3/11/03.

# Purpose: clean up pompep, March 2003.

print "What file do you want to clean? ";
chomp ($infile_name = <STDIN>);

# Setting up cleanup:

$outfile_name = $infile_name . ".cleaned_FASTA";
open (INFILE, $infile_name);
open (OUTFILE, ">$outfile_name");

# One more check:

print "The cleaned version of $infile_name will be $outfile_name.\n";

while (<INFILE>) 
{
    chomp ($line_to_output = $_);
    if ($line_to_output =~ /^>(\S+)\s*.*$/) 
    {
        $first_half_of_header  = $1;
        $first_half_of_header  =~ s/;/    /g;
        $first_half_of_header  =~ s/,/    /g;
        print OUTFILE ">$first_half_of_header\n";
    }
    elsif ($line_to_output =~ /^>(\S+)/)
    {
        $header  = $1;
        $header  =~ s/;/    /g;
        $header  =~ s/,/    /g;
        print OUTFILE ">$header\n";
    } 
    else 
    {
        $line_to_output =~ s/\d//g;
        $line_to_output =~ s/ //g;
        $line_to_output =~ s/\-//g;
        $line_to_output =~ s/\*//g;
        unless ($line_to_output =~ /^\n$/) 
        {
            print OUTFILE "$line_to_output\n";
        }
    }
}

close INFILE;
close OUTFILE;
