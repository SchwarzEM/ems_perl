#!/usr/bin/perl -w

# clean_prom_DNA.pl
# Erich Schwarz <emsch@its.caltech.edu>, 5/20/03

# Purpose: clean up various DNA text files.

print "File to be cleaned: ";
chomp($infile_name = <STDIN>);

$outfile_name = $infile_name . ".output";
open (INFILE, $infile_name);
open (OUTFILE, ">$outfile_name");

print "The cleaned version of $infile_name will be $outfile_name.\n";

while (<INFILE>) 
{
    chomp($input_line = $_);
    if ($input_line =~ /^>/) 
    {
        print OUTFILE "$input_line\n";
    }
    elsif ($input_line =~ /\S+/) 
    {
        $input_line =~ s/\d//g;
        $input_line =~ s/ //g;
        while ($input_line =~ /^([\S]{68})(.*)/)
        {
            print OUTFILE "$1\n";
            $input_line = $2;
        }
        if ($input_line =~ /^[\S]{1,67}$/) 
        {
            print OUTFILE "$input_line\n";
        }
    }
}
