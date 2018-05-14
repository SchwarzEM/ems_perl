#!/usr/bin/perl
use strict;

# fix_MSWord_txt.pl
# Erich Schwarz <emsch@its.caltech.edu>, 9/27/02.
# Purpose: fix '\r in place of \n' bug of MS word 'text' files.

print "Input file?  ";

my $input = <STDIN>;
chomp ($input);

my $output = $input . ".output";
my $protein_record = "";

open (INPUT, "$input")    || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    my $input_line = $_;
    my @protein_records = split /\r/, $input_line;
    foreach $protein_record (@protein_records) 
    {
        print OUTPUT "$protein_record\n";
    }
}
