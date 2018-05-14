#!/usr/bin/perl -w

# fix_MS_returns.pl
# Erich Schwarz <emsch@its.caltech.edu>, 8/8/03
# Purpose: clean up "\r" in defective Word-to-text.

print "input file: ";
$input = <STDIN>;
chomp($input);

$output = $input . ".cleaned_Word_txt_one_hopes";

open (INPUT,  "$input")   || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    $input_line = $_;
    $input_line =~ s/\r/\n/g;
    print OUTPUT "$input_line";
}

close INPUT;
close OUTPUT;

