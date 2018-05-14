#!/usr/bin/perl -w

# mask_uppercase_dna.pl
# Purpose: mask uppercase nucleotides as 'X'.

# Erich Schwarz, 10/18/04, <emsch@its.caltech.edu>

unless (($ARGV[0]) or (! $ARGV[1])) { 
    print "Format: ./mask_uppercase_dna.pl [input ASCII text file]\n";
    die;
}

$input = $ARGV[0];
chomp($input);
$output = $input . ".masked";

open (INPUT, "$input") || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    unless ($input_line =~ /^>/) 
    {
        $input_line =~ s/[A-Z]/X/g;
    }
    print OUTPUT "$input_line\n";
}

close INPUT;
close OUTPUT;
