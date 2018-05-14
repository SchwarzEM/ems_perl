#!/usr/bin/perl -w

# quik_hak_03apr2002_A.pl

print "input file: ";
$input = <STDIN>;
chomp($input);
$output = $input . ".output";

open (INPUT, "$input") || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    $input_line =~ s/Phenotype\t/Detailed_description\t/;
    $input_line =~ s/Remark\t/Detailed_description\t/;
    print OUTPUT "$input_line\n";
}

close INPUT;
close OUTPUT;
