#!/usr/bin/perl -w

# Program: runon_to_clean_tfa.pl
#    Erich Schwarz, 4/24/03
# 
# Purpose: get endless-line sequence files to have "\n" after 68 characters of a line lacking ">".

if ($#ARGV != 0) 
{
    print "Required: input nrdb file!\n";
    print "What will input file be? ";
    chomp($infile = <STDIN>);
} 
else 
{
    chomp($infile = $ARGV[0]);
}

$outfile = ($infile . ".cleaned");
print "The input file is $infile; the output file $outfile\n";

open (INFILE, "$infile") || die "File $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outfile\n";


while (<INFILE>)
{
    chomp($input_line = $_);
    if ($input_line =~ /^>/)
    {
        print OUTFILE "$input_line\n";
    }
    else
    {
         while ($input_line =~ /^([\S]{68})(.*)/) 
         {
             print OUTFILE "$1\n";
             $input_line = $2;
             unless ($input_line =~ /^([\S]{68})(.*)/) 
             {
                 print OUTFILE "$input_line\n";
             }
         }
    }
}
