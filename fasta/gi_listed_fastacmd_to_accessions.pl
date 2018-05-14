#!/usr/bin/perl -w

# gi_listed_fastacmd_to_accessions.pl
# Erich Schwarz <emsch@its.caltech.edu>, 11/26/02.
# Purpose: given a *.fastacmd FASTA file, from a now-obsolete gi number list, extract (one hopes) up-to-date accession numbers.

print "input file: ";
$input = <STDIN>;
chomp($input);
$rough_output = $input . ".rough_output";
$output       = $input . ".accs";

open (INPUT, "$input") || die;
open (ROUGH_OUTPUT, ">$rough_output") || die;

# Example of input lines:

# >gi|10179844|gb|AAG13909.1|AF263245_5 mycarose O-acyltransferase [Micromonospora megalomicea subsp. nigra]
# >gi|1076118|pir||JC4001 macrolide 3-O-acyltransferase (EC 2.3.1.-) - Streptomyces sp >gi|551631 ... [Streptomyces thermotolerans]
# >gi|10956333|ref|NP_052782.1| pXO1-86 [Bacillus anthracis] ...

# Example of desired output:

# ...

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    if ($input_line =~ /^>gi\|\d+\|\w+\|+(\S+)\.\d+\|+/)
    {
        print ROUGH_OUTPUT "$1\n";
    }
    elsif ($input_line =~ /^>gi\|\d+\|\w+\|+(\S+)\.\d+\s+/)
    {
        print ROUGH_OUTPUT "$1\n";
    }
    elsif ($input_line =~ /^>gi\|\d+\|\w+\|+(\S+)\|+/)
    {
        print ROUGH_OUTPUT "$1\n";
    }
    elsif ($input_line =~ /^>gi\|\d+\|\w+\|+(\S+)\s+/)
    {
        print ROUGH_OUTPUT "$1\n";
    }
}

close INPUT;
close ROUGH_OUTPUT;

system "sort $rough_output | uniq - > $output";
system "rm $rough_output";

