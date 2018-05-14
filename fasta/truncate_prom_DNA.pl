#!/usr/bin/perl -w

# truncate_prom_DNA.pl
# Erich Schwarz <emsch@its.caltech.edu>, 10/5/03.

# Purpose: take long DNA files with a fixed-length promoter + extraneous seq., and discard latter.

print "What DNA sequence file do you want to trim? ";
chomp ($infile_name = <STDIN>);

print "How many nucleotides do you want to keep? ";
chomp ($length = <STDIN>);
$length =~ s/\s//g;
print "$length\n";

$outfile_name = $infile_name . ".cleaned_FASTA";
open (INFILE, $infile_name);
open (OUTFILE, ">$outfile_name");

$output_line = "";

while (<INFILE>) 
{
    chomp ($input_line = $_);
    if ($input_line =~ /^>/)
    {
        unless ($output_line eq "") 
        {
            if ($output_line =~ /^([\S]{$length})(.*)/)
            {
                $output_line = $1;
            }
            while ($output_line =~ /^([\S]{68})(.*)/)
            {
                print OUTFILE "$1\n";
                $output_line = $2;
            }
            if ($output_line =~ /^[\S]{1,67}$/)
            {
                print OUTFILE "$output_line\n";
            }
        }
        print OUTFILE "$input_line\n";
    }
    else 
    {
        $input_line =~ s/\d//g;
        $input_line =~ s/ //g;
        $input_line =~ s/\-//g;
        $input_line =~ s/\*//g;
#        if ($input_line =~ /\S/) 
#        {
            $output_line = ("$output_line" . "$input_line");
#        }
    }
}

            
if ($output_line =~ /^([\S]{$length})(.*)/)
{
    $output_line = $1;
}
while ($output_line =~ /^([\S]{68})(.*)/)
{
    print OUTFILE "$1\n";
    $output_line = $2;
}
if ($output_line =~ /^[\S]{1,67}$/)
{
    print OUTFILE "$output_line\n";
}

close INFILE;
close OUTFILE;
