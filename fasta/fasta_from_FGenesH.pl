#!/usr/bin/perl -w

# fasta_from_FGenesH.pl
# Erich Schwarz, 10/19/03

# Purpose: extract FASTA (both protein and nucleotide) from an FGenesH output, either manually or in batch mode.

if ($#ARGV != 0) 
{
    print "FGenesH file to extract protein and DNA FASTA subsets from? ";
    $input_fgenesh = <STDIN>;
} 
else 
{
    $input_fgenesh = $ARGV[0];
}
chomp($input_fgenesh);

if ($input_fgenesh =~ /^(\d+)\D+(Contig\d+)\D+/) 
{
    $serial_no = $1;
    $contig_no = $2;
}
elsif ($input_fgenesh =~ /^(\d+)\D+/) 
{
    $serial_no = $1;
    $contig_no = "";
}
else
{
    $serial_no = $input_fgenesh;
    $contig_no = "";
}

$output_prot = $input_fgenesh . ".prot.tfa";
$output_nucl = $input_fgenesh . ".nucl.tfa";

open (INPUT_FGENESH,      "$input_fgenesh") || die "Can't open $input_fgenesh. $!\n";
open (OUTPUT_PROTEIN,     ">$output_prot") || die "Can't open $output_prot. $!\n";
open (OUTPUT_NUCLEOTIDES, ">$output_nucl") || die "Can't open $output_nucl. $!\n";

$reading_protein     = "no";
$reading_nucleotides = "no";
$reading_ok = "no";

while (<INPUT_FGENESH>) 
{
    $fasta_line = $_;

    if ($fasta_line =~ /FGENESH 2.0 Prediction of potential genes in C_elegans genomic DNA/)
    {
        $reading_ok = "yes";
    }
    elsif (($fasta_line =~ /^>FGENESH:/) && ($reading_ok eq "no")) 
    {
        print "$input_fgenesh was not generated from FGENESH 2.0 and C. elegans defaults.\n";
        die;
    }
    elsif ($fasta_line =~ /^>FGENESH:.mRNA.\s+(\d+)\s+(.*)/) 
    {
        print OUTPUT_NUCLEOTIDES ">$serial_no";
        print OUTPUT_NUCLEOTIDES "_";
        print OUTPUT_NUCLEOTIDES "$contig_no";
        print OUTPUT_NUCLEOTIDES "_";
        print OUTPUT_NUCLEOTIDES "$1";
        print OUTPUT_NUCLEOTIDES ".dna    $2\n";
        $reading_protein         = "no";
        $reading_nucleotides     = "yes";
    }
    elsif ($fasta_line =~ /^>FGENESH:\s+(\d+)\s+(.*)/)
    {    
        print OUTPUT_PROTEIN ">$serial_no";
        print OUTPUT_PROTEIN "_";
        print OUTPUT_PROTEIN "$contig_no";
        print OUTPUT_PROTEIN "_";
        print OUTPUT_PROTEIN "$1    $2\n";
        $reading_protein     = "yes";
        $reading_nucleotides = "no";
    }
    elsif ($reading_nucleotides eq "yes" && ($fasta_line =~ /^\w+/) )
    {
        print OUTPUT_NUCLEOTIDES "$fasta_line";
    }
    elsif ($reading_protein eq "yes" && ($fasta_line =~ /^\w+/) ) 
    {
        print OUTPUT_PROTEIN "$fasta_line";
    }
}

close INPUT_FGENESH;
close OUTPUT_PROTEIN;
close OUTPUT_NUCLEOTIDES;
