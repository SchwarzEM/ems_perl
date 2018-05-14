#!/usr/bin/perl -w

## Usage: extractgff.pl source feature strand name.gff name.fasta outputfile
##
## with source          source according to gff version 2 standard
##      feature         feature according to gff version 2 standard
##      strand          + or - for the strand
##      name.gff        gff filename. GFF format required.
##      name.fasta      fasta filename of sequence. FASTA format required.
##      outputfile      name of outputfile
##
## Example: extractgff.pl curated exon + CHROMOSOME_I.gff CHROMOSOME_I.dna outfile.txt
##
##          Searches CHROMOSOME_I.gff for the source keyword 'curated', the feature 'exon' and
##          the strand specification '+'. If found, the script finds the sequence in the file
##          CHROMOSOME_I.dna. The outputfile is presented in FASTA format and contains the 
##          start and end positions, as specified in the gff format, as part of the name of the 
##          sequence (starting with a "B" and an "E"). 
##
## For documentation on the gff format, see http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
##
## Send comments to mueller@its.caltech.edu.

$source = $ARGV[0];
$feature = $ARGV[1];
$strand = $ARGV[2];
$gff = $ARGV[3];
$fasta = $ARGV[4];
$outfile = $ARGV[5];

if ($strand =~ /[^+-]/) {die "No correct strand specified."};
open (GENES, "$gff") or die "Cannot open $gff.";
open (INFASTA, "$fasta") or die "Cannot open $fasta.";
open (OUTFASTA, ">$outfile") or die "Cannot open $outfile.";

## skip first line
$line1st = <INFASTA>;
$line = <INFASTA>;

$all = "";
while ($line ne "") {
    chop($line);
    $all .= $line;
    $line = <INFASTA>;
}

$line = <GENES>;
## skip first lines of comments
while (index($line,"#") == 0) {$line = <GENES>;}

$line =~  s/[ \t]+/ /g;
while ($line ne "") {
    @words = split(/ /, $line);
    $aux = substr($line1st, 1, length($line1st) - 1);
    chomp($aux);
    if ($words[0] ne $aux) {
	print ("Warning: Sequence names in $gff and $fasta are not the same!\n");
    }
    if (($words[1] eq $source) && ($words[2] eq $feature) && ($words[6] eq $strand)) {
	$chr = substr($all, $words[3] - 1,$words[4] - $words[3]);
	chomp($words[0]);
	print OUTFASTA (">$words[0]");
	print OUTFASTA ("B$words[3]");
	print OUTFASTA ("E$words[4]\n");
	$charactercount = length($chr);
	for ($j = 0; $j < $charactercount; $j +=80) {
	    $chr1 = substr($chr, $j, 80);
	    if ($chr1 ne "") {print OUTFASTA ("$chr1\n");}
	}
    }
    $line = <GENES>;
    $line =~  s/[ \t]+/ /g;
}
close (GENES) or die "Cannot close $gff.";
close (INFASTA) or die "Cannot close $fasta."; 
close (OUTFASTA) or die "cannot close $outfile.";

