#!/usr/bin/env perl

# deseq_csv2tsv.pl -- Erich Schwarz <ems394@cornell.edu>, 7/29/2015.
# Purpose: convert typical CSV from DESeq2 to less typical and more usable *.tsv.txt.

use strict;
use warnings;
use autodie;

my $infile = q{};
my $tag    = q{};

$infile = $ARGV[0] if $ARGV[0];
$tag    = $ARGV[1] if $ARGV[1];

if ( (! -r $infile ) or ( $tag !~ /\S/xms ) ) { 
    die "Format: deseq_csv2tsv.pl [readable infile] [tag with at least one non-space] (print to STDOUT or >outfile)\n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;

    if ( $input !~ /\A (?: [^,]+ [,]){6} [^,]+ \z/xms ) { 
        die "From input file $infile, cannot parse commas in line: $input\n";
    }
    $input =~ s/[,]/\t/g;

    if ( $input =~ /\A \"\" \t \"baseMean\" \t \"log2FoldChange\" \t \"lfcSE\" \t \"stat\" \t \"pvalue\" \t "padj\"/xms ) {
        $input =~ s/\A\"\"/Gene/;
        $input =~ s/([^\"\s])\"/$1$tag/g;
        $input =~ s/\"//g;
    }
    elsif ( $input =~ /\A \" \S+ \" \t /xms ) {
        $input =~ s/\"//g;
    }
    else {
        die "From input file $infile, cannot correctly parse format of text line: $input\n";
    }
    print "$input\n";
}
close $INFILE;

