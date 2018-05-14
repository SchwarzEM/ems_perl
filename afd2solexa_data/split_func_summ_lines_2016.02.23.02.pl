#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile    = q{};
my $selection = q{};
my $print     = 'batch_best';

$infile    = $ARGV[0] if $ARGV[0];
$selection = $ARGV[1] if $ARGV[1];

if ( ( $selection ne 'genotype_best' ) and ( $selection ne 'batch_best' ) ) {
    die "Format: split_func_summ_lines_2016.02.23.02.pl [infile] ['genotype_best' or 'batch_best']\n";
}

if (! -r $infile) {
    warn "Cannot read infile: $infile\n";
    die "Format: split_func_summ_lines_2016.02.23.01.pl [infile] ['genotype_best' or 'batch_best']\n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A GO_term /xms ) { 
        print "$input\n";
    }
    elsif ( $input =~ /\A [^\t]+ \t ([^\t]+) /xms ) {
        my $sig_text = $1;
        $print       = 'batch_best';
        if ( $sig_text =~ /\A [^\]]+ \[ ( [^\[]+ ) \] /xms ) { 
            my $condition = $1;
            if ( $condition !~ /batch/ ) {
                $print = 'genotype_best';
            }
        }
        else {
            die "Cannot parse significance field of input line: $input\n";
        }
        if ( $selection eq $print ) {
             print "$input\n";
        }
    }
    else {
        die "Cannot parse input line: $input\n";
    }
} 
close $INFILE;

