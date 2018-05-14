#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile    = q{};
my $selection = q{};
my $print     = 'batch_only';

$infile    = $ARGV[0] if $ARGV[0];
$selection = $ARGV[1] if $ARGV[1];

if ( ( $selection ne 'genotype_ok' ) and ( $selection ne 'batch_only' ) ) {
    die "Format: split_func_summ_lines_2016.02.23.01.pl [infile] ['genotype_ok' or 'batch_only']\n";
}

if (! -r $infile) {
    warn "Cannot read infile: $infile\n";
    die "Format: split_func_summ_lines_2016.02.23.01.pl [infile] ['genotype_ok' or 'batch_only']\n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A GO_term  /xms ) { 
        print "$input\n";
    }
    elsif ( $input =~ /\A [^\t]+ \t ([^\t]+) /xms ) {
        my $sig_text = $1;
        $print       = 'batch_only';
        while ( $sig_text =~ /\[ ([^\[]+) \] /xmsg ) { 
            my $condition = $1;
            if ( $condition !~ /batch/ ) {
                $print = 'genotype_ok';
            }
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


