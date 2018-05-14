#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};
my $type   = q{};

$infile = $ARGV[0] if $ARGV[0];
$type   = $ARGV[1] if $ARGV[1];

if ( (! -r $infile ) or ( ( $type ne '5spp' ) and ( $type ne '7spp' ) ) ) {
    die "Format: filter_singletons_18sep2016.pl [readable nigoni data infile] [either '5spp' or '7spp']\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    my @orthos = ();

    if ( $input =~ /\A Gene \t /xms ) { 
        print "$input\n";
    }

    elsif ( $input =~ /\A (?: [^\t]* \t){7} 
                          ([^\t]*) \t [^\t]* \t 
                          ([^\t]*) \t [^\t]* \t 
                          ([^\t]*) \t [^\t]* \t 
                          ([^\t]*) \t /xms ) {
        my $ofind_4spp = $1;
        my $ofind_5spp = $2;
        my $ofind_6spp = $3;
        my $ofind_8spp = $4;

        my $print = 'no';
        if ( $type eq '5spp' ) {
            push @orthos, $ofind_4spp, $ofind_5spp, $ofind_6spp;
        }
        elsif ( $type eq '7spp' ) {
            push @orthos, $ofind_4spp, $ofind_5spp, $ofind_6spp, $ofind_8spp;
        }        
        foreach my $ortho (@orthos) {
            # Instance:
            # OG0041267(1 genes,1 taxa):   Cnig_chr_I.g1(nigoni)
            if ( $ortho =~ /\A OG\d+ \( \d+ \s+ genes, (\d+) \s+ taxa \) [:] \s+ /xms ) {
                my $taxon_count = $1;
                if ( $taxon_count >= 2 ) {
                    $print = 'yes';
                }
            }
            else {
                die "Cannot parse OrthoFinder group \"$ortho\" in line: $input\n";
            }
        }
        if ( $print eq 'yes' ) {
            print "$input\n";
        }
    }
    else {
        die "Cannot parse input line: $input\n";
    }
}
close $INFILE;

