#!/usr/bin/env perl

# meme_pfm2naive_pfm.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/22/2010.
# Purpose: convert a meme PFM to a naive PFM, web-pastable into Tomtom.

use strict;
use warnings;

my $data_ref;

my $typecount    = 0;
my $orig_length  = 0;
my $length       = 0;
my $motif_number = 0;

my $data_line    = q{};
my @values       = ();

while (my $input = <>) { 
    chomp $input;
    # letter-probability matrix: alength= 4 w= 15 nsites= 50 E= 9.6e-020
    if ( $input =~ / \A letter-probability \s+ matrix: \s+ alength= \s+ (\d+) \s+ w= \s+ (\d+) /xms ) { 
        $typecount = $1;
        $length    = $2;
        $motif_number++;
        $orig_length = $length + 1;
        $typecount--;
    }
    #  0.600000  0.140000  0.120000  0.140000
    elsif ( ( $length > 0 ) and ( $input =~ / \A \s* ( \S+ (?: \s+ \S+){$typecount} ) \s* /xms ) ) { 
        $data_line = $1;
        @values = split /\s+/, $data_line;
        foreach my $value (@values) { 
            my $residue = $orig_length - $length;
            push @{ $data_ref->{$motif_number}->{$residue} }, $value;
        }
        $length--;
    }
    elsif ( $input =~ /\S/xms ) { 
        die "Can't parse: $input\n";
    }
}

$typecount++;
foreach my $motif ( sort { $a <=> $b } keys %{ $data_ref } ) { 
    foreach my $i (1..$typecount) { 
        $i--;
        my @output_values = ();
        foreach my $res ( sort { $a <=> $b} keys %{ $data_ref->{$motif} } ) { 
            @values = @{ $data_ref->{$motif}->{$res} };
            push @output_values, $values[$i];
        }
        my $output_line = join "  ", @output_values;
        print "$output_line\n";
    }
}

