#!/usr/bin/env perl

# comp_2_GOlists.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/10/2010.
# Purpose: given two or more TSVs from refined GO term lists of FUNC, find terms shared by all of them, and assign them their least significant (largest) observed p-value.

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $term      = q{};
my $p_value   = q{};
my $seen_ref;
my @infiles   = @ARGV;
my $termcount = @infiles;

foreach my $infile (@infiles) { 
    open my $INFILE, '<', $infile or die "Can't open input file $infile: $!";
    while ( my $input = <$INFILE> ) { 
        chomp $input;
        if ( $input =~ /\A ([^\t]+ \t [^\t]+ \t GO:\d+) .* \t ([^\t]+) /xms ) { 
            $term    = $1;
            $p_value = $2;
            if (! looks_like_number($p_value) ) { 
                die "Final term of this text is not a number: $input\n";
            }
            $seen_ref->{$term}->{'infile'}->{$infile} = 1;
            if ( (! exists $seen_ref->{$term}->{'p_value'}                )
                 or (     ( exists $seen_ref->{$term}->{'p_value'}      ) 
                      and ( $p_value >= $seen_ref->{$term}->{'p_value'} ) ) ) { 
                $seen_ref->{$term}->{'p_value'} = $p_value;
            }
        }
    }
    close $INFILE or die "Can't close filehandle to input file $infile: $!";
}

my @termlist = sort keys %{ $seen_ref };
foreach my $term1 (@termlist) { 
    my $sightings = keys %{ $seen_ref->{$term1}->{'infile'} };
    if ( $sightings == $termcount ) { 
        my $p_val1 = $seen_ref->{$term1}->{'p_value'};
        print "$term1\t$p_val1\n";
    }
}

