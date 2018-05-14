#!/usr/bin/env perl

# extract_Evals_pctIs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/3/2009.
# Purpose: given a BLAST report, for each hit, get an E-value and the first 'Identity = ' line.

use strict;
use warnings;
use Getopt::Long;

my $blast1 = q{};

my $query           = q{};
my $count_next_E    = 0;
my $count_next_pctI = 0;
my $positives_ref   = ();

GetOptions ( "blast1|b1=s" => \$blast1, );

if ( $blast1 !~ /\S/ ) { 
    die "Missing BLAST1 argument.\n";
}
if (! -r $blast1) { 
    die "Can't read BLAST1 $blast1.\n";
}

open my $BLAST1, '<', $blast1 or die "Can't read index BLAST report $blast1: $!";
while (my $input = <$BLAST1>) { 
    record_E_vals_pctIs($input);
}
close $BLAST1 or die "Can't close filehandle to index BLAST report $blast1: $!";

sub record_E_vals_pctIs {
    my $input = $_[0];
    chomp $input;

    # E.g.: "Query= NODE_113_length_1388_cov_69.486313"
    if ( $input =~ /\A Query= \s+ (\S+) /xms ) {
        $query           = $1;
        $count_next_E    = 0;
        $count_next_pctI = 0;
    }

    # E.g.: Identities = 75/78 (96%)
    if ( $count_next_pctI and ( $input =~ / ( \A \s+ Identities \s+ = \s+ \d+\/\d+ .+ ) \z /xms ) ) { 
        my $pct_Id = $1;
        $positives_ref->{$query}->{'percent_identity'} = $pct_Id;
        $count_next_pctI = 0;
    }

    # E.g.: "XERC_ECOL6  [...]  442   e-125"
    if ( $count_next_E and ( $input =~ / \A .+ \s+ (\S+) \s* \z/xms ) ) { 
        my $E_value = $1;

        # Fix things like 'e-128' to be Perl-valid floating numbers.
        if ( $E_value =~ /\A e/xms) { 
            $E_value = 1 . $E_value;
        }
        $positives_ref->{$query}->{'E-value'} = $E_value;
        $count_next_E    = 0;
        $count_next_pctI = 1;
    }

    # E.g.: "Sequences producing significant alignments: [...] (bits) Value"
    if ( $input =~ / \A 
                     Sequences\s
                     producing\s
                     significant\s
                     alignments:\s+
                     \(bits\)\s
                     Value\s* \z 
                   /xms ) {
        $count_next_E = 1;
    }
}

my @positives = sort keys %{ $positives_ref };
foreach my $positive (@positives) { 
    my $E_value = q{};
    my $pct_Id  = q{};

    if ( exists $positives_ref->{$positive}->{'E-value'} ) { 
        $E_value = $positives_ref->{$positive}->{'E-value'};
        $E_value = sprintf("%.2g", $E_value);
    }
    if ( exists $positives_ref->{$positive}->{'percent_identity'} ) {
        $pct_Id = $positives_ref->{$positive}->{'percent_identity'};
    }
    print $positive, 
        "\t", 
        $E_value, 
        "\t", 
        $pct_Id,
        "\n",
        ;
}

