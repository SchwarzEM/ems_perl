#!/usr/bin/env perl

# extract_Eval_ratios.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/3/2009.
# Purpose: given two BLAST reports, take all hits in 'blast1' and get their E-value ratios vs. 'blast2'.
# Note: to avoid division-by-zero, assign each '0' E-value a nominal minimum score of 1e-1000.

use strict;
use warnings;
use Getopt::Long;

my $blast1 = q{};
my $blast2 = q{};

my $query         = q{};
my $count_next_E  = 0;
my $positives_ref = ();

GetOptions ( "blast1|b1=s" => \$blast1,
             "blast2|b2=s" => \$blast2, );

if ( ( $blast1 !~ /\S/ ) or ( $blast2 !~ /\S/ ) ) { 
    die "Missing BLAST1 ($blast1) or BLAST2 ($blast2) argument.\n";
}
if ( (! -r $blast1) or (! -r $blast2) ) { 
    die "Can't read BLAST1 $blast1 or BLAST2 $blast2.\n";
}

open my $BLAST1, '<', $blast1 or die "Can't read index BLAST report $blast1: $!";
while (my $input = <$BLAST1>) { 
    record_E_values($input,$blast1,'index');
}
close $BLAST1 or die "Can't close filehandle to index BLAST report $blast2: $!";

open my $BLAST2, '<', $blast2 or die "Can't read comparison BLAST report $blast2: $!";
while (my $input = <$BLAST2>) { 
    record_E_values($input,$blast2,'comparison');
}
close $BLAST2 or die "Can't close filehandle to comparison BLAST report $blast2: $!";

sub record_E_values {
    my ($input, $sourcefile, $type) = @_;
    chomp $input;
    # E.g.: "Query= NODE_113_length_1388_cov_69.486313"
    if ( $input =~ /\A Query= \s+ (\S+) /xms ) {
        $query = $1;
        $count_next_E = 0;
    }
    # E.g.: "XERC_ECOL6  [...]  442   e-125"
    if ( $count_next_E and ( $input =~ / \A .+ \s+ (\S+) \s* \z/xms ) ) { 
        my $E_value = $1;

        # Fix things like 'e-128' to be Perl-valid floating numbers.
        if ( $E_value =~ /\A e/xms) { 
            $E_value = 1 . $E_value;
        }

        # Guard against division-by-zero with comparison E-values.
        # Through trial and error, found the lowest no. Perl handles is ~< 1e-300:
        if ( ($E_value eq '0.0') or ($E_value == 0) or ($E_value < 1e-300) ) { 
            if ($type eq 'comparison') {             
                $E_value = 1e-300;
            }
        }

        $positives_ref->{$query}->{$sourcefile} = $E_value;
        $count_next_E = 0;
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
    my $blast1_value = q{};
    my $blast1_text  = q{};
    my $blast2_value = q{};
    my $blast2_text  = q{};

    if ( exists $positives_ref->{$positive}->{$blast1} ) { 
        $blast1_value = $positives_ref->{$positive}->{$blast1};
        $blast1_value = sprintf("%.2g", $blast1_value);
        $blast1_text = $blast1_value;
    }
    if (! exists $positives_ref->{$positive}->{$blast1} ) {
        $blast1_value = 10;
        $blast1_text = 'none';
    }
    if ( exists $positives_ref->{$positive}->{$blast2} ) { 
        $blast2_value = $positives_ref->{$positive}->{$blast2};
        $blast2_value = sprintf("%.2g", $blast2_value);
        $blast2_text = $blast2_value;
    }
    if (! exists $positives_ref->{$positive}->{$blast2} ) {
        $blast2_value = 10;
        $blast2_text = 'none';
    }
    my $blast1_2_ratio = ($blast1_value / $blast2_value);
    $blast1_2_ratio = sprintf("%.2g", $blast1_2_ratio);
    my $blast2_1_ratio = ($blast2_value / ( $blast1_value + 1e-300 ) );
    $blast2_1_ratio = sprintf("%.2g", $blast2_1_ratio);

    # If for some reason I want full reports, this 'if' is easy to reverse:
    if ( $blast1_text ne 'none') { 
        print $positive, 
              "\t", 
              $blast1_text, 
              "\t", 
              $blast2_text,
              "\t",
              $blast2_1_ratio,
              "\t",
              $blast1_2_ratio, 
              "\n",
              ;
    }
}

