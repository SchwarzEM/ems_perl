#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $data_ref;

my @ok_types = ( '24.PI/L3i',
                 'L3i/24.PI',                

                 '24.HCM/L3i',
                 'L3i/24.HCM',

                 '24.PI/24.HCM',
                 '24.HCM/24.PI',

                 '5.D/24.PI',
                 '24.PI/5.D',

                 '5.D/L3i',   
                 'L3i/5.D',

                 '12.D/5.D',
                 '5.D/12.D',

                 '17.D/12.D',
                 '12.D/17.D',

                 '19.D/17.D',
                 '17.D/19.D',

                 '18D.Alb.4hr/18D.noAlb.4hr',
                 '18D.noAlb.4hr/18D.Alb.4hr',

                 '18D.Cry5B.4hr/18D.HEPES.4hr',
                 '18D.HEPES.4hr/18D.Cry5B.4hr',

                 '18D.Cry5B.24hr/18D.HEPES.24hr',
                 '18D.HEPES.24hr/18D.Cry5B.24hr',

                 '18D.SB.plusCry5B/18D.SB.plusHEPES',
                 '18D.SB.plusHEPES/18D.SB.plusCry5B',
);

my %ok_match = ();
foreach my $ok_type (@ok_types) {
    $ok_match{$ok_type} = 1;
}

while (my $input = <>) {
    chomp $input;
    if ( $input =~ / log10 \( (\S+) \) \t ([^\t]+) \t (\d+) \t (\S+) \t (\S+) \t (\S+) /xms ) { 
        my $comparison      = $1;
        my $protein_feature = $2;
        my $gene_count      = $3;
        my $log_expr_dist   = $4;
        my $p_value         = $5;
        my $q_value         = $6;

        # $gene_count is *forced* to be \d+, i.e., a non-negative integer; so we do not need to enforce that further in the next line.

        if ( ( looks_like_number($log_expr_dist) and looks_like_number($p_value) and looks_like_number($q_value) ) and ( $log_expr_dist > 0 ) ) {
            if (! exists $ok_match{$comparison} ) {
                die "Do not recognize the comparison in: $input\n";
            }
            my $annot = "$protein_feature\t$gene_count\t$p_value\t$q_value\n";
            if ( $q_value <= 0.05 ) { 
                push @{ $data_ref->{'type'}->{$comparison} }, $annot;
            }
        }
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

foreach my $ok_type (@ok_types) {
    if ( exists $data_ref->{'type'}->{$ok_type} ) { 
        print "Significant protein features, among genes upregulated in $ok_type\tGenes with feature\tp-value\tq-value\n";
        print @{ $data_ref->{'type'}->{$ok_type} };
        print "\t\n";
    }
    else { 
        warn "Can't find annotations for genes upregulated in $ok_type!\n";
    }
}

