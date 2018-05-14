#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $data_ref;

my @ok_types = ( '24.PI/L3i',
                 'L3i/24.PI',                

                 '24HCM/L3i',
                 'L3i/24HCM',

                 '24.PI/24HCM',
                 '24HCM/24.PI',

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
    if ( $input =~ / log10 \( (\S+) \) \t ([^\t]+) \t (\S+) \t (\S+) /xms ) { 
        my $comparison      = $1;
        my $protein_feature = $2;
        my $log_expr_dist   = $3;
        my $p_value         = $4;
        if ( ( looks_like_number($log_expr_dist) and looks_like_number($p_value) ) and ( $log_expr_dist > 0 ) ) {
            if (! exists $ok_match{$comparison} ) {
                die "Do not recognize the comparison in: $input\n";
            }
            my $annot = "$protein_feature\t$p_value\n";
            push @{ $data_ref->{'type'}->{$comparison} }, $annot;
        }
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

foreach my $ok_type (@ok_types) {
    if ( exists $data_ref->{'type'}->{$ok_type} ) { 
        print "Significant protein features, among genes upregulated in $ok_type:\t\n";
        print @{ $data_ref->{'type'}->{$ok_type} };
    }
    else { 
        warn "Can't find annotations for genes upregulated in $ok_type!\n";
    }
}

