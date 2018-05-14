#!/usr/bin/env perl

# tabulate_motif_statistics_23apr2013.pl -- Erich Schwarz <ems394@cornell.edu>, 4/24/2013.
# Purpose: given lots of data, export it into namesafed files, one file per data set.  Goal is to avoid overloading Word.

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
        my $ok_label = $ok_type;
        $ok_label =~ s/\//.vs./;
        my $output = 'Acey_2012.10.24.' . $ok_label . '.motifs_summary_24apr2013.txt';
        $output = safename($output);
        open my $OUTPUT, '>', $output or die "Can't open output file: $output\n";
        print $OUTPUT "Significant protein features, among genes upregulated in $ok_type:\t\n";
        print $OUTPUT @{ $data_ref->{'type'}->{$ok_type} };
        close $OUTPUT or die "Can't close filehandle to output file: $output\n";
    }
    else { 
        warn "Can't find annotations for genes upregulated in $ok_type!\n";
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


