#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

while (my $infile = <>) {
    chomp $infile;
    if ( $infile =~ /\A (\S+) \. stats \. raw \. txt \z/xms ) {
        my $stem    = $1;
        my $outfile = "$stem.stats.filt1.txt";
        $outfile    = safename($outfile);

        my @text = ();
        my $print = 0;  

        open my $INFILE, '<', $infile;
        while (my $input = <$INFILE>) {
            chomp $input;
            if ( $input =~ /\A (.+) \t \S+ \t (\S+) \z/xms ) {
                my $info = $1;
                my $qval = $2;
                my $output = "$info\t$qval";
                if ( $qval eq 'q-value' ) {
                    push @text, $output;
                }
                elsif ( ( looks_like_number($qval) ) and ( $qval <= 0.01 ) ) {
                    $print = 1;
                    push @text, $output;
                }
            }
            else {
                die "From infile $infile, cannot parse input: $input\n";
            }
        }
        close $INFILE;

        if ($print) {
            open my $OUTFILE, '>', $outfile;
            foreach my $output (@text) {
                print $OUTFILE "$output\n";
            }
            close $OUTFILE;
        }

    }
    else {
        die "Cannot parse putative file: $infile\n";
    }
}

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

