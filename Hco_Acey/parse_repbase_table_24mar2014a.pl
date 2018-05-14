#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $name    = q{};
my $type    = q{};
my $species = q{};
my $source  = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (\S+ \.RB_ \S+) \t ([^\t]+) \t ([^\t]+) \z/xms ) {
            $name    = $1;
            $type    = $2;
            $species = $3;
            $source  = 'RepBase';
            $type    =~ s/\A\s+//;
            $type    =~ s/\s+\z//;
            $species =~ s/\A\s+//;
            $species =~ s/\s+\z//;
        }
        # Acey_RScout_0001
        elsif ( $input =~ /\A > (Acey_RScout_\d+) \s /xms ) { 
            $name    = $1;
            $type    = 'ab initio prediction';
            $species = 'Ancylostoma ceylanicum';
            $source  = 'This study';
        }
        elsif ( $input =~ /\A > (Cele_R=\d+) \s /xms ) {
            $name    = $1;
            $type    = 'ab initio prediction';   
            $species = 'Caenorhabditis elegans';
            $source  = 'This study';
        }
        elsif ( $input =~ /\A > (Cbri_R=\d+) \s /xms ) {
            $name    = $1;
            $type    = 'ab initio prediction';
            $species = 'Caenorhabditis briggsae';
            $source  = 'This study';
        }
        elsif ( $input =~ /\A > (Nam_R=\d+) \s /xms ) {
            $name    = $1;
            $type    = 'ab initio prediction';
            $species = 'Necator americanus';
            $source  = 'This study';
        }
        else { 
            die "Cannot parse sequence line: $input\n";
        }
        print "$name\t$species\t$type\t$source\n";
    }
}

