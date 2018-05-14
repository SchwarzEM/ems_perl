#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $prefix = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        # For dealing with really wacky 'species' like "Mammalian expression vector pCMV-Script"
        if ( $input =~ /\A > (\S+ .* \s vector \s (p\S+)) \s* \z/xms ) {
            my $header = $1;
            my $prefix = $2 . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # "Eutheria" is totally a species.  Yup.  So are "Anopheles", etc.
        elsif ( $input =~ /\A > (\S+ .* ([A-Z][a-z]{3}) [a-z]*) \s* \z/xms ) {
            my $header = $1;
            my $prefix = $2 . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # For stuff like "Murine leukemia virus" and "Autographa californica MNPV"
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z]) [a-z]+ \s+ ([a-z]{3}) [a-z]* \s+ (\S{3}) \S*) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2; 
            my $prefix_middle = $3;
            my $prefix_end    = $4;
            my $prefix        = $prefix_start . $prefix_middle . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # "Malus x domestica", "Saccharum hybrid cultivar", etc.
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z]) [a-z]+ \s+ ([a-z]) [a-z]* \s+ ([a-z]{3}) [a-z]*) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_middle = $3;
            my $prefix_end    = $4;
            my $prefix        = $prefix_start . $prefix_middle . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # "Arabidopsis lyrata subsp. lyrata"
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z]) [a-z]+ \s+ ([a-z]{3}) [a-z]* \s+ subsp. \s+ \S+) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            my $prefix       = $prefix_start . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # "Orpinomyces sp. OUS1" and "Capitella sp. 1"
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z]) [a-z]+ \s+ sp\. \s+ (\S+)) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            my $prefix       = $prefix_start . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # "Caenorhabditis sp. 11 MAF-2010" and "Fundulus sp. 'Laguna de Labradores'"
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z]) [a-z]+ \s+ sp\. \s+ (\S+) .*) \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            $prefix_end =~ s/\W//g;
            if ( $prefix_end =~ /\a (\w{3}) \w* /xms ) { 
                $prefix_end = $1;
            }
            my $prefix       = $prefix_start . '.sp.' . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # "Ovis sp."
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z]) [a-z]+ \s+ (sp)\.) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            my $prefix       = $prefix_start . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # The standard, well-behaved case (we allow \S in the species name to allow '-' names).
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z]) [a-z]+ \s+ ([a-z]{3}) \S*) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            my $prefix       = $prefix_start . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # Darn you, "Glycine max" and "Trichoplusia ni"!
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z][a-z]{2}) [a-z]* \s+ ([a-z]{1,3}) ) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            my $prefix       = $prefix_start . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }   
        # "Zea".
        elsif ( $input =~ /\A > (\S+ .* ([A-Z][a-z]{2}) [a-z]*) \s* \z/xms ) {
            my $header = $1;
            my $prefix = $2 . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # unidentified sea urchin"
        elsif ( $input =~ /\A > (\S+ \s .* unidentified \s+ ([a-z]{3}) [a-z]* \s ([a-z]{3}) [a-z]* ) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            my $prefix       = 'unid_' . $prefix_start . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # "Xenopus (Silurana) tropicalis"
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z]) [a-z]+ \s+ \( [A-Z] [a-z]+ \) \s+ ([a-z]{3}) \S*) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            my $prefix       = $prefix_start . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        # "Avian endogenous retrovirus EAV-HP"
        elsif ( $input =~ /\A > (\S+ .* \s+ ([A-Z][a-z]{3}) [a-z]* \s+ .* \s+ retrovirus \s+ (\S+) ) \s* \z/xms ) {
            my $header        = $1;
            my $prefix_start  = $2;
            my $prefix_end    = $3;
            my $prefix       = $prefix_start . '.RV.' . $prefix_end . '.RB_';
            print '>', $prefix, "$header\n";
        }
        else {
            die "Can't parse FASTA header: $input\n";
        }
    }
    else { 
        if ( $input =~ /\S/xms ) { 
            print "$input\n";
        }
    }
}

