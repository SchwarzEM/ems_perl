#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A Gene /xms ) { 
        $input =~ s/Phobius/Secreted/;
        $input =~ s/ASPR_36align/36_aligned/;
    }
    elsif (  $input =~ /\A Acey_2012.08.05_ /xms ) {
        $input =~ s/ORTHOMCL896.14spp/+/;
        $input =~ s/ASPR_36align/+/;
        if ( $input =~ /\A ([^\t]* \t [^\t]* \t [^\t]* \t) ([^\t]*) (\t [^\t]* \t [^\t]* \t [^\t]*) /xms ) {
            my $text1 = $1;
            my $text2 = $2;
            my $text3 = $3;
            if ( $text2 ne 'SigP' ) {
                $text2 = q{};
            }
            if ( $text2 eq 'SigP' ) { 
                $text2 = q{+};
            }
            $input = $text1 . $text2 . $text3;
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }
    print "$input\n";
}

    
