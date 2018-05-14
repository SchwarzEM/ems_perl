#!/usr/bin/env perl

use strict;
use warnings;

my $header = q{};

while (my $input = <>) {
    chomp $input;
    if (! $header) { 
        print "$input\n";
        $header = $input;
    }
    else { 
        # edgeR = 3rd column; Phobius = 9th column; GO terms = 15th column
        if ( $input =~ / (?: [^\t]* \t){2} ([^\t]*) \t (?: [^\t]* \t){5} ([^\t]*) \t /xms ) { 
            my $edgeR   = $1;
            my $phobius = $2;

            # Sanity check:
            if ( ( $edgeR =~ /\S/xms ) and ( $edgeR !~ /\S+ [ ] to [ ] \S+ [ ] \[ (\+|-): [ ]/xms ) ) {
                die "Can't parse putative edgeR value \"$edgeR\" in input: $input\n";
            }

            # Accept any gene whose products:
            #     are unambiguously secreted
            #     are conserved in Haemonchus, but absent from mouse or human
            #     are significant upreg. in 24.PI and/or 5.D
            #     and have no sig. downreg. anywhere

            if (     ( ( $input =~ /\S+\(haemonchus/xms ) and ( $input !~ /ENSG\d+\(human\)/xms ) and ( $input !~ /ENSMUSG\d+\(mouse\)/xms ) )
                 and ( ( $edgeR =~ /ACEY\.L3i [ ] to [ ] ACEY\.24\.PI [ ] \[\+: [ ]/xms ) or ( $edgeR =~ /ACEY\.24\.PI [ ] to [ ] ACEY\.5\.D [ ] \[\+: [ ]/xms ) )
                 and ( $edgeR !~ / \S+ [ ] to [ ] \S+ [ ] \[\-/xms ) 
                 and ( $phobius eq 'SigP' ) 
               ) {
                print "$input\n";
            }
        }
    }
}

