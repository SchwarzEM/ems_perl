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
        if ( $input =~ / (?: [^\t]* \t){2} ([^\t]*) \t (?: [^\t]* \t){5} ([^\t]*) \t (?: [^\t]* \t){5} ([^\t]*) \t /xms ) { 
            my $edgeR   = $1;
            my $phobius = $2;
            my $go_term = $3;

            # Sanity check:
            if ( ( $edgeR =~ /\S/xms ) and ( $edgeR !~ /\S+ [ ] to [ ] \S+ [ ] \[ (\+|-): [ ]/xms ) ) {
                die "Can't parse putative edgeR value \"$edgeR\" in input: $input\n";
            }
            if ( ( $go_term =~ /\S/xms ) and ( $go_term !~ /GO:\d+/xms ) ) { 
                die "Can't parse putative GO term value \"$go_term\" in input: $input\n";
            }

            if (     ( ( $input !~ /ENSG\d+\(human\)/xms ) and ( $input !~ /ENSMUSG\d+\(mouse\)/xms ) )
                 and ( ( $edgeR =~ /ACEY\.L3i [ ] to [ ] ACEY\.24\.PI [ ] \[\+: [ ]/xms ) or ( $edgeR =~ /ACEY\.24\.PI [ ] to [ ] ACEY\.5\.D [ ] \[\+: [ ]/xms ) )
                 and ( $edgeR !~ / \S+ [ ] to [ ] \S+ [ ] \[\-/xms ) 
                 and ( $phobius eq 'SigP' )
                 and ( $go_term =~ /(GO:0004867)/xms )    # for serine-type endopeptidase inhibitor activity == GO:0004867
               ) {
                print "$input\n";
            }
        }
    }
}
