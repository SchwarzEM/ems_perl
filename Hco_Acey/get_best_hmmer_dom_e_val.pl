#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

while (my $input = <>) { 
    chomp $input;

    # Sample of wanted input text:
    # Acey_s0065.g3635.t1  -          TF352155.hmmer_aln.PF00059.16 -              4e-26   92.4  10.5   1.5e-17   64.8   0.5   2.7   3   0   0   3   3   2   2 -
    # This:    ***                                                                            and this: ******  

    if ( ( $input !~ /\A [#]/xms ) and ( $input =~ /\A (\S+) (?: \s+ \S+){6} \s+ (\S+) /xms ) ) {
        my $query               = $1;
        my $best_domain_e_value = $2;
        if ( exists $data_ref->{'query'}->{$query} ) { 
            die "Redundant data line with second instance of query $query: $input\n";
        }
        $data_ref->{'query'}->{$query}->{best_domain_e_value} = $best_domain_e_value;
    }
}

my @queries = sort { $data_ref->{'query'}->{$a}->{best_domain_e_value} <=> $data_ref->{'query'}->{$b}->{best_domain_e_value} } 
                  keys %{ $data_ref->{'query'} };

foreach my $query (@queries) { 
    my $best_domain_e_value = $data_ref->{'query'}->{$query}->{best_domain_e_value};
    print "$query\t$best_domain_e_value\n";
}


