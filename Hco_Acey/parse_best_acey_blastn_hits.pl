#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    my $acey_gene = q{};
    my $rep_ele   = q{};
    my $e_value   = q{};
    if ( $input =~ /\A (\S+)\.t\d+ \s+ (\S+) \s+ (?: \S+ \s+){8} (\S+) /xms ) { 
        $acey_gene = $1;
        $rep_ele   = $2;
        $e_value   = $3;
        if ( exists $data_ref->{'acey_gene'}->{$acey_gene}->{'e_value'} ) { 
            my $past_e_value = $data_ref->{'acey_gene'}->{$acey_gene}->{'e_value'};
            if ( $e_value < $past_e_value ) { 
                $data_ref->{'acey_gene'}->{$acey_gene}->{'e_value'} = $e_value;
                $data_ref->{'acey_gene'}->{$acey_gene}->{'rep_ele'}  = $rep_ele;
            }
        }
        else { 
            $data_ref->{'acey_gene'}->{$acey_gene}->{'e_value'} = $e_value;
            $data_ref->{'acey_gene'}->{$acey_gene}->{'rep_ele'} = $rep_ele;
        }
    }
}

my $header = "Gene\tBest_Repeat_BlastN [E-value]";

my @genes = sort keys %{ $data_ref->{'acey_gene'} };
foreach my $gene (@genes) { 
    if ($header) { 
        print "$header\n";
        $header = q{};
    }
    my $rep_ele = $data_ref->{'acey_gene'}->{$gene}->{'rep_ele'};
    my $e_value = $data_ref->{'acey_gene'}->{$gene}->{'e_value'};
    print "$gene\t$rep_ele [$e_value]\n";
}

