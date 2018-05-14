#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $data_ref;
my %types    = ( dauer_entry => 1,
                 dauer_exit  => 1,
                 herm_L4     => 1,
                 L1          => 1,
                 L2          => 1,
                 L3          => 1,
                 male_L4     => 1, 
                 YA          => 1, );

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \/ (\w+) \. [^\/\s]+ \/ [^\/\s]+ \.fastq\.gz /xms ) { 
        my $stage = $1;
        if ( exists $types{$stage} ) { 
             $data_ref->{'stage'}->{$stage}->{'input'}->{$input} = 1;
        }
    }
}

print "\n";
my @stages = sort keys %{ $data_ref->{'stage'} };
foreach my $stage (@stages) {
    my @inputs = sort keys %{ $data_ref->{'stage'}->{$stage}->{'input'} };
    foreach my $input (@inputs) { 
        print "$stage\t$input\n";
    }
    print "\n";
}


