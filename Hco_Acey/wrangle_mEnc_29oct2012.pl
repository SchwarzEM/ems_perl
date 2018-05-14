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

print '#!/bin/bash', "\n\n";

my @stages = sort keys %{ $data_ref->{'stage'} };
foreach my $stage (@stages) {
    my @inputs = sort keys %{ $data_ref->{'stage'}->{$stage}->{'input'} };
    print "    zcat @inputs | quality_trim_fastq.pl -q 64 -u 38 -n -m 38 -i - -o $stage.modencode.38nt_29oct2012.se.fq ;\n";
    print "    fastq2fa_simple.pl $stage.modencode.38nt_29oct2012.se.fq > $stage.modencode.38nt_29oct2012.se.fa ;\n";
    print "    rm $stage.modencode.38nt_29oct2012.se.fq ;\n";
    print "    gzip -9 $stage.modencode.38nt_29oct2012.se.fa ;\n";
    print "\n";
}

