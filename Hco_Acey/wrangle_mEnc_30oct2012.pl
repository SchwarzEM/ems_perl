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

my $dir_to_open     = $ARGV[0];
my @subdirs_to_open = ();

opendir my $DIR, $dir_to_open or die "Can't open handle to directory $dir_to_open;\n";
my @subdirs = readdir $DIR;
foreach my $subdir (@subdirs) { 
    $subdir = $dir_to_open . q{/} . $subdir;
    push @subdirs_to_open, $subdir;
}
closedir $DIR or die "Can't close handle to directory $dir_to_open;\n";

foreach my $subdir (@subdirs_to_open) { 
    if ( $subdir =~ /\A \S+ \/ (\w+) \. [^\/\s]+ \z /xms ) { 
        my $stage = $1;
        if ( exists $types{$stage} ) { 
             my $input_file_set = $subdir . '/*.fastq.gz';
             $data_ref->{'stage'}->{$stage}->{'input'}->{$input_file_set} = 1;
        }
    }
}

print '#!/bin/bash', "\n\n";

my @stages = sort keys %{ $data_ref->{'stage'} };
foreach my $stage (@stages) {
    my @inputs = sort keys %{ $data_ref->{'stage'}->{$stage}->{'input'} };
    print "    zcat @inputs | quality_trim_fastq.pl -q 64 -u 38 -n -m 38 -i - -o $stage.modencode.38nt_30oct2012.se.fq ;\n";
    print "    fastq2fa_simple.pl $stage.modencode.38nt_30oct2012.se.fq > $stage.modencode.38nt_30oct2012.se.fa ;\n";
    print "    rm $stage.modencode.38nt_30oct2012.se.fq ;\n";
    print "    gzip -9 $stage.modencode.38nt_30oct2012.se.fa ;\n";
    print "\n";
}

