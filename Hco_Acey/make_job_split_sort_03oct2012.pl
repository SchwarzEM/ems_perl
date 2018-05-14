#!/usr/bin/env perl

use strict;
use warnings;

my $infile  = $ARGV[0];
my $split   = $ARGV[1];
my @vals    = ();
my @all_pes = ();
my @all_ses = ();

foreach my $val (1..$split) { 
    $val = sprintf "%02i", $val;
    push @vals, $val;
}

print '#!/bin/bash', "\n\n";

print "    split_fasta.or.q.pl -f -i $infile -s $split ;\n";
print "\n";
print "    program_done_e-ping.pl -p done_10split_split_sort_03oct2012.pl ;\n";
print "\n";

foreach my $val (@vals) { 
    my $out_pe = "$infile.out.$val.pe.fa";
    my $out_se = "$infile.out.$val.se.fa";
    print '    paired_vs_unp_fastq.or.a.pl --fasta --r1 "#0\/1" --r2 "#0\/2"',
          " -i $infile.out.$val -p $out_pe -u $out_se ;\n",
          ;
    push @all_pes, $out_pe;
    push @all_ses, $out_se;
}

print "\n";
print "    program_done_e-ping.pl -p done_10sort_split_sort_03oct2012.pl ;\n";
print "\n";
print "    cat @all_pes > $infile.all.out.pe.fa ;\n";
print "    cat @all_ses > $infile.all.out.se.fa ;\n";
print '    paired_vs_unp_fastq.or.a.pl --fasta --r1 "#0\/1" --r2 "#0\/2"',
      " split -i $infile.all.out.se.fa -p $infile.all.out2.pe.fa -u $infile.all.out2.se.fa ;\n",
      ;
print "\n";
print "    program_done_e-ping.pl -p done_all_split_sort_03oct2012.pl ;\n\n";
