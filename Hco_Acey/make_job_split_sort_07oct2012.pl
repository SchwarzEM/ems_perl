#!/usr/bin/env perl

use strict;
use warnings;

my $infile    = $ARGV[0];
my $split     = $ARGV[1];
my @vals      = ();
my @all_jumbs = ();
my @all_pes   = ();
my @all_ses   = ();

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
    my $in_jumb = "$infile.out.$val"; 
    my $out_pe  = "$infile.out.$val.pe.fa";
    my $out_se  = "$infile.out.$val.se.fa";
    print '    timeout -m 125000000 -o timeout.report.07oct2012.', $val, '.txt paired_vs_unp_fastq.or.a.pl --fasta --r1 "#0\/1" --r2 "#0\/2"',
          " -i $in_jumb -p $out_pe -u $out_se ;\n",
          ;
    push @all_jumbs, $in_jumb; 
    push @all_pes, $out_pe;
    push @all_ses, $out_se;
}

print "\n";
print "    rm @all_jumbs ;\n";
print "\n";
print "    program_done_e-ping.pl -p done_10sort_split_sort_03oct2012.pl ;\n";
print "\n";
print "    cat @all_pes > $infile.all.done1.pe.fa ;\n";
print "    rm @all_pes ;\n";
print "\n";
print "    cat @all_ses > $infile.all.out1.se.fa ;\n";
print "    rm @all_ses ;\n";
print "\n";
print '    timeout -m 125000000 -o timeout.report.07oct2012b.txt paired_vs_unp_fastq.or.a.pl --fasta --r1 "#0\/1" --r2 "#0\/2"',
      " -i $infile.all.out1.se.fa -p $infile.all.done2.pe.fa -u $infile.all.done2.se.fa ;\n",
      ;
print "\n";
print "    program_done_e-ping.pl -p done_all_split_sort_07oct2012.pl ;\n\n";

