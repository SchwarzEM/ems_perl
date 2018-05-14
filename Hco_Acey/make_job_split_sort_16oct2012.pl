#!/usr/bin/env perl

use strict;
use warnings;

my $infile    = $ARGV[0];
my $split     = $ARGV[1];
my $type      = $ARGV[2];
my $type_let  = q{};
my @vals      = ();
my @all_jumbs = ();
my @all_pes   = ();
my @all_ses   = ();

if ( (! $type) or ( ( $type ne 'fasta' ) and ( $type ne 'fastq' ) ) ) {
    die "Format: ./make_job_split_sort_16oct2012.pl [infile] [split no.] [fasta/fastq]\n";
}

if ( $type eq 'fasta' ) {
    $type_let = 'f';
}
if ( $type eq 'fastq' ) {
    $type_let = 'q';
}

foreach my $val (1..$split) { 
    $val = sprintf "%02i", $val;
    push @vals, $val;
}

print '#!/bin/bash', "\n\n";

print "    split_fasta.or.q.pl -$type_let -i $infile -s $split ;\n";
print "\n";
print "    e_ping -p done_10split_split_sort_16oct2012 ;\n";
print "\n";

foreach my $val (@vals) {
    my $in_jumb = "$infile.out.$val"; 
    my $out_pe  = "$infile.out.$val.pe.fa";
    my $out_se  = "$infile.out.$val.se.fa";
    print '    timeout -m 125000000 -o timeout.report.16oct2012.', $val, '.txt paired_vs_unp_fastq.or.a.pl --', $type, ' --r1 "#0\/1" --r2 "#0\/2"',
          " -i $in_jumb -p $out_pe -u $out_se ;\n",
          ;
    push @all_jumbs, $in_jumb; 
    push @all_pes, $out_pe;
    push @all_ses, $out_se;
}

print "\n";
print "    rm @all_jumbs ;\n";
print "\n";
print "    e_ping -p done_$split.way_sort_split_sort_16oct2012 ;\n";
print "\n";
print "    cat @all_pes > $infile.all.done1.pe.fa ;\n";
print "    rm @all_pes ;\n";
print "\n";
print "    cat @all_ses > $infile.all.out1.se.fa ;\n";
print "    rm @all_ses ;\n";
print "\n";
print '    timeout -m 125000000 -o timeout.report.16oct2012b.txt paired_vs_unp_fastq.or.a.pl --', $type, ' --r1 "#0\/1" --r2 "#0\/2"',
      " -i $infile.all.out1.se.fa -p $infile.all.done2.pe.fa -u $infile.all.done2.se.fa ;\n",
      ;
print "\n";
print "    e_ping -p done_all_split_sort_16oct2012 ;\n\n";

