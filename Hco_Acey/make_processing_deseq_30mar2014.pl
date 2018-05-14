#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $report = 'DESeq_summary_30mar2014';

my $header_not_printed = 1;

my $type = q{};

while (my $input = <>) {
    chomp $input;
    # [Sample inputs:]
    # DESeq.comp.ACEY.L3i.vs.ACEY.24HCM.30mar2013.raw.txt
    # [or]
    # DESeq.comp.ACEY.L3i.vs.ACEY.24HCM.30mar2013.post_12.D_trip.raw.txt
    if ( $input =~ /\A (DESeq\.comp\. (\S+) \. vs \. (\S+) \.30mar2013 (.*)) \.raw\.txt \z/xms ) { 
        my $stem   = $1;
        my $cond_a = $2;
        my $cond_b = $3;
        $type      = $4;

        if ( $type =~ /trip/xms ) { 
            $type = 'post_12.D';
        }
        else {
            $type = 'post_17.D';
        }
        my $header = '#!/bin/bash' . "\n\n" . "    rm $report.$type.txt ;\n" . "    touch $report.$type.txt ;\n\n";
        print $header if $header_not_printed;
        $header_not_printed = 0;
        my $output = "$stem.filt.txt";

        print "    cat $input", ' | grep Acey | perl -ne \' s/["]//g; s/[ ]+/\t/g; print \' | cut -f 2,7,9 > ', "$output ;\n";

        print "    echo \"$cond_a to $cond_b: upregulated, q <= 0.001 [analysis type $type]\" >> $report.$type.txt ;\n";
        print "    cat $output",
              ' | perl -ne \' $input = $_; chomp $input; if ( $input =~ /\A[^t]+\t(\S+)\t(\S+)\z/xms )',
              ' { $log_reg = $1; $q_val = $2; if ( ( $q_val ne \'NA\' ) and ( $log_reg > 0 ) and ( $q_val <= 0.001 ) )',
              ' { print "$input\n"; } } \' | wc -l',
              " >> $report.$type.txt ;\n",
              ;
        print "    echo \"$cond_a to $cond_b: downregulated, q <= 0.001 [analysis type $type]\" >> $report.$type.txt ;\n";
        print "    cat $output",
              ' | perl -ne \' $input = $_; chomp $input; if ( $input =~ /\A[^t]+\t(\S+)\t(\S+)\z/xms )',
              ' { $log_reg = $1; $q_val = $2; if ( ( $q_val ne \'NA\' ) and ( $log_reg < 0 ) and ( $q_val <= 0.001 ) )',
              ' { print "$input\n"; } } \' | wc -l',
              " >> $report.$type.txt ;\n",
              ;
        print "    echo >> $report.$type.txt ;\n";
        print "\n";
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

print "    e_ping -p done.DESeq.$report.$type ;\n\n";


