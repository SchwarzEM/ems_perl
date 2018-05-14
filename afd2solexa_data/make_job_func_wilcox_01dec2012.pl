#!/usr/bin/env perl

use strict;
use warnings;

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \.func_input\.txt \z/xms ) { 
        my $stem       = $1;
        my $output_dir = $stem . '.func_wilcoxon_outdir';
        $output_dir    = safename($output_dir);

        # Print header -- only once, and only if there's something else to print:
        print $header if $header;
        $header = q{};

        # Then print one stanza per $input:
        print "    mkdir $output_dir ;\n";
        print "    nohup func_wilcoxon -t /home/schwarzc/func_work/go_201011-termdb-tables -i $input -o $output_dir &\n";
        print "\n";
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

