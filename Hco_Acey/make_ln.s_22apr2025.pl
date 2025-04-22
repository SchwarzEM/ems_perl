#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %rename = (
    'HW1'  => 'HWDMSO_rep1',
    'HW2'  => 'HWDMSO_rep2',
    'HW4'  => 'HWDMSO_rep3',
    'HW5'  => 'HWDMSO_rep4',
    'HW6'  => 'HWEMO_8n_rep1',
    'HW7'  => 'HWEMO_8n_rep2',
    'HW8'  => 'HWEMO_8n_rep3',
    'HW9'  => 'HWEMO_8n_rep4',
    'HW10' => 'HWEMO_8n_rep5',
    'HW11' => 'HWA9_3_rep1',
    'HW12' => 'HWA9_3_rep2',
    'HW13' => 'HWA9_3_rep3',
    'HW14' => 'HWA9_3_rep4',
    'HW15' => 'HWA9_3_rep5',
    'HW16' => 'HWAG_30_rep1',
    'HW17' => 'HWAG_30_rep2',
    'HW18' => 'HWAG_30_rep3',
    'HW20' => 'HWAG_30_rep4',
    'HW21' => 'HWAG_20_rep1',
    'HW22' => 'HWAG_20_rep2',
    'HW23' => 'HWAG_20_rep3',
    'HW24' => 'HWAG_20_rep4',
    'HW25' => 'HWAG_20_rep5',
    'WW1'  => 'WWDMSO_rep1',
    'WW2'  => 'WWDMSO_rep2',
    'WW3'  => 'WWDMSO_rep3',
    'WW4'  => 'WWDMSO_rep4',
    'WW5'  => 'WWDMSO_rep5',
    'WW6'  => 'WWEMO_75n_rep1',
    'WW7'  => 'WWEMO_75n_rep2',
    'WW8'  => 'WWEMO_75n_rep3',
    'WW9'  => 'WWEMO_75n_rep4',
    'WW10' => 'WWEMO_75n_rep5',
    'WW11' => 'WWA9_10_rep1',
    'WW12' => 'WWA9_10_rep2',
    'WW13' => 'WWA9_10_rep3',
    'WW14' => 'WWA9_10_rep4',
    'WW15' => 'WWA9_10_rep5',
    'WW16' => 'WWAG_30_rep1',
    'WW17' => 'WWAG_30_rep2',
    'WW18' => 'WWAG_30_rep3',
    'WW19' => 'WWAG_30_rep4',
    'WW20' => 'WWAG_30_rep5',
);

while ( my $input = <> ) {
    chomp $input;
    # Sample input: 01.HWRawData/HW10_1.fq.gz
    if ( $input =~ / (\S+ \/) ([A-Z]{2}\d+) (_\d\.fq\.gz)  /xms ) {
        my $lead_text = $1;
        my $old_tag   = $2;
        my $back_text = $3;

        my $new_tag = q{};
        my $output  = q{};

        if ( exists $rename{$old_tag} ) {
            $new_tag = $rename{$old_tag};
            $output = $lead_text . $new_tag . $back_text;
        }

        if ($output) {
            print "    ln -s $input $output ;\n";
        }
        else {
            warn "Could not rename tag in: $input\n";
        }
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

