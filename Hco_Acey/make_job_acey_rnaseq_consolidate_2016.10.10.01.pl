#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) {
    # sample inputs:
    # 12990_Alb_ele.4hr.L4_retro_06aug2012.no_cont.min40.R1.se.fq.gz
    # 12664_ACEY.5.D_retro_21may2012.no_cont.min40.se.fq.gz

    chomp $input;
    if ( $input =~ /\A (\S+) \.no_cont\.min40\.pe\.R1\.fq\.gz \z/xms ) {
        my $stem = $1;
        $data_ref->{'stem_pe'}->{$stem}->{'pe'}->{$input} = 1;
    }
    elsif ($input =~ /\A (\S+) \.no_cont\.min40 (?: \.R1){0,1} \.se\.fq\.gz \z/xms ) {
        my $stem = $1;
        $data_ref->{'stem_se'}->{$stem}->{'se'}->{$input} = 1;
    }
    else { 
        die "Cannot parse input: $input\n";
    }
}

my @stem_pes = sort keys %{ $data_ref->{'stem_pe'} };
foreach my $stem_pe (@stem_pes) {
    my @pe_reads = ();
    @pe_reads = sort keys %{ $data_ref->{'stem_pe'}->{$stem_pe}->{'pe'} };

    my @se_reads = ();
    if ( exists $data_ref->{'stem_se'}->{$stem_pe}->{'se'} ) {
        @se_reads = sort keys %{ $data_ref->{'stem_se'}->{$stem_pe}->{'se'} };
    }
    my $output = "$stem_pe.no_cont.min40.all.fq";
    $output    = safename($output);

    my $used_reads_dir = "used_reads_2016.10.10_dir";
    $used_reads_dir    = safename($used_reads_dir);

    print $header if $header;
    print "mkdir $used_reads_dir ;\n\n" if $header;
    $header = q{};

    print "zcat @pe_reads @se_reads > $output ;\n";
    print "mv -i @pe_reads @se_reads $used_reads_dir ;\n";
    print "\n";
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
