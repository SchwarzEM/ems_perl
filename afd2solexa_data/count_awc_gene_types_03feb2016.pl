#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my @awc_types = qw( AWC_pool AWC_cell_1 AWC_cell_2 AWC_cell_3 AWC_cell_4 AWC_cell_5 AWC_cell_any );

my $header = "AWC_type\tGene_count\tProtcod_count\tnrRNA_count\tHkeep_count\tGPCR_count";

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (WBGene\S+) 
                       \t (\S+) \t (\S+) \t (\S+) \t (\S+) \t (\S+) \t (\S+) 
                       \t ([^\t]+) 
                       \t ([^\t]*) 
                       \t ([^\t]*) /xms ) {
        my $gene    = $1;
        my $awc_pool = $2;
        my $awc_cell_1  = $3;
        my $awc_cell_2  = $4;
        my $awc_cell_3  = $5;
        my $awc_cell_4  = $6;
        my $awc_cell_5  = $7;
        my $coding  = $8;
        my $hkeep   = $9;
        my $gpcrs   = $10;

        my $awc_cell_any = 0;

        if ( $awc_pool >= 0.1 ) {
            $data_ref->{'AWC_pool'}->{'gene'}->{$gene} = 1;

            if ( $coding =~ /protein/xms ) {
                $data_ref->{'AWC_pool'}->{'protcod'}->{$gene} = 1;
            }
            if ( $coding =~ /\A ncRNA \z/xms ) {
                $data_ref->{'AWC_pool'}->{'ncrna'}->{$gene} = 1;
            }

            if ($hkeep) {
                $data_ref->{'AWC_pool'}->{'hkeep'}->{$gene} = 1;
            }
            if ($gpcrs) {
                $data_ref->{'AWC_pool'}->{'gpcrs'}->{$gene} = 1;
            }

            $awc_cell_any = 1;
        }

        if ( $awc_cell_1 >= 0.1 ) {
            $data_ref->{'AWC_cell_1'}->{'gene'}->{$gene} = 1;

            if ( $coding =~ /protein/xms ) {
                $data_ref->{'AWC_cell_1'}->{'protcod'}->{$gene} = 1;
            }
            if ( $coding =~ /\A ncRNA \z/xms ) {
                $data_ref->{'AWC_cell_1'}->{'ncrna'}->{$gene} = 1;
            }

            if ($hkeep) {
                $data_ref->{'AWC_cell_1'}->{'hkeep'}->{$gene} = 1;
            }
            if ($gpcrs) {
                $data_ref->{'AWC_cell_1'}->{'gpcrs'}->{$gene} = 1;
            }

            $awc_cell_any = 1;
        }

        if ( $awc_cell_2 >= 0.1 ) {    
            $data_ref->{'AWC_cell_2'}->{'gene'}->{$gene} = 1;

            if ( $coding =~ /protein/xms ) {
                $data_ref->{'AWC_cell_2'}->{'protcod'}->{$gene} = 1;
            }
            if ( $coding =~ /\A ncRNA \z/xms ) {
                $data_ref->{'AWC_cell_2'}->{'ncrna'}->{$gene} = 1;
            }

            if ($hkeep) {
                $data_ref->{'AWC_cell_2'}->{'hkeep'}->{$gene} = 1;
            }
            if ($gpcrs) {
                $data_ref->{'AWC_cell_2'}->{'gpcrs'}->{$gene} = 1;
            }

            $awc_cell_any = 1;
        }

        if ( $awc_cell_3 >= 0.1 ) {    
            $data_ref->{'AWC_cell_3'}->{'gene'}->{$gene} = 1;

            if ( $coding =~ /protein/xms ) {
                $data_ref->{'AWC_cell_3'}->{'protcod'}->{$gene} = 1;
            }
            if ( $coding =~ /\A ncRNA \z/xms ) {
                $data_ref->{'AWC_cell_3'}->{'ncrna'}->{$gene} = 1;
            }

            if ($hkeep) {
                $data_ref->{'AWC_cell_3'}->{'hkeep'}->{$gene} = 1;
            }
            if ($gpcrs) {
                $data_ref->{'AWC_cell_3'}->{'gpcrs'}->{$gene} = 1;
            }

            $awc_cell_any = 1;
        }

        if ( $awc_cell_4 >= 0.1 ) {    
            $data_ref->{'AWC_cell_4'}->{'gene'}->{$gene} = 1;

            if ( $coding =~ /protein/xms ) {
                $data_ref->{'AWC_cell_4'}->{'protcod'}->{$gene} = 1;
            }
            if ( $coding =~ /\A ncRNA \z/xms ) {
                $data_ref->{'AWC_cell_4'}->{'ncrna'}->{$gene} = 1;
            }

            if ($hkeep) {
                $data_ref->{'AWC_cell_4'}->{'hkeep'}->{$gene} = 1;
            }
            if ($gpcrs) {
                $data_ref->{'AWC_cell_4'}->{'gpcrs'}->{$gene} = 1;
            }

            $awc_cell_any = 1;
        }

        if ( $awc_cell_5 >= 0.1 ) {    
            $data_ref->{'AWC_cell_5'}->{'gene'}->{$gene} = 1;

            if ( $coding =~ /protein/xms ) {
                $data_ref->{'AWC_cell_5'}->{'protcod'}->{$gene} = 1;
            }
            if ( $coding =~ /\A ncRNA \z/xms ) {
                $data_ref->{'AWC_cell_5'}->{'ncrna'}->{$gene} = 1;
            }

            if ($hkeep) {
                $data_ref->{'AWC_cell_5'}->{'hkeep'}->{$gene} = 1;
            }
            if ($gpcrs) {
                $data_ref->{'AWC_cell_5'}->{'gpcrs'}->{$gene} = 1;
            }

            $awc_cell_any = 1;
        }

        if ($awc_cell_any) {
            $data_ref->{'AWC_cell_any'}->{'gene'}->{$gene} = 1;

            if ( $coding =~ /protein/xms ) {
                $data_ref->{'AWC_cell_any'}->{'protcod'}->{$gene} = 1;
            }
            if ( $coding =~ /\A ncRNA \z/xms ) {
                $data_ref->{'AWC_cell_any'}->{'ncrna'}->{$gene} = 1;
            }

            if ($hkeep) {
                $data_ref->{'AWC_cell_any'}->{'hkeep'}->{$gene} = 1;
            }
            if ($gpcrs) {
                $data_ref->{'AWC_cell_any'}->{'gpcrs'}->{$gene} = 1;
            }

            $awc_cell_any = 0;
        }
    }
    elsif ( $input !~ /\A Gene /xms ) {
        die "Cannot parse input line: $input\n";
    }
}

foreach my $awc_type (@awc_types) {
    if ( exists $data_ref->{$awc_type}->{'gene'} ) {
        my @genes = sort keys %{ $data_ref->{$awc_type}->{'gene'} };
        my $gene_count = @genes;

        my @protcods = sort keys %{ $data_ref->{$awc_type}->{'protcod'} };
        my $protcod_count = @protcods;

        my @ncrnas = sort keys %{ $data_ref->{$awc_type}->{'ncrna'} };
        my $ncrna_count = @ncrnas;

        my @hkeeps = sort keys %{ $data_ref->{$awc_type}->{'hkeep'} };
        my $hkeep_count = @hkeeps;

        my @gpcrs = sort keys %{ $data_ref->{$awc_type}->{'gpcrs'} };
        my $gpcr_count = @gpcrs;

        print "$header\n" if $header;
        $header = q{};
        print "$awc_type\t$gene_count\t$protcod_count\t$ncrna_count\t$hkeep_count\t$gpcr_count\n";
    }
}

