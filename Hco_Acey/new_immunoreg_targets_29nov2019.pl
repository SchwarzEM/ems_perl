#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Statistics::Descriptive;
use Scalar::Util qw(looks_like_number);

while (my $input = <>) {
    chomp $input;
    my @data = split /\t/, $input;

    if ( $input =~ /\A Gene \t/xms ) {
        print "$input\n";
    }

    elsif (      ( exists $data[45] )
             and ( exists $data[31] )
             and ( exists $data[32] )
             and ( exists $data[33] )
       	     and ( exists $data[46] )
       	     and ( exists $data[23] )
       	     and ( exists $data[16] )
       	     and ( exists $data[11] ) ) {

        my $non_dex_int_vs_dex_int_logFC = $data[45];
        my @int_TPMs                     = @data[31..33];
        my $non_dex_int_vs_dex_int_FDR   = $data[46];
        my $hum_perc_id                  = $data[23];
        my $ofind_summary                = $data[16];
        my $phobius                      = $data[11];

        my $stat1 = Statistics::Descriptive::Full->new();
        $stat1->add_data(@int_TPMs);
        my $max_int_TPM = $stat1->max();

        if ( ( $non_dex_int_vs_dex_int_logFC =~ /\S+/xms ) 
                  and looks_like_number($non_dex_int_vs_dex_int_logFC)
                  and ( $non_dex_int_vs_dex_int_logFC >= 3.322 )   # 2 ^ 3.322 = 10

             and ( $non_dex_int_vs_dex_int_FDR =~ /\S+/xms ) 
                 and looks_like_number($non_dex_int_vs_dex_int_FDR) 
                 and ( $non_dex_int_vs_dex_int_FDR <= 0.05 )

              and ( $max_int_TPM =~ /\S+/xms ) and ( $max_int_TPM >= 10 )

              and ( $hum_perc_id =~ /\S+/xms ) and ( $hum_perc_id < 90 )

              and ( $ofind_summary =~ /\S+/xms ) 
                  and ( $ofind_summary =~ /haemonchus/xms ) 
                  and ( $ofind_summary =~ /necator/xms )

              and ( ( $phobius =~ /\ASigP\z/xms ) or ( $phobius =~ /\bSigP;/xms ) ) 
           ) {
            if ( $input =~ /Immunoregulated_R01/ ) {
                print "$input\n";
            }
            else {
                $input =~ s/\A(\S+\t[^\t]*\t[^\t]*\t)(\t.+)\z/$1Immunoregulated_R01_new$2/;
                print "$input\n";
            }
        }
        elsif ( $input =~ /Immunoregulated_R01/ ) {
            $input =~ s/Immunoregulated_R01/Immunoregulated_R01_prev/g;
            print "$input\n";
        }
    }
}

