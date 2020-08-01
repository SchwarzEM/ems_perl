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

    elsif (      ( exists $data[53] )
             and ( exists $data[39] )
             and ( exists $data[40] )
             and ( exists $data[41] )
             and ( exists $data[32] )
       	     and ( exists $data[54] )
       	     and ( exists $data[31] )
       	     and ( exists $data[22] )
             and ( exists $data[25] )
       	     and ( exists $data[17] ) ) {

        my $non_dex_int_vs_dex_int_logFC = $data[53];      # intestines_19dpi_noDEX.vs.intestines_19dpi_DEX.logFC [54]
        my @int_TPMs                     = @data[39..41];  # intestines_19dpi_noDEX_biorep_1.TPM [40],
                                                           #     intestines_19dpi_noDEX_biorep_2.TPM [41], intestines_19dpi_noDEX_biorep_3.TPM [42]
        my $L3i_TPM                      = $data[32];      # L3i.TPM [33]
        my $non_dex_int_vs_dex_int_FDR   = $data[54];      # intestines_19dpi_noDEX.vs.intestines_19dpi_DEX.FDR [55]
        my $hum_perc_id                  = $data[31];      # Human_perc.id [32]
        my $ofind_summary                = $data[22];      # OFind_Summary [23]
        my $homol_deg                    = $data[25];      # Homol_DEG [26]
        my $phobius                      = $data[17];      # Phobius [18]

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

              and ( $L3i_TPM =~ /\S+/xms ) and ( $L3i_TPM < 50 )

              and ( $hum_perc_id =~ /\S+/xms ) and ( $hum_perc_id < 90 )

              and ( $ofind_summary =~ /\S+/xms ) 
                  and ( $ofind_summary =~ /haemonchus/xms ) 
                  and ( $ofind_summary =~ /necator/xms )

              and ( $homol_deg =~ /\S+/xms )

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
    }
}

