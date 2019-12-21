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

    elsif (      ( exists $data[48] )
             and ( exists $data[34] )
             and ( exists $data[35] )
             and ( exists $data[36] )
             and ( exists $data[27] )
       	     and ( exists $data[49] )
       	     and ( exists $data[26] )
       	     and ( exists $data[17] )
             and ( exists $data[20] )
       	     and ( exists $data[12] ) ) {

        my $non_dex_int_vs_dex_int_logFC = $data[48];      # intestines_19dpi_noDEX.vs.intestines_19dpi_DEX.logFC [49]
        my @int_TPMs                     = @data[34..36];  # intestines_19dpi_noDEX_biorep_1.TPM [35],
                                                           #     intestines_19dpi_noDEX_biorep_2.TPM [36], intestines_19dpi_noDEX_biorep_3.TPM [37]
        my $L3i_TPM                      = $data[27];      # L3i.TPM [28]
        my $non_dex_int_vs_dex_int_FDR   = $data[49];      # intestines_19dpi_noDEX.vs.intestines_19dpi_DEX.FDR [50]
        my $hum_perc_id                  = $data[26];      # Human_perc.id [27]
        my $ofind_summary                = $data[17];      # OFind_Summary [18]
        my $homol_deg                    = $data[20];      # Homol_DEG [21]
        my $phobius                      = $data[12];      # Phobius [13]

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

