#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Statistics::Descriptive;
use Scalar::Util qw(looks_like_number);

while (my $input = <>) {
    chomp $input;
    my @data = split /\t/, $input;

    my $int_vs_non_int_logFC = q{};
    my @int_TPMs             = ();
    my $L3i_TPM              = q{};
    my $int_vs_non_int_FDR   = q{};
    my $hum_perc_id          = q{};
    my $ofind_summary        = q{};
    my $homol_deg            = q{};
    my $phobius              = q{};

    if ( $input =~ /\A Gene \t/xms ) {
        print "$input\n";
    }

    elsif (      ( exists $data[51] )
             and ( exists $data[39] )
             and ( exists $data[40] )
             and ( exists $data[41] )
             and ( exists $data[32] )
       	     and ( exists $data[52] )
       	     and ( exists $data[31] )
       	     and ( exists $data[22] )
             and ( exists $data[25] )
       	     and ( exists $data[17] )) {

        $int_vs_non_int_logFC = $data[51];      # intestines_19dpi_noDEX.vs.non.intestines_19dpi_noDEX.logFC [52]
        @int_TPMs             = @data[39..41];  # intestines_19dpi_noDEX_biorep_1.TPM [40],
                                                #     intestines_19dpi_noDEX_biorep_2.TPM [41], intestines_19dpi_noDEX_biorep_3.TPM [42]
        $L3i_TPM              = $data[32];      # L3i.TPM [33]
        $int_vs_non_int_FDR   = $data[52];      # intestines_19dpi_noDEX.vs.non.intestines_19dpi_noDEX.FDR [53]
        $hum_perc_id          = $data[31];      # Human_perc.id [32]
        $ofind_summary        = $data[22];      # OFind_Summary [23]
        $homol_deg            = $data[25];      # Homol_DEG [26]
        $phobius              = $data[17];      # Phobius [18]

        my $stat1 = Statistics::Descriptive::Full->new();
        $stat1->add_data(@int_TPMs);
        my $max_int_TPM = $stat1->max();

        if ( ( $int_vs_non_int_logFC =~ /\S+/xms ) 
                  and looks_like_number($int_vs_non_int_logFC)
                  and ( $int_vs_non_int_logFC >= 3.322 )   # 2 ^ 3.322 = 10

             and ( $int_vs_non_int_FDR =~ /\S+/xms ) 
                 and looks_like_number($int_vs_non_int_FDR) 
                 and ( $int_vs_non_int_FDR <= 0.05 )

              and ( $max_int_TPM =~ /\S+/xms ) and ( $max_int_TPM >= 10 )

              and ( $L3i_TPM =~ /\S+/xms ) and ( $L3i_TPM < 50 )

              and ( $hum_perc_id =~ /\S+/xms ) and ( $hum_perc_id < 90 )

              and ( $ofind_summary =~ /\S+/xms ) 
                  and ( $ofind_summary =~ /haemonchus/xms ) 
                  and ( $ofind_summary =~ /necator/xms )

              and ( $homol_deg =~ /\S+/xms )

              and ( ( $phobius =~ /\ASigP\z/xms ) or ( $phobius =~ /\bSigP;/xms ) ) 
           ) {
            if ( $input =~ /Intestinal_R01/ ) {
                print "$input\n";
            }
            else {
                $input =~ s/\A(\S+\t[^\t]*\t)(\t.+)\z/$1Intestinal_R01_new$2/;
                print "$input\n";
            }
        }
    }
}

