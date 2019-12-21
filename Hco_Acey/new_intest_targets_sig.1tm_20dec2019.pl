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

    elsif (      ( exists $data[46] )
             and ( exists $data[34] )
             and ( exists $data[35] )
             and ( exists $data[36] )
             and ( exists $data[27] )
       	     and ( exists $data[47] )
       	     and ( exists $data[26] )
       	     and ( exists $data[17] )
             and ( exists $data[20] )
       	     and ( exists $data[12] )) {

        $int_vs_non_int_logFC = $data[46];      # intestines_19dpi_noDEX.vs.non.intestines_19dpi_noDEX.logFC [47]
        @int_TPMs             = @data[34..36];  # intestines_19dpi_noDEX_biorep_1.TPM [35],
                                                #     intestines_19dpi_noDEX_biorep_2.TPM [36], intestines_19dpi_noDEX_biorep_3.TPM [37]
        $L3i_TPM              = $data[27];      # L3i.TPM [28]
        $int_vs_non_int_FDR   = $data[47];      # intestines_19dpi_noDEX.vs.non.intestines_19dpi_noDEX.FDR [48]
        $hum_perc_id          = $data[26];      # Human_perc.id [27]
        $ofind_summary        = $data[17];      # OFind_Summary [18]
        $homol_deg            = $data[20];      # Homol_DEG [21]
        $phobius              = $data[12];      # Phobius [13]

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

              and ( ( $phobius =~ /SigP[+]TM\(1x\)/xms ) and ( $phobius !~ /SigP;/xms ) )
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

