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

    elsif (      ( exists $data[45] )
             and ( exists $data[33] )
             and ( exists $data[34] )
             and ( exists $data[35] )
             and ( exists $data[26] )
       	     and ( exists $data[46] )
       	     and ( exists $data[25] )
       	     and ( exists $data[16] )
             and ( exists $data[19] )
       	     and ( exists $data[11] )) {

        $int_vs_non_int_logFC = $data[45];      # intestines_19dpi_noDEX.vs.non.intestines_19dpi_noDEX.logFC [46]
        @int_TPMs             = @data[33..35];  # intestines_19dpi_noDEX_biorep_1.TPM [34], 
                                                #     intestines_19dpi_noDEX_biorep_2.TPM [35], intestines_19dpi_noDEX_biorep_3.TPM [36]
        $L3i_TPM              = $data[26];      # L3i.TPM [27]
        $int_vs_non_int_FDR   = $data[46];      # intestines_19dpi_noDEX.vs.non.intestines_19dpi_noDEX.FDR [47]
        $hum_perc_id          = $data[25];      # Human_perc.id [26]
        $ofind_summary        = $data[16];      # OFind_Summary [17]
        $homol_deg            = $data[19];      # Homol_DEG [20]
        $phobius              = $data[11];      # Phobius [12]

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

