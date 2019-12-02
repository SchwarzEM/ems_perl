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
    my $phobius              = q{};

    if ( $input =~ /\A Gene \t/xms ) {
        print "$input\n";
    }

    elsif (      ( exists $data[43] )
             and ( exists $data[31] )
             and ( exists $data[32] )
             and ( exists $data[33] )
             and ( exists $data[24] )
       	     and ( exists $data[44] )
       	     and ( exists $data[23] )
       	     and ( exists $data[16] )
       	     and ( exists $data[11] )) {

        $int_vs_non_int_logFC = $data[43];
        @int_TPMs             = @data[31..33];
        $L3i_TPM              = $data[24];
        $int_vs_non_int_FDR   = $data[44];
        $hum_perc_id          = $data[23];
        $ofind_summary        = $data[16];
        $phobius              = $data[11];

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
                  and ( $ofind_summary !~ /haemonchus/xms ) 
                  and ( $ofind_summary =~ /necator/xms )

              and ( ( $phobius =~ /\ASigP\z/xms ) or ( $phobius =~ /\bSigP;/xms ) ) 
           ) {
            if ( $input =~ /Intestinal_R01/ ) {
		$input =~ s/Intestinal_R01/Intestinal_R01_prev/g;
                print "$input\n";
            }
            else {
                $input =~ s/\A(\S+\t[^\t]*\t)(\t.+)\z/$1Intestinal_hookworm-specific$2/;
                print "$input\n";
            }
        }
    }
}

