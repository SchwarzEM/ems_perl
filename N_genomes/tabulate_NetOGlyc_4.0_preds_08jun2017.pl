#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Gene\tOGlyc_4.0";

my %bad2good_names = ( 
    "01" => "Csp5_scaffold_00069.alt.g2455",
    "02" => "brenneri_WS250.HM2_0139.g27725",
    "03" => "brenneri_WS250.HM2_0011.g6190",
    "04" => "CBG20255",
    "05" => "Csp5_scaffold_04733.alt.g27387",
    "06" => "Cnig_chr_I.g2419",
    "07" => "Csp5_scaffold_00375.alt.g8079",
    "08" => "FL81_17792",
    "09" => "FL81_17791",
    "10" => "brenneri_WS250.HM2_0010.g5865",
    "11" => "FL81_26171",
    "12" => "D2062.7",
    "13" => "FL81_17821",
    "14" => "Csp5_scaffold_04733.alt.g27386",
    "15" => "FL81_19279",
    "16" => "brenneri_WS250.HM2_0011.g6193",
    "17" => "brenneri_WS250.HM2_0011.g6141",
    "18" => "brenneri_WS250.HM2_0033.g13691",
    "19" => "tropicalis_2016.08.11_015.g12484",
    "20" => "Cnig_chr_III.g11664",
    "21" => "R06F6.7",
    "22" => "brenneri_WS250.HM2_0011.g6140",
    "23" => "wallacei_2016.07.09_01.g1181",
    "24" => "CBG25787",
    "25" => "FL81_17790",
    "26" => "FL81_15843",
    "27" => "Cnig_chr_III.g11663",
    "28" => "Cnig_chr_IV.g14916",
    "29" => "Cnig_chr_I.g2786",
    "30" => "F46F5.6",
    "31" => "tropicalis_2016.08.11_034.g18194",
    "32" => "F46F5.9",
    "33" => "brenneri_WS250.HM2_0011.g6194",
    "34" => "Cnig_chr_I.g2781",
    "35" => "brenneri_WS250.HM2_0011.g6145",
    "36" => "brenneri_WS250.HM2_0011.g6195",
    "37" => "Cnig_chr_IV.g13698",
    "38" => "Csp5_scaffold_00036.alt.g1521",
    "39" => "ZK484.5",
    "40" => "Csp5_scaffold_00375.alt.g8081",
    "41" => "CBG25393",
    "42" => "CBG05758",
    "43" => "Cnig_chr_II.g5744",
    "44" => "Csp5_scaffold_00107.alt.g3418",
    "45" => "CBG25574",
    "46" => "F32B5.4",
    "47" => "Csp5_scaffold_02717.alt.g22336",
    "48" => "FL81_05289",
    "49" => "Cnig_chr_II.g5743",
    "50" => "CBG03021",
    "51" => "D2062.6",
    "52" => "wallacei_2016.07.09_02.g5526",
    "53" => "FL81_22033",
    "54" => "Cnig_chr_I.g2782",
    "55" => "Csp5_scaffold_00420.alt.g8676",
    "56" => "wallacei_2016.07.09_02.g5525",
    "57" => "FL81_17794",
    "58" => "CBG25392",
    "59" => "wallacei_2016.07.09_02.g5524",
    "60" => "brenneri_WS250.HM2_0011.g6139",
    "61" => "Cnig_chr_III.g11661",
    "62" => "Csp5_scaffold_00069.alt.g2453",
    "63" => "Cnig_chr_I.g2784",
    "64" => "Cnig_chr_I.g2783",
    "65" => "Csp5_scaffold_00069.alt.g2454",
    "66" => "CBG25391",
    "67" => "CBG25390",
    "68" => "FL81_22031",
    "69" => "CBG25394",
);

my %ok_names = ( 
    "brenneri_WS250.HM2_0010.g5865" => 1,
    "brenneri_WS250.HM2_0011.g6139" => 1,
    "brenneri_WS250.HM2_0011.g6140" => 1,
    "brenneri_WS250.HM2_0011.g6141" => 1,
    "brenneri_WS250.HM2_0011.g6145" => 1,
    "brenneri_WS250.HM2_0011.g6190" => 1,
    "brenneri_WS250.HM2_0011.g6193" => 1,
    "brenneri_WS250.HM2_0011.g6194" => 1,
    "brenneri_WS250.HM2_0011.g6195" => 1,
    "brenneri_WS250.HM2_0033.g13691" => 1,
    "brenneri_WS250.HM2_0139.g27725" => 1,
    "CBG03021" => 1,
    "CBG05758" => 1,
    "CBG20255" => 1,
    "CBG25390" => 1,
    "CBG25391" => 1,
    "CBG25392" => 1,
    "CBG25393" => 1,
    "CBG25394" => 1,
    "CBG25574" => 1,
    "CBG25787" => 1,
    "Cnig_chr_I.g2419" => 1,
    "Cnig_chr_I.g2781" => 1,
    "Cnig_chr_I.g2782" => 1,
    "Cnig_chr_I.g2783" => 1,
    "Cnig_chr_I.g2784" => 1,
    "Cnig_chr_I.g2786" => 1,
    "Cnig_chr_II.g5743" => 1,
    "Cnig_chr_II.g5744" => 1,
    "Cnig_chr_III.g11661" => 1,
    "Cnig_chr_III.g11663" => 1,
    "Cnig_chr_III.g11664" => 1,
    "Cnig_chr_IV.g13698" => 1,
    "Cnig_chr_IV.g14916" => 1,
    "Csp5_scaffold_00036.alt.g1521" => 1,
    "Csp5_scaffold_00069.alt.g2453" => 1,
    "Csp5_scaffold_00069.alt.g2454" => 1,
    "Csp5_scaffold_00069.alt.g2455" => 1,
    "Csp5_scaffold_00107.alt.g3418" => 1,
    "Csp5_scaffold_00375.alt.g8079" => 1,
    "Csp5_scaffold_00375.alt.g8081" => 1,
    "Csp5_scaffold_00420.alt.g8676" => 1,
    "Csp5_scaffold_02717.alt.g22336" => 1,
    "Csp5_scaffold_04733.alt.g27386" => 1,
    "Csp5_scaffold_04733.alt.g27387" => 1,
    "D2062.6" => 1,
    "D2062.7" => 1,
    "F32B5.4" => 1,
    "F46F5.6" => 1,
    "F46F5.9" => 1,
    "R06F6.7" => 1,
    "ZK484.5" => 1,
    "FL81_05289" => 1,
    "FL81_15843" => 1,
    "FL81_17790" => 1,
    "FL81_17791" => 1,
    "FL81_17792" => 1,
    "FL81_17794" => 1,
    "FL81_17821" => 1,
    "FL81_19279" => 1,
    "FL81_22031" => 1,
    "FL81_22033" => 1,
    "tropicalis_2016.08.11_015.g12484" => 1,
    "tropicalis_2016.08.11_034.g18194" => 1,
    "wallacei_2016.07.09_01.g1181" => 1,
    "wallacei_2016.07.09_02.g5524" => 1,
    "wallacei_2016.07.09_02.g5525" => 1,
    "wallacei_2016.07.09_02.g5526" => 1,
);

my %isoforms = ( 
    "01" => "Csp5_scaffold_00069.alt.g2455.t4",
    "02" => "brenneri_WS250.HM2_0139.g27725.t1",
    "03" => "brenneri_WS250.HM2_0011.g6190.t2",
    "04" => "CBG20255",
    "05" => "Csp5_scaffold_04733.alt.g27387.t4",
    "06" => "Cnig_chr_I.g2419.t3",
    "07" => "Csp5_scaffold_00375.alt.g8079.t1",
    "08" => "FL81_17792-RA",
    "09" => "FL81_17791-RA",
    "10" => "brenneri_WS250.HM2_0010.g5865.t1",
    "11" => "FL81_26171-RA",
    "12" => "D2062.7",
    "13" => "FL81_17821-RA",
    "14" => "Csp5_scaffold_04733.alt.g27386.t1",
    "15" => "FL81_19279-RA",
    "16" => "brenneri_WS250.HM2_0011.g6193.t1",
    "17" => "brenneri_WS250.HM2_0011.g6141.t1",
    "18" => "brenneri_WS250.HM2_0033.g13691.t1",
    "19" => "tropicalis_2016.08.11_015.g12484.t1",
    "20" => "Cnig_chr_III.g11664.t1",
    "21" => "R06F6.7",
    "22" => "brenneri_WS250.HM2_0011.g6140.t2",
    "23" => "wallacei_2016.07.09_01.g1181.t1",
    "24" => "CBG25787",
    "25" => "FL81_17790-RA",
    "26" => "FL81_15843-RA",
    "27" => "Cnig_chr_III.g11663.t1",
    "28" => "Cnig_chr_IV.g14916.t1",
    "29" => "Cnig_chr_I.g2786.t1",
    "30" => "F46F5.6",
    "31" => "tropicalis_2016.08.11_034.g18194.t2",
    "32" => "F46F5.9",
    "33" => "brenneri_WS250.HM2_0011.g6194.t1",
    "34" => "Cnig_chr_I.g2781.t1",
    "35" => "brenneri_WS250.HM2_0011.g6145.t1",
    "36" => "brenneri_WS250.HM2_0011.g6195.t1",
    "37" => "Cnig_chr_IV.g13698.t1",
    "38" => "Csp5_scaffold_00036.alt.g1521.t1",
    "39" => "ZK484.5",
    "40" => "Csp5_scaffold_00375.alt.g8081.t1",
    "41" => "CBG25393",
    "42" => "CBG05758",
    "43" => "Cnig_chr_II.g5744.t1",
    "44" => "Csp5_scaffold_00107.alt.g3418.t1",
    "45" => "CBG25574",
    "46" => "F32B5.4",
    "47" => "Csp5_scaffold_02717.alt.g22336.t1",
    "48" => "FL81_05289-RA",
    "49" => "Cnig_chr_II.g5743.t1",
    "50" => "CBG03021",
    "51" => "D2062.6",
    "52" => "wallacei_2016.07.09_02.g5526.t1",
    "53" => "FL81_22033-RA",
    "54" => "Cnig_chr_I.g2782.t1",
    "55" => "Csp5_scaffold_00420.alt.g8676.t1",
    "56" => "wallacei_2016.07.09_02.g5525.t1",
    "57" => "FL81_17794-RA",
    "58" => "CBG25392",
    "59" => "wallacei_2016.07.09_02.g5524.t1",
    "60" => "brenneri_WS250.HM2_0011.g6139.t1",
    "61" => "Cnig_chr_III.g11661.t1",
    "62" => "Csp5_scaffold_00069.alt.g2453.t1",
    "63" => "Cnig_chr_I.g2784.t1",
    "64" => "Cnig_chr_I.g2783.t1",
    "65" => "Csp5_scaffold_00069.alt.g2454.t1",
    "66" => "CBG25391",
    "67" => "CBG25390",
    "68" => "FL81_22031-RA",
    "69" => "CBG25394",
);

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\d+) \t netOGlyc-4.0.0.13 \t CARBOHYD \t (\d+) \t (\d+) \t \S+ \t [.] \t [.] \t [#]POSITIVE \s* \z/xms ) {
        my $bad_name  = $1;
        my $start_res = $2;
        my $end_res   = $3;

        if ( $start_res != $end_res ) {
            die "Failed to identify single residue site in $input\n";
        }
        if (! exists $bad2good_names{$bad_name} ) { 
            die "Cannot parse sequence index number (\"$bad_name\") in: $input\n";
        }
        my $seq_name = $bad2good_names{$bad_name};
        my $isoform  = $isoforms{$bad_name};

        if (! exists $ok_names{$seq_name} ) {
            if ( $seq_name ne 'FL81_26171' ) { 
                die "Cannot recognize sequence name (\"$seq_name\") in: $input\n";
            }
        }
        else { 
            $data_ref->{'seq_name'}->{$seq_name}->{'o_glyc_site'}->{$start_res} = 1;
            $data_ref->{'seq_name'}->{$seq_name}->{'isoform'} = $isoform;
        }
    }
    elsif ( $input =~ /POSITIVE/xms ) { 
        die "Failed to parse: $input\n";
    }
}

my @output_seqs = sort keys %{ $data_ref->{'seq_name'} };
foreach my $output_seq (@output_seqs) {
    my @residues = sort { $a <=> $b } keys %{ $data_ref->{'seq_name'}->{$output_seq}->{'o_glyc_site'} };
    my $residue_text = join q{,}, @residues;
    my $isoform = $data_ref->{'seq_name'}->{$output_seq}->{'isoform'};

    print "$header\n" if $header;
    $header = q{};

    print "$output_seq\t$residue_text [$isoform]\n";
}

