#!/usr/bin/env perl

# Purpose: cut off lower values of TRA-like traits.

use strict;
use warnings;

my $P_COUNT     = 14;
my $RS_COUNT    =  4;
my $RandS_COUNT = 25;

# Gold standards:

# FEM_APIME	401	15	1.36351240796291	3	31
# FEM_BOMTE	416	14	0.823480971707547	6	32
# FEM_MELCO	411	15	0.999933337777482	6	35
# FEM_NASVI	405	13	0.866608892740484	6	37
# TRA_ANAOB	417	15	0.441163495191318	2	16
# TRA_BACOL	422	13	0.433318889370354	3	21
# TRA_CERCA	429	14	0.499982143494875	2	13
# TRA_LUCCU	377	16	0.666638890046248	6	32
# TRSF_DROGR	171	15	0.999933337777482	9	37
# TRSF_DROHY	201	15	0.93744140991188	8	31
# TRSF_DROME	197	12	1.71404085130696	12	50
# TRSF_DROMO	207	15	0.789432135150781	8	33
# TRSF_DROVI	199	11	0.647020763484501	9	41

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \S+ \t \d+ \t (\d+) \t \S+ \t (\d+) \t (\d+) /xms ) { 
        my $obs_p_count    = $1;
        my $obs_rs_count   = $2;
        my $obs_r_and_s_co = $3;
        if (     ( $obs_p_count    >= $P_COUNT     ) 
             and ( $obs_rs_count   >= $RS_COUNT    ) 
             and ( $obs_r_and_s_co >= $RandS_COUNT ) ) { 
            print "$input\n";
        }
    }
}

