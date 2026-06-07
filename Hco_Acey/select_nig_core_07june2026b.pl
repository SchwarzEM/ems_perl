#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    if (     ( $input =~ /HKU_nigoni_JU1422/xms ) 
         and ( $input =~ /nigoni_EG5268/xms )
         and ( $input =~ /nigoni_JU1421/xms )
         and ( $input =~ /nigoni_JU2484/xms )
         and ( $input =~ /nigoni_JU2617/xms )
         and ( $input =~ /nigoni_JU4356/xms )
         and ( $input =~ /nigoni_VSL2202/xms )
         and ( $input =~ /nigoni_YR106/xms )
         and ( $input =~ /nigoni_ZF1220/xms )
         and ( $input =~ /orig_nigoni_JU1422/xms )
         and (    ($input !~ /briggsae_AF16/xms )
               or ($input !~ /briggsae_ED3036/xms )
               or ($input !~ /briggsae_HK104/xms )
               or ($input !~ /briggsae_JU349/xms )
               or ($input !~ /briggsae_QR24/xms )
               or ($input !~ /briggsae_QX1410/xms )
               or ($input !~ /briggsae_VX34/xms )
             )
         and (    ($input =~ /briggsae_AF16/xms )
               or ($input =~ /briggsae_ED3036/xms )
               or ($input =~ /briggsae_HK104/xms )
               or ($input =~ /briggsae_JU349/xms )
               or ($input =~ /briggsae_QR24/xms )
               or ($input =~ /briggsae_QX1410/xms )
               or ($input =~ /briggsae_VX34/xms )
             )
       ) {
        print "$input\n";
    }
}


