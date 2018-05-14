#!/usr/bin/env perl

use strict;
use warnings;

print "Gene",
      "\t",
      "Secreted_100\tSecreted_150\tSecreted_200\tSecreted_Any",
      "\t",
      "Cons_Secreted_100\tCons_Secreted_150\tCons_Secreted_200\tCons_Secreted_Any",
      "\t",
      "Hco_Cons_Secreted_100\tHco_Cons_Secreted_150\tHco_Cons_Secreted_200\tHco_Cons_Secreted_Any",
      "\t",
      "Hco_Cons_Only_Secreted_100\tHco_Cons_Only_Secreted_150\tHco_Cons_Only_Secreted_200\tHco_Cons_Only_Secreted_Any",
      "\n",
      ;

while (my $input = <>) { 
    chomp $input;
    # Gene    Max_prot_size   Start/stop codons       Phobius Omcl_Summary [14 spp.]

    # We require 'Complete' for everything, particularly the <100-residue peptides.
    if ( $input =~ /\A (Acey_2012.08.05_\S+) \t (\d+) \t Complete \t SigP \t ([^\t]*) \z/xms ) { 
        my $gene         = $1;
        my $size         = $2;
        my $conservation = $3;

        my $annot_100 = q{};
        my $annot_150 = q{};
        my $annot_200 = q{};
        my $annot_any = q{};

        my $cons_annot_100 = q{};
        my $cons_annot_150 = q{};
        my $cons_annot_200 = q{};
        my $cons_annot_any = q{};

        my $hco_cons_annot_100 = q{};
        my $hco_cons_annot_150 = q{};
        my $hco_cons_annot_200 = q{};
        my $hco_cons_annot_any = q{};

        my $hco_cons_only_annot_100 = q{};
        my $hco_cons_only_annot_150 = q{};
        my $hco_cons_only_annot_200 = q{};
        my $hco_cons_only_annot_any = q{};

        $annot_100 = 'Secreted_100' if ( $size <= 100 );
        $annot_150 = 'Secreted_150' if ( $size <= 150 );
        $annot_200 = 'Secreted_200' if ( $size <= 200 );
        $annot_any = 'Secreted_Any';

        if ( $conservation =~ /\A ORTHOMCL\d+.+ \z/xms ) { 
            $cons_annot_100 = 'Cons_Secreted_100' if ( $size <= 100 );
            $cons_annot_150 = 'Cons_Secreted_150' if ( $size <= 150 );
            $cons_annot_200 = 'Cons_Secreted_200' if ( $size <= 200 );
            $cons_annot_any = 'Cons_Secreted_Any';
        }

        if ( $conservation =~ / haemonchus /xms ) {
            $hco_cons_annot_100 = 'Hco_Cons_Secreted_100' if ( $size <= 100 );
            $hco_cons_annot_150 = 'Hco_Cons_Secreted_150' if ( $size <= 150 );
            $hco_cons_annot_200 = 'Hco_Cons_Secreted_200' if ( $size <= 200 );
            $hco_cons_annot_any = 'Hco_Cons_Secreted_any';
        }

        if (     ( $conservation =~ / ,3[ ]taxa .+ ancylostoma .+ haemonchus .+ haemonchus_aug /xms ) 
              or ( $conservation =~ / ,2[ ]taxa .+ ancylostoma .+ haemonchus /xms ) 
           ) { 
            $hco_cons_only_annot_100 = 'Hco_Cons_Only_Secreted_100' if ( $size <= 100 );
            $hco_cons_only_annot_150 = 'Hco_Cons_Only_Secreted_150' if ( $size <= 150 );
            $hco_cons_only_annot_200 = 'Hco_Cons_Only_Secreted_200' if ( $size <= 200 );
            $hco_cons_only_annot_any = 'Hco_Cons_Only_Secreted_Any';
        }

        print "$gene",
              "\t",
              "$annot_100\t$annot_150\t$annot_200\t$annot_any",
              "\t",
              "$cons_annot_100\t$cons_annot_150\t$cons_annot_200\t$cons_annot_any",
              "\t",
              "$hco_cons_annot_100\t$hco_cons_annot_150\t$hco_cons_annot_200\t$hco_cons_annot_any",
              "\t",
              "$hco_cons_only_annot_100\t$hco_cons_only_annot_150\t$hco_cons_only_annot_200\t$hco_cons_only_annot_any",
              "\n",
              ;
    }
}

