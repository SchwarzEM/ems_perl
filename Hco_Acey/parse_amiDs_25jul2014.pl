#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %prefix2class = (
    Acey => 'Nematodes',
    Acpis => 'Arthropods',
    Acyrt => 'Arthropods',
    Aeaeg => 'Arthropods',
    Amel => 'Arthropods',
    Angam => 'Arthropods',
    Apis => 'Arthropods',
    Asuu => 'Nematodes',
    Atcep => 'Arthropods',
    Bflo => 'Metazoa',
    Bmal => 'Nematodes',
    Bomor => 'Arthropods',
    Bxyl => 'Nematodes',
    Cbri => 'Nematodes',
    Cele => 'Nematodes',
    Cerat => 'Arthropods',
    Cfam => 'Vertebrates',
    Cgig => 'Metazoa',
    Cint => 'Metazoa',
    Ctel => 'Metazoa',
    Cuqui => 'Arthropods',
    Daple => 'Arthropods',
    Dapul => 'Arthropods',
    Depon => 'Arthropods',
    Diaph => 'Arthropods',
    Dimm => 'Nematodes',
    Dmel => 'Arthropods',
    Dnov => 'Vertebrates',
    Drer => 'Vertebrates',
    Eeur => 'Vertebrates',
    Ggal => 'Vertebrates',
    Hbac => 'Nematodes',
    Hcon => 'Nematodes',
    Hemel => 'Arthropods',
    Hpol => 'Nematodes',
    Hrob => 'Metazoa',
    Hsap => 'Vertebrates',
    Hydra => 'Metazoa',
    Ixsca => 'Arthropods',
    Lafr => 'Vertebrates',
    Lcha => 'Vertebrates',
    Lgig => 'Metazoa',
    Lloa => 'Nematodes',
    Lone => 'Nematodes',
    Maur => 'Vertebrates',
    Mdom => 'Vertebrates',
    Mesca => 'Arthropods',
    Mhap => 'Nematodes',
    Mmus => 'Vertebrates',
    Name => 'Nematodes',
    Navit => 'Arthropods',
    Nbra => 'Nematodes',
    Oana => 'Vertebrates',
    Phum => 'Arthropods',
    Plano => 'Arthropods',
    Pmar => 'Vertebrates',
    Ppac => 'Nematodes',
    Pred => 'Nematodes',
    Rhpro => 'Arthropods',
    Skow => 'Metazoa',
    Soinv => 'Arthropods',
    Spur => 'Metazoa',
    Stmar => 'Arthropods',
    Tcan => 'Nematodes',
    Teurt => 'Arthropods',
    Tgut => 'Vertebrates',
    Trcas => 'Arthropods',
    Trich => 'Metazoa',
    Tspi => 'Nematodes',
    Wban => 'Nematodes',
    Xtro => 'Vertebrates',
);

while (my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A > /xms ) {
        my $seqname = q{};
        my $prefix  = q{};
        if ( $input =~ /\A > (([A-Za-z]+) _ \S+) /xms ) { 
            $seqname = $1;
            $prefix  = $2;
            if (! exists $prefix2class{$prefix} ) {
                die "Cannot parse prefix \"$prefix\" in: $input\n";
            }
            else {
                my $type = $prefix2class{$prefix};
                print "$type\t$seqname\n";
            }
        }
        else { 
            die "Cannot parse sequence name in: $input\n";
        }
    }
}


