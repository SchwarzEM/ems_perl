#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Condition\tUpreg. genes\tDownreg. genes";

my @conds = (
    'Col_WT vs. lgo-2',
    'Col_WT vs. atml1-3',
    'atml1-3 vs. lgo-2',
    'LGOoe vs. lgo-2',
    'LGOoe vs. atml1-3',
    'LGOoe vs. Col_WT',
    'LGOoe atml1-3 vs. lgo-2',
    'LGOoe atml1-3 vs. atml1-3',
    'LGOoe atml1-3 vs. Col_WT',
    'LGOoe vs. LGOoe atml1-3',
    'batch_two vs. batch_one',  
    'batch_three vs. batch_one',
    'batch_three vs. batch_two',
);

foreach my $cond (@conds) {
    $data_ref->{'cond_ok'}->{$cond} = 1;
}

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) (up|down) reg\. \t (\d+) \z/xms ) { 
        my $cond  = $1;
        my $dir   = $2;
        my $genes = $3;
        $cond =~ s/\A\s+//;
        $cond =~ s/\s+\z//;
        $cond =~ s/ATML1::LGO/LGOoe/g;
        $data_ref->{'cond'}->{$cond}->{$dir} = $genes;
    }
    else {
        die "Cannot parse: $input\n";
    }
}

foreach my $cond (@conds) {
    my $up_genes   = $data_ref->{'cond'}->{$cond}->{'up'};
    my $down_genes = $data_ref->{'cond'}->{$cond}->{'down'};

    print "$header\n" if $header;
    $header = q{};

    $up_genes   = commify($up_genes);
    $down_genes = commify($down_genes);

    print "$cond\t$up_genes\t$down_genes\n";
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}


