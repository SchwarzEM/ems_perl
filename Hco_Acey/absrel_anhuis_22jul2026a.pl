#!/usr/bin/env perl

use strict;
use autodie;
use warnings;

my $data_ref;

$data_ref->{'taxon_type'}->{'Baylor'} = 'Anhui';
$data_ref->{'taxon_type'}->{'Keiser'} = 'Anhui';

$data_ref->{'taxon_type'}->{'Aroian'}  = 'non_Anhui';
$data_ref->{'taxon_type'}->{'Ilik2'}    = 'non_Anhui';
$data_ref->{'taxon_type'}->{'Mag3'}     = 'non_Anhui';
$data_ref->{'taxon_type'}->{'Oita'}     = 'non_Anhui';
$data_ref->{'taxon_type'}->{'obscurus'} = 'non_Anhui';

my $header = "Anhui_max_p\tNon_Anhui_max_p\tAnhui_min_p\tNon_Anhui_min_p";

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t [^\t]+ \t [^\t]+ \t [^\t]+ \t (\S+) \t (\S+) \z/xms ) {
        my $pangene = $1;
        my $taxon   = $2;
        my $orig_p  = $3;
        my $adj_p   = $4;

        if ( exists $data_ref->{'taxon_type'}->{$taxon} ) {
            $data_ref->{'pangene'}->{$pangene}->{'taxon'}->{$taxon}->{'orig_p'} = $orig_p;
            $data_ref->{'pangene'}->{$pangene}->{'taxon'}->{$taxon}->{'adj_p'}  = $adj_p;

            my $taxon_type = $data_ref->{'taxon_type'}->{$taxon};

            if (! $data_ref->{'pangene'}->{$pangene}->{'max_p'}->{$taxon_type} ) {
                $data_ref->{'pangene'}->{$pangene}->{'max_p'}->{$taxon_type} = $adj_p;
            }
            elsif ( $adj_p > $data_ref->{'pangene'}->{$pangene}->{'max_p'}->{$taxon_type} ) {
                $data_ref->{'pangene'}->{$pangene}->{'max_p'}->{$taxon_type} = $adj_p;
            }

            if (! $data_ref->{'pangene'}->{$pangene}->{'min_p'}->{$taxon_type} ) {
                $data_ref->{'pangene'}->{$pangene}->{'min_p'}->{$taxon_type} = $adj_p;
            }
            elsif ( $data_ref->{'pangene'}->{$pangene}->{'min_p'}->{$taxon_type} > $adj_p ) {
                $data_ref->{'pangene'}->{$pangene}->{'min_p'}->{$taxon_type} = $adj_p;
            }
        }
    }
    else {
        die "Cannot parse input: $input\n";
    } 
}

my @pangenes = sort keys %{ $data_ref->{'pangene'} };
foreach my $pangene (@pangenes) {
    print "$header\n" if $header;
    $header = q{};

    my $Anhui_max_p     = 'NA';
    my $Anhui_min_p     = 'NA';
    my $non_Anhui_max_p = 'NA';
    my $non_Anhui_min_p = 'NA';

    if ( exists $data_ref->{'pangene'}->{$pangene}->{'max_p'}->{'Anhui'} ) {
        $Anhui_max_p = $data_ref->{'pangene'}->{$pangene}->{'max_p'}->{'Anhui'};
    }
    if ( exists $data_ref->{'pangene'}->{$pangene}->{'min_p'}->{'Anhui'} ) {
        $Anhui_min_p = $data_ref->{'pangene'}->{$pangene}->{'min_p'}->{'Anhui'};
    }

    if ( exists $data_ref->{'pangene'}->{$pangene}->{'max_p'}->{'non_Anhui'} ) {
        $non_Anhui_max_p = $data_ref->{'pangene'}->{$pangene}->{'max_p'}->{'non_Anhui'};
    }
    if ( exists $data_ref->{'pangene'}->{$pangene}->{'min_p'}->{'non_Anhui'} ) {
        $non_Anhui_min_p = $data_ref->{'pangene'}->{$pangene}->{'min_p'}->{'non_Anhui'};
    }
    print "$Anhui_max_p\t$non_Anhui_max_p\t$Anhui_min_p\t$non_Anhui_min_p\n";

}
