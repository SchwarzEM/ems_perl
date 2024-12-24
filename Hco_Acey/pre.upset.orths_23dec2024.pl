#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile  = q{};
my @sp_list = qw( s_hermaphroditum s_carpocapsae c_elegans h_bacteriophora_Hawdon a_suum t_spiralis );
my %targets = ();

$infile = $ARGV[0] if $ARGV[0];

if (! $infile ) {
    die "Format: pre.upset.orths_23dec2024.pl [orthofinder file] > [pre-upset data table]\n";
}

foreach my $species ( @sp_list ) {
    $targets{$species} = 1;
}

my $header = join "\t", @sp_list;
$header    = "Name\t$header";

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    my $ogroup = q{};
    my $desc   = q{};
    my $output = q{};
    if ( $input =~ /\A (OG\d+) .+? [:] \s+ (.+) \z/xms ) { 
        $ogroup = $1;
        $desc   = $2;
        $output = $ogroup;
        foreach my $species ( @sp_list ) {
            if ( $desc =~ /$species/xms ) {
                $output = $output . "\t1";
            }
            else {
                $output = $output . "\t0";
            }
        }
        if ($header) {
            print "$header\n";
            $header = q{};
        }
        print "$output\n";
    }
    else {
        die "From infile $infile, cannot parse: $input\n";
    }
}
close $INFILE;
