#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @infiles = ();
@infiles = @ARGV if @ARGV;

my $data_ref;

if (! @infiles) {
    die "Format: pool_oma_acc2names_12nov2018.pl [most to least preferred acc2names files] > [prioritized sum]\n";
}

foreach my $infile (@infiles) {
    open my $INFILE, '<', $infile;
    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (.+) \z/xms ) { 
            my $acc  = $1;
            my $name = $2;
            if (! exists $data_ref->{'acc'}->{$acc}->{'name'} ) { 
                $data_ref->{'acc'}->{$acc}->{'name'} = $name;
            }
        }
        else {
            die "From input file $infile, cannot parse: $input\n";
        }
    }
    close $INFILE;
}

my @accs = sort keys %{ $data_ref->{'acc'} };
foreach my $acc (@accs) {
    my $name = $data_ref->{'acc'}->{$acc}->{'name'};
    print "$acc\t$name\n";
}

