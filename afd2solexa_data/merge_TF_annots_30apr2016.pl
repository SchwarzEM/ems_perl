#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header  = "Gene\tTF";

my $data_ref;

foreach my $infile (@ARGV) {
    if ( $infile =~ /\A (\S+)_TFs\.WS\d+\.txt \z/xms ) { 
        my $source  = $1;
        open my $INFILE, '<', $infile;
        while (my $input = <$INFILE>) {
            chomp $input;
            if ( $input =~ /(WBGene\d+\S+)/xms ) {
                my $gene = $1;
                $data_ref->{'gene'}->{$gene}->{'source'}->{$source} = 1;
            }
            else { 
                die "From input file $infile, cannot parse: $input\n";
            }
        }
        close $INFILE;
    }
    else {
        die "Cannot parse input file name $infile\n";
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my @sources     = sort keys %{ $data_ref->{'gene'}->{$gene}->{'source'} };
    my $source_text = join q{+}, @sources;
    my $annot       = q{TF[} . $source_text . q{]};

    print "$header\n" if $header;
    $header = q{};
    print "$gene\t$annot\n";
}

