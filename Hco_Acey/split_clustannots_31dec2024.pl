#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};
$infile    = $ARGV[0] if $ARGV[0];

my $data_ref;

my $head_text = "Gene\tClust_annot";

if (! $infile ) {
    die "Format: split_clustannots_31dec2024.pl [annotations for multiple cluster memberships] > (set of single-cluster annotation files)\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    if ( ( $input !~ /\A Gene /xms ) and ( $input =~ /\A (\S+) \t (\S [^\t]+ \S) \z/xms ) ) {
        my $gene    = $1;
        my $cluster = $2;
        $data_ref->{'cluster'}->{$cluster}->{'gene'}->{$gene} = 1;
    }
}
close $INFILE;

my @clusters = sort keys %{ $data_ref->{'cluster'} };

foreach my $cluster (@clusters) {
    my $clustname = $cluster;
    $clustname    =~ s/\s/_/g;

    my @genes = sort keys %{ $data_ref->{'cluster'}->{$cluster}->{'gene'} };

    my $outfile = "$clustname.$infile";
    $outfile    = safename($outfile);

    my $header = $head_text;

    open my $OUTFILE, '>', $outfile;
    foreach my $gene (@genes) {
        print $OUTFILE "$header\n" if $header;
        $header = q{};
        print $OUTFILE "$gene\t$cluster\n";
    }
    close $OUTFILE;
}

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

