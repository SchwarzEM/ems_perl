#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $names = q{};
my $gff3  = q{};

$names = $ARGV[0] if $ARGV[0];
$gff3  = $ARGV[1] if $ARGV[1];

if ( (! $names ) or (! $gff3 ) ) {
    die "Format: rename_gff3_seqs.pl [old-to-new sequence name table] [GFF3] > [sequence-renamed GFF3]\n";
}

my %old2new = ();

open my $NAMES, '<', $names;
while ( my $input = <$NAMES> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $old = $1;
        my $new = $2;
        $old2new{$old} = $new;
    }
    else {
        die "From nametable $names, cannot parse: $input\n";
    }
}
close $NAMES;

open my $GFF3, '<', $gff3;
while ( my $input = <$GFF3> ) {
    chomp $input;
    if ( $input =~ /\A [#] /xms ) {
        print "$input\n";
    }
    elsif ( $input =~ /\A (\S+) \t (.*) \z/xms ) {
        my $old  = $1;
        my $text = $2;
        if ( exists $old2new{$old} ) {
            my $new = $old2new{$old};
            print "$new\t$text\n";
        }
        else {
            die "From GFF3 $gff3, cannot identify new sequence name in: $input\n";
        }
    }
    else {
        die "From GFF3 $gff3, cannot parse: $input\n";
    }
}
close $GFF3;

