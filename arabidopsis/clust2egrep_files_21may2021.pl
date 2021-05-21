#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $clust = q{};
$clust    = $ARGV[0] if $ARGV[0];

my $i = 0;

my $annot = q{};
$annot = $ARGV[1] if $ARGV[1];

if (     (! $clust ) 
      or (! $annot ) 
   ) {
    die "Format: clust2egrep_files_21may2021.pl [clust gene table] [annotation file] => single-column egrep files, and subset annotation files\n";
}

open my $CLUST, '<', $clust;
while (my $input = <$CLUST>) {
    if ( $i == 0 ) {
        chomp $input;
        my @cols = split '\t', $input;
        $i = @cols;
    }
}
close $CLUST;

foreach my $j (1..$i) {
    # Get automatically well-formatted, ASCII-sortable index numbers.
    my $DIGITS    = length($i);
    my $sf_format = '%0' . $DIGITS . 'u';
    my $k         = sprintf($sf_format, $j) or die "Can't zero-pad block number $j\n";

    my $grep_column_file = "clust_grep.$k.txt";
    $grep_column_file = safename($grep_column_file);

    open my $GREP_COLUMN_FILE, '>', $grep_column_file;

    my @clust_lines = `cut -f $j $clust | tail --lines=+3`;

    print $GREP_COLUMN_FILE "^Gene\n";

    foreach my $clust_line (@clust_lines) {
        chomp $clust_line;
        if ( $clust_line =~ /\S+/xms ) {
            print $GREP_COLUMN_FILE "^$clust_line\n";
        }
    }

    close $GREP_COLUMN_FILE;

    my $annot_subset_file = "$annot.subset.$k.txt";
    $annot_subset_file = safename($annot_subset_file);

    system "egrep -f $grep_column_file $annot > $annot_subset_file";
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

