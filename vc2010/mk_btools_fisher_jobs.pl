#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $motif_bed_list = q{};
my $tr_bed         = q{};
my $genome_file    = q{};
my $tag            = q{};

$motif_bed_list = $ARGV[0] if $ARGV[0];
$tr_bed         = $ARGV[1] if $ARGV[1];
$genome_file    = $ARGV[2] if $ARGV[2];
$tag            = $ARGV[3] if $ARGV[3];

if ( (! $motif_bed_list ) or (! $tr_bed ) or (! $genome_file ) or (! $tag ) ) {
    die "Format: mk_btools_fisher_jobs.pl",
        " [list of FIMO motif hit BED files]",
        " [tandem repeat BED]",
        " [BEDtools genome file]",
        " [specifying tag, like \"2023.05.31.01\"]",
        " > [bedtools fisher line-commands]\n",
        ;
}

if ( (! -r $motif_bed_list ) or (! -r $tr_bed ) or (! $genome_file ) ) {
    die "Cannot read either $motif_bed_list or $tr_bed or $genome_file or combination of them.\n";
}

open my $LIST, '<', $motif_bed_list;
while ( my $motif_bed = <$LIST> ) {
    chomp $motif_bed;
    if (! -r $motif_bed ) {
        die "Cannot read FIMO motif hit BED file: $motif_bed\n";
    }

    # For each motif BED, extract its stem name for use in the output and error files.
    my $stem = q{};
    my $motif_bed_basename = basename($motif_bed);
    if ( $motif_bed_basename =~ /\A (\S+)_fimo\.bed\z/xms ) {
        $stem = $1;
    }
    else {
        die "Cannot extract stem name from: $motif_bed (basename, $motif_bed_basename)\n";
    }

    my $output = "$stem.fisher.$tag.txt";
    my $error  = "$stem.fisher.$tag.err";

    $output = safename($output);
    $error  = safename($error);

    print "bedtools fisher -f 1.0 -a $motif_bed -b $tr_bed -g $genome_file 1>$output 2>$error ;\n";
}
close $LIST;


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

