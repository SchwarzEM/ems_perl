#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use Scalar::Util qw(looks_like_number);

my $inlist = q{};
my $window = q{};
my $tmpdir = 'vcftools_tmpdir';

$inlist = $ARGV[0] if $ARGV[0];
$window = $ARGV[1] if $ARGV[1];
$tmpdir = $ARGV[2] if $ARGV[2];

$tmpdir = safename($tmpdir);

if ( (! $inlist ) or (! $window ) or (! looks_like_number($window) ) or ( $window != int($window) ) or ( $window < 1 ) ) {
    die "Format: make_snpdensities_13mar2026a.pl\n",
        "    [input *.vcf.gz file list]\n",
        "    [window size; positive integer]\n",
        "    [tmp dir; default name vcftools_tmpdir]\n",
        "    > [vcftools --gzvcf X --SNPdensity commands]\n",
        ;
} 

open my $INLIST, '<', $inlist;
while ( my $infile = <$INLIST> ) {
    chomp $infile;
    my $basename = basename($infile);

    # Sample infile:
    # Ilik2_gen.ngm.Anhui_illu.sens.2026.01.29.pic_flag.sort.Ilik2_chrI.fb.2026.03.04.01.vcf.gz
    if ( $basename =~ / \A (\S+_gen\.ngm\.\S+_illu) \. .+ \. (\S+?_chr[IVX]+) \. .+ \.vcf\.gz \z/xms ) {
        my $gen = $1;
        my $chr = $2;
        print "vcftools --gzvcf $infile --remove-indels --SNPdensity $window --out $gen.$chr --temp $tmpdir ;\n" ;
    }
    else {
        die "Can't parse basename: $basename\n";
    }
}
close $INLIST;

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

