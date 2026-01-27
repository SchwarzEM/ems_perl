#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $json_stem    = q{};
my $infile_list1 = q{};
my $infile_list2 = q{};
my $i            = 0;
my $j            = 0;

$json_stem    = $ARGV[0] if $ARGV[0];
$infile_list1 = $ARGV[1] if $ARGV[1];
$infile_list2 = $ARGV[2] if $ARGV[2];
$i            = $ARGV[3] if $ARGV[3];

my @targets = ();

if ( (! $json_stem ) or (! $infile_list1 ) or (! $infile_list2 ) ) {
    die "Format: make_jsons_24nov2025a.pl [json_stem] [single ligand file list]  [constant multiple target file list] [optional starting positive integer] ;\n";
}

if ( (! looks_like_number($i) ) or ( $i < 1 ) or ( $i != int($i) ) ) {
    warn "Format: make_jsons_24nov2025a.pl [json_stem] [single ligand file list]  [constant multiple target file list] [optional starting positive integer] ;\n";
    die "The [optional starting integer] must be (1) positive and (2) an integer, >= 1.\n";
}

# Once $i has been tested, reduce it by 1 so that loop numbering will work as if it started at 0.
$i--;

open my $LIST2, '<', $infile_list2;
while ( my $infile2 = <$LIST2> ) {
    chomp $infile2;
    push @targets, $infile2;
}
close $LIST2;

open my $LIST1, '<', $infile_list1;
while ( my $infile = <$LIST1> ) { 
    chomp $infile;
    $i++;
    $j = sprintf "%03i", $i;
    my $outfile = $json_stem. '_' . "$j.make_json.sh";
    $outfile    = safename($outfile);
    open my $OUTFILE, '>', $outfile;
    print $OUTFILE '$PROJECT/ems_perl/Hco_Acey/afserv_json_22nov2025.pl';
    print $OUTFILE " $json_stem" . '_' . "$j";
    print $OUTFILE " $infile";
    print $OUTFILE " @targets ;";
    print $OUTFILE "\n";
    close $OUTFILE;
}
close $LIST1;

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

