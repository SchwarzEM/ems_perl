#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $json_stem    = q{};
my $infile_list1 = q{};
my $infile_list2 = q{};

$json_stem    = $ARGV[0] if $ARGV[0];
$infile_list1 = $ARGV[1] if $ARGV[1];
$infile_list2 = $ARGV[2] if $ARGV[2];

my $i = 0;
my $j = 0;

my @targets = ();

if ( (! $json_stem ) or (! $infile_list1 ) or (! $infile_list2 ) ) {
    die "Format: make_jsons_24nov2025a.pl [json_stem] [varying single ligand file list] [constant multiple target file list] ;\n";
}

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
    $j = sprintf "%02i", $i;
    my $outfile = "$json_stem.$j.json";
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

