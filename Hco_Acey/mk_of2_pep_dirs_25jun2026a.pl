#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $list   = q{};
my $prefix = q{};
my $i      = 0;

$list   = $ARGV[0] if $ARGV[0];
$prefix = $ARGV[1] if $ARGV[1];

if ( (! $list ) or (! $prefix ) ) {
    die "Format: mk_of2_pep_dirs_25jun2026a.pl [list of FASTAs] [new dir prefix] => (mdkir and copy FASTAs to a numbered series of dirs)\n";
}

open my $LIST, '<', $list;
while ( my $input = <$LIST> ) {
    chomp $input;
    my @files = split '; ', $input;
    $i++;
    my $outdir = "$prefix.$i";
    $outdir    = safename($outdir);
    print "mkdir $outdir ;\n";
    foreach my $file (@files) {
        if (! -e $file ) {
            die "Cannot find file \"$file\" in: $input\n";
        }
        print "cp -ip $file $outdir ;\n";
    }
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

