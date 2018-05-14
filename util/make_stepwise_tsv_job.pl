#!/usr/bin/env perl

# make_stepwise_tsv_job.pl -- Erich Schwarz <ems394@cornell.edu>, 2/11/2012.

use strict;
use warnings;
use Getopt::Long;

my $start_file     = q{};
my @building_files = ();
my $target_file    = q{};
my $i              = 1;
my $help;

GetOptions ( 'start=s'    => \$start_file,
             'build=s{,}' => \@building_files,
             'target=s'   => \$target_file,
             'help'       => \$help,        );

if ($help or (! $start_file) or (! @building_files) or (! $target_file)  ) { 
    die "Format: make_stepwise_tsv_job.pl\n",
        "            -s|--start [starting file for add_tab_annots.pl chain]\n",
        "            -b|--go_assoc    [list of more building files, in order to add]\n",
        "            -t|--no_IEA      [final target file name]\n",
        "            -h|--help        [print this message]\n",
        "            [print script to <STDOUT>]\n",
        ;
}

$target_file = safename($target_file);

my $first_build_file = shift @building_files;

print '#!/bin/bash', "\n\n";

print "    add_tab_annots.pl -f -i $start_file -a $first_build_file > tmp$i.txt ;\n";
print "    fill_tsv_tabs.pl tmp$i.txt > ";
$i++;
print "tmp$i.txt ;\n\n";

foreach my $build_file (@building_files) { 
    print "    add_tab_annots.pl -f -i tmp$i.txt -a $build_file > ";
    $i++;
    print "tmp$i.txt ;\n";

    print "    fill_tsv_tabs.pl tmp$i.txt > ";
    $i++;
    print "tmp$i.txt ;\n\n";
}

print "    mv tmp$i.txt $target_file ;\n";
print "\n";

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


