#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir

my $start_dir = getcwd;
my $lib_no    = q{};
my $init_desc = q{};
my $lib_desc  = q{};
my $recheck1  = 0;

my $header      = '#!/bin/bash' . "\n\n";
my $output_file = q{};

while (my $input = <>) {
    chomp $input;

    # Sample three input lines:
    # 
    # Lane #1 : (17751) Index #13 Acey intestines 19dpi -DEX biorep 1
    #      https://jumpgate.caltech.edu/library/17751
    # https://jumpgate.caltech.edu/runfolders/volvox02/161005_SN787_0580_BH3WGFBCXY/Unaligned/Project_17751_index13/Sample_17751/

    if ( (! $lib_no) and (! $init_desc) and (! $lib_desc) and (! $recheck1) 
          and ( $input =~ /\A Lane [ ] [#]\d+ [ ] : [ ] \( (\d+) \) [ ] Index [ ] [#]\d+ [ ] (\S.*\S) \s* \z/xms ) ) { 
        $lib_no    = $1;
        $init_desc = $2;

        $lib_desc = $init_desc;
        $lib_desc =~ s/\s+/_/g;
        $lib_desc =~ s/non[-]intestines/non.intestines/g;
        $lib_desc =~ s/[-]/no/g;
    }
    elsif ( $lib_no and $init_desc and $lib_desc and (! $recheck1) 
            and ( $input =~ /\A [ ]+ https[:]\/\/jumpgate\.caltech\.edu\/library\/ (\d+) \s* \z/xms ) ) { 
        my $putative_lib_no1 = $1;
        if ( $putative_lib_no1 != $lib_no ) {
            die "Inconsistent library number in: $input\n";
        }
        $recheck1 = 1;
    }
    elsif ( $lib_no and $init_desc and $lib_desc and $recheck1 
            and ( $input =~ /\A https:\/\/jumpgate\.caltech\.edu
                            \/runfolders\/volvox02
                            \/ [a-zA-Z0-9]+ _ [a-zA-Z0-9]+ _ [a-zA-Z0-9]+ _ [a-zA-Z0-9]+
                            \/Unaligned\/Project_ (\d+) _index \d+ \/Sample_ (\d+) \/ \z/xms ) ) {
        my $putative_lib_no1 = $1;
        my $putative_lib_no2 = $2;
        if ( ( $putative_lib_no1 != $lib_no ) or ( $putative_lib_no2 != $lib_no ) ) {
            die "Inconsistent library number in: $input\n";
        }

        # If we get through all of these, we can generate a stanza of commands.
        my $target_url    = $input;
        my $target_subdir = $lib_no . q{_} . $lib_desc;

        $target_subdir = catdir($start_dir, $target_subdir);
        $target_subdir = safename($target_subdir);

        $output_file = "$lib_no" . '_logfile_09oct2016.txt';
        $output_file = safename($output_file);

        # print header only once, at the top of the output script
        print $header if $header;
        $header = q{};

        print "cd $start_dir ;\n";
        print "mkdir $target_subdir ;\n";
        print "cd $target_subdir ;\n";
        print "wget --user=gec --password=gecilluminadata --output-file $output_file ";
        print "--recursive --level=1 --no-parent --no-directories --no-check-certificate --accept .fastq.gz ";
        print "$target_url ;\n";
        print "\n";

        # zero out these values:
        $lib_no    = q{};
        $init_desc = q{};
        $lib_desc  = q{};
        $recheck1  = 0;
    }
    else { 
        die "Cannot parse input line: $input\n";
    }
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

