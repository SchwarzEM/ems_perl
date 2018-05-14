#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my @infile_names = ();
my $go_data_file = q{};
my $help;

GetOptions(
    "input=s{,}" => \@infile_names,
    "go=s"       => \$go_data_file,
    "help"       => \$help,
);


# Take first argument as source of file names; do streaming edit of rest.

if ($help or (! @infile_names) or (! $go_data_file) ) { 
    die 'Format: make_job_func_wilcox_05dec2012.pl --input|-i [list of target files, either in 1+ files or \'-\' stream] --go|-g [/PATH/go_DATE-termdb-tables dir] --help|-h [print this message]', "\n", ;
}

my $header = '#!/bin/bash' . "\n\n";

my $INFILE;
foreach my $infile (@infile_names) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INFILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INFILE,  '<', $infile or die "Cannot open input file $infile: $!";
    }
    while (my $input = <$INFILE>) { 

        chomp $input;
        if ( $input =~ /\A (\S+) \.func_input\.txt \z/xms ) { 
            my $stem       = $1;
            my $output_dir = $stem . '.func_wilcoxon_outdir';
            $output_dir    = safename($output_dir);

            # Print header -- only once, and only if there's something else to print:
            print $header if $header;
            $header = q{};

            # Then print one stanza per $input:
            print "    mkdir $output_dir ;\n";
            print "    nohup func_wilcoxon -t $go_data_file -i $input -o $output_dir &\n";
            print "\n";
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }
    close $INFILE or die "Cannot close filehandle to input file $infile: $!";
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

