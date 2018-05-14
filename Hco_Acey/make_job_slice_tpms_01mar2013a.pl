#!/usr/bin/env perl

# make_job_slice_tpms_23feb2013.pl -- Erich Schwarz <ems394@cornell.edu>, 2/23/2013.
# Purpose: given a list of RSEM 1.2.0 gene results files, make a script that will efficiently slice them.

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my @input_files = ();
my $prefix      = q{};
my $slice       = q{};
my $help;

GetOptions ( 'input_files=s{,}' => \@input_files,
             'prefix=s'         => \$prefix,
             'slice=i'          => \$slice,
             'help'             => \$help, );

if ( $help or (! $slice) or (! @input_files) ) { 
    die "Format: make_job_slice_tpms_23feb2013.pl\n",
        "    --input|-i  [input files or '-']\n",
        "    --prefix|-p [prefix to append to outputs -- optional]\n",
        "    --slice|-s  [slice number to use]\n",
        "    --help|-h   [print this message]\n",
        ;
}

# Accept either a stream from '-' or a standard file.
my $INPUT;
foreach my $infile (@input_files) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT, '<', $infile or die "Can't open input file $infile. $!\n";
    }

    # Record the incoming FASTA data.
    while (my $input = <$INPUT>) {

        chomp $input;
        my $stem = basename $input;
        $stem =~ s/\.genes\.results//;
        $stem =~ s/\.modencode\z//;
        $stem =~ s/Alb_ele/Alb/;
    
        print "    /mnt/home/emsch/perl.svn/trunk/afd2solexa_data/extract_rsem_slices.pl -p -o $prefix",
              "$stem -i $input -s $slice ;\n",
              ;
    }
    close $INPUT or die "Can't close filehandle to input file $infile: $!\n";
}


