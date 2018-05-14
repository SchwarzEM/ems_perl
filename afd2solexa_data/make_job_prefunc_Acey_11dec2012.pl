#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $infile_list = q{};
my $gene2go     = q{};
my $help;

GetOptions ( 'infile_list=s' => \$infile_list,
             'gene2go=s'     => \$gene2go,
             'help'         => \$help,     );

if ( $help or (! $infile_list) or (! $gene2go) ) { 
    die "Format: make_job_prefunc_Acey_11dec2012.pl\n",
        "    --infile_list|-i  [input file, or stream '-', with list of files]\n",
        "    --gene2go|-g      [gene2go table; e.g., \"NAME.blast2go.annot.filt.txt\"]\n",
        "    --help|-h         [print this message]\n",
        ;
}

# E.g., $gene2go = '../../augustus_preds_24oct2012/blast2go/Acey_2012.10.24.pep.max_isos.blast2go.annot.filt.txt';

if (! -e $gene2go) { 
    die "Can't find gene2go $gene2go\n";
}

my $INPUT_FILE;
if ($infile_list eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile_list or die "Can't open input file $infile_list. $!\n";
}
while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    if ( $input =~ /\A (\S+)\.txt \z/xms ) { 
        my $stem   = $1;
        my $output = $stem . '.func_input.txt';
        $output    = safename($output);
        print "    qual2func_table_25nov2012.pl -g $gene2go -t $input > $output ;\n";
    }
    else { 
        die "Can't parse input $input\n";
    }
}
close $INPUT_FILE or die "Can't close filehandle to input file $infile_list. $!\n";

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

