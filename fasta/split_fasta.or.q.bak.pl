#!/usr/bin/env perl

# split_fasta.or.q.pl -- Erich Schwarz <emsch@caltech.edu>, 10/3/2012.
# Purpose: given an integer N and a FASTA/Q file, make N subfiles with ~equal numbers of entries; makes huge jumbled FASTA/Q readfiles easier to sort by divide-and-conquer.

use strict;
use warnings;
use Getopt::Long;

my $fasta;
my $fastq;
my $infile   = q{};
my $splitval  = 0;
my $splitfrac = 1;
my $lines_infile = 0;
my $help;

GetOptions ( 'infile=s'   => \$infile,
             'splitval=i' => \$splitval,
             'f'          => \$fasta,
             'q'          => \$fastq,
             'help'       => \$help, );

&die_loudly if $help;

if ( (! $fasta) and (! $fastq) ) { 
    warn "Must specify either FASTA or FASTQ input with -f or -q\n\n";
    &die_loudly;
}

if ( ( $splitval != int($splitval) ) or ( $splitval <= 1 ) ) { 
    warn "Split value must be an integer >=2, not \"$splitval\"\n\n";
    &die_loudly;
}

if ( $splitval > 99 ) {
    warn "Split value $splitval -- over 99? are you high?\n\n";
    &die_loudly;
}

if ( (! $infile) or (! $splitval) or ( $fasta and $fastq ) or ( (! $fasta) and (! $fastq) ) ) { 
    &die_loudly;
}

$splitfrac = 2 if $fasta;
$splitfrac = 4 if $fastq;

open my $INFILE, '<', $infile or die "Can't open input file $infile\n";
while (my $input = <$INFILE>) { 
    $lines_infile++;
}
close $INFILE or die "Can't close filehandle to input file $infile\n";

my $stanzas_infile = ($lines_infile / $splitfrac);
$stanzas_infile = int $stanzas_infile;

my $stanzas_outfile = ($stanzas_infile / $splitval);
$stanzas_outfile = int $stanzas_outfile;

if ( ( $stanzas_outfile <= 0 ) or ( $stanzas_outfile == $stanzas_infile ) ) { 
    die "Cannot do effective split: $stanzas_infile input stanzas, $stanzas_outfile output stanzas\n";
}

my $outfile_number    = 1;
my $outfile_suffix    = sprintf "%02i", $outfile_number;
my $outfile           = $infile . q{.out.} . $outfile_suffix;
$outfile              = safename($outfile);
my $printable_stanzas = 0;

open $INFILE, '<', $infile or die "Can't open input file $infile\n";
open my $OUTFILE, '>', $outfile or die "Can't open first output subset file $outfile\n";
while (my $input1 = <$INFILE>) {
    my $input2 = <$INFILE>;
    my $input3;
    my $input4;
    if ($fastq) { 
        $input3 = <$INFILE>;
        $input4 = <$INFILE>;
    }
    if (                 ( ( $input1 !~ /\S/xms ) or ( $input2 !~ /\S/xms ) )
         or ( $fastq and ( ( $input3 !~ /\S/xms ) or ( $input4 !~ /\S/xms ) ) ) ) {  
        die "Anomalous stanza:\n$input1$input2$input3$input4\n";
    }
    $printable_stanzas++;
    if ( $printable_stanzas > $stanzas_outfile ) { 
        $outfile_number++;
        # Note that this limit may lead to unevenly-sized outputs, but being able to enforce an output quantity is more predictable and useful.
        if ( $outfile_number <= $splitval ) { 
            close $OUTFILE;
            $outfile_suffix    = sprintf "%02i", $outfile_number;
            $outfile           = $infile . q{.out.} . $outfile_suffix;
            $outfile           = safename($outfile);
            $printable_stanzas = 1;
            open $OUTFILE, '>', $outfile or die "Can't open next output subset file $outfile\n";
        }
    }
    print $OUTFILE $input1;
    print $OUTFILE $input2;
    print $OUTFILE $input3 if $fastq;
    print $OUTFILE $input4 if $fastq;
}
close $INFILE or die "Can't close filehandle to input file $infile\n";
close $OUTFILE or die "Can't close filehandle to final output subset file $outfile\n";

sub die_loudly {
    die "Format: split_fastq.pl\n",
        "    -f [FASTA input, default]\n",
        "    -q [FASTQ input; mutually exclusive with -f]\n",
        "    --infile|-i [single FASTA/Q input file]\n",
        "    --splitval|-s [number of partitions to make; 2-99]\n",
        "    --help|-h\n",
        ;
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

