#!/usr/bin/env perl

# filter_minimum_ACGT_fastq.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/5/2011.
# Purpose:  Take four lines at a time.  Print only those reads with a minimum number of continuous ACGT (or acgt) in them.

use strict;
use warnings;
use Getopt::Long;

my $input  = q{};

my $head1  = q{};
my $nt_seq = q{};
my $head2  = q{};
my $quals  = q{};

my %opts   = ();

GetOptions(
    "input=s"     => \$opts{'input'},
    "output=s"    => \$opts{'output'},
    "minimum:i"   => \$opts{'minimum'}, 
    "help"        => \$opts{'help'},
);

if (    (! $opts{'input'}   ) 
     or (! $opts{'output'}  ) 
     or (! $opts{'minimum'} ) ) {
    $opts{'help'} = 1;
}

if ( $opts{'help'} ) {
    print "\n";
    print "usage: filter_minimum_ACGT_fastq.pl\n";
    print "       -i|--input      <in>       input file name (fastq); opt. stream input via '-'\n";
    print "       -o|--output     <out>      output file name (fastq), required\n";
    print "       -m|--minimum    <integer>  Minimum length of continuous ACGT in the read\n";
    print "\n";
    print "       -h|--help                  help - print this message\n";
    print "\n";
    exit;
}


# Accept either a stream from '-' or a standard file.
my $INFILE;
if ($opts{'input'} eq '-') {
    # Special case: get the stdin handle
    $INFILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INFILE,  '<', $opts{'input'}  or die "Cannot open input file $opts{'input'}: $!";
}

open my $OUTFILE, '>', $opts{'output'} or die "Cannot open output file $opts{'output'}: $!";

# If no positive-integer minimum length given, make it zero:
if ( (! $opts{'minimum'} ) or ( $opts{'minimum'} < 1 ) ) {
    die "Not able to accept minimum length of $opts{'minimum'}!\n";
}

while (<$INFILE>) { 
    # Tricky Perl syntax: 
    #     it forces me to start with $_ already loaded from <$INFILE>,
    #     but afterwards needs prompting for linefeeds.
    $head1 = $_;
    chomp $head1;
    $nt_seq = <$INFILE>;
    chomp $nt_seq;
    $head2  = <$INFILE>;
    chomp $head2;
    $quals  = <$INFILE>;
    chomp $quals;

    if ( ( $head1 !~ /\A @ \S+ /xms ) or ( $head2 !~ /\A \+ /xms ) ) { 
        warn "Can't parse one or both of these headers:\n";
        warn "$head1\n";
        warn "$head2\n";
        die;
    }

    if ( $nt_seq =~ / [ACGTacgt]{$opts{'minimum'},} /xms ) { 
        print $OUTFILE "$head1\n";
        print $OUTFILE "$nt_seq\n";
        print $OUTFILE "$head2\n";
        print $OUTFILE "$quals\n";
    }
}

