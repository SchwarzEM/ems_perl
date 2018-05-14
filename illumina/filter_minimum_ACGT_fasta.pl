#!/usr/bin/env perl

# filter_minimum_ACGT_fasta.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/6/2011.
# Purpose:  Take two lines at a time.  Print only those reads with a minimum number of continuous ACGT (or acgt) in them.

use strict;
use warnings;
use Getopt::Long;

my $input  = q{};

my $head1  = q{};
my $nt_seq = q{};

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
    print "usage: filter_minimum_ACGT_fasta.pl\n";
    print "       -i|--input      <in>       input file name (FASTA one-line read); opt. stream input via '-'\n";
    print "       -o|--output     <out>      output file name (FASTA one-line read), required\n";
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

    if ( $head1 !~ /\A > \S+ /xms ) { 
        die "Can't parse this header: $head1\n";
    }

    if ( $nt_seq =~ / [ACGTacgt]{$opts{'minimum'},} /xms ) { 
        print $OUTFILE "$head1\n";
        print $OUTFILE "$nt_seq\n";
    }
}

