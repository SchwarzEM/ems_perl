#!/usr/bin/perl
#
# Mario Stanke
my $usage = "$0 -- filter a multiple fasta file to keep only those sequences that have a length between min and max\n";
$usage .= "Usage: $0 refseq.fa min max\n";
$usage .= "\n";
 
if ($#ARGV != 2) {
    die "Unknown option\n\n$usage";
}
 
my $multiseqfilename = $ARGV[0];
my $min = $ARGV[1];
my $max = $ARGV[2];
my $firstseq = 1;

#
# Go through the multiseq.fa file and keep only those sequences that are longer than $reflen.
#

$/=">";
open (MULTI, "<$multiseqfilename") || die "Couldn't open $multiseqfilename\n";
open (LONG, ">${multiseqfilename}.long") || die "Couldn't open ${multiseqfilename}.long\n";
while (<MULTI>) {
    s/>$//;
    /(.*)\n/;
    $seqname=$1;
    $sequencepart=$'; #'
    $sequence = $sequencepart;
    $sequence =~ s/\n//g;
    if (length $sequence >= $min && length $sequence <= $max){
	print LONG ">$seqname\n$sequencepart";
	$firstseq=0;
    }
}
close (LONG);
close (MULTI);

if ($firstseq==1) { # have lo longer sequences at all
    system("rm ${multiseqfilename}.long");
}
