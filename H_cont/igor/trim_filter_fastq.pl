#!/usr/bin/env perl

# trim_filter_fastq.pl -- Igor Antoshechkin <igor.antoshechkin@caltech.edu>, 7/1/2010.
# Purpose: "converts qseq files to fasta or fastq, and can also calculate the average quality score for each read and trim a pre-determined number of bases from either end".

use strict;
use warnings;
use Getopt::Long;

my %opts = (
    trim3 => 0,
    trim5 => 0,
);
my $program_name = $0 =~ /([^\/]+)$/ ? $1 : '';

GetOptions(
    "help"        => \$opts{help},
    "input:s"     => \$opts{input},
    "output:s"    => \$opts{output},
    "fasta"       => \$opts{fasta},
    "fastq"       => \$opts{fastq},
    "qual"        => \$opts{qual},
    "trim3:s"     => \$opts{trim3},
    "trim5:s"     => \$opts{trim5},
    "meanqual:s"  => \$opts{meanqual},
    "threshold:s" => \$opts{threshold},
);

if (   !$opts{input}
    || !$opts{output}
    || !( $opts{fasta} || $opts{fastq} || $opts{qual} ) )
{
    $opts{help} = 1;
}

if ( $opts{help} ) {
    print
      "usage: $program_name -input <fastq> -output <out> -fasta|fastq|qual\n";
    print "       -help                     help - print this message\n";
    print "       -input     <in>           input file name (fastq)\n";
    print "       -output    <out>          output file name\n";
    print "       -fasta                    write fasta\n";
    print "       -fastq                    write fastq\n";
    print "       -qual                     write quality\n";
    print "       -trim3     <int>          trim bases from 3 prime end\n";
    print "       -trim5     <int>          trim bases from 5 prime end\n";
    print
      "       -meanqual  <float>        avarage quality cutoff: default 0\n";
    print
"       -threshold <float>        quality threshold - reads that have any scores below this value will be filtered out\n";
    exit;
}

open( OUT, ">$opts{output}" ) || die "cannot open $opts{output}: $!\n";

my $count     = 0;
my $filtered  = 0;
my $qc        = 0;
my $threshold = 0;
my $id        = '';
my ( $f1, $f2, $seq, $qualstr, @quality );

open( IN, "<$opts{input}" ) || die " cannot open $opts{input}: $!\n";
while (<IN>) {
    chomp;
    if (/^@(.+)/) {
        if ($id) {
            my $passed = 1;
            if ( $opts{meanqual} || $opts{threshold} ) {
                @quality = split( /|/, $qualstr );
                if ( $opts{meanqual} ) {
                    my $mean = 0;
                    foreach (@quality) {
                        $mean += ord($_) - 64;
                    }
                    $mean = $mean / ( $#quality + 1 );
                    if ( $mean < $opts{meanqual} ) {
                        $qc++;
                        $passed = 0;
                    }
                }
                if ( $opts{threshold} && $passed ) {
                    foreach (@quality) {
                        if ( ord($_) - 64 < $opts{threshold} ) {
                            $passed = 0;
                            $threshold++;
                            last;
                        }
                    }
                }
            }
            if ($passed) {
                if ( $opts{fasta} ) {
                    print OUT ">$id\n$seq\n";
                }
                if ( $opts{fastq} ) {
                    print OUT "\@$id\n$seq\n";
                    print OUT "+$id\n$qualstr\n";
                }
                if ( $opts{qual} ) {
                    print OUT ">$id\n";
                    print OUT join( " ", map { ord($_) - 64 } @quality ), "\n";
                }
            }
        }

        $id = $1;
        $f1 = 1;
        $count++;
        if ( $count % 1000000 == 0 ) {
            warn "$count reads processed\n";
        }
    }
    elsif ($f1) {
        $f1 = 0;
        $seq =
          substr( $_, $opts{trim5}, length($_) - $opts{trim5} - $opts{trim3} );
    }
    elsif (/^\+/) {
        $f2 = 1;
    }
    elsif ($f2) {
        $f2 = 0;
        $qualstr =
          substr( $_, $opts{trim5}, length($_) - $opts{trim5} - $opts{trim3} );
    }

}

if ($id) {
    my $passed = 1;
    if ( $opts{meanqual} || $opts{threshold} ) {
        @quality = split( /|/, $qualstr );
        if ( $opts{meanqual} ) {
            my $mean = 0;
            foreach (@quality) {
                $mean += ord($_) - 64;
            }
            $mean = $mean / ( $#quality + 1 );
            if ( $mean < $opts{meanqual} ) {
                $qc++;
                $passed = 0;
            }
        }
        if ( $opts{threshold} && $passed ) {
            foreach (@quality) {
                if ( ord($_) - 64 < $opts{threshold} ) {
                    $passed = 0;
                    $threshold++;
                    last;
                }
            }
        }
    }
    if ($passed) {
        if ( $opts{fasta} ) {
            print OUT ">$id\n$seq\n";
        }
        if ( $opts{fastq} ) {
            print OUT "\@$id\n$seq\n";
            print OUT "+$id\n$qualstr\n";
        }
        if ( $opts{qual} ) {
            print OUT ">$id\n";
            print OUT join( " ", map { ord($_) - 64 } @quality ), "\n";
        }
    }
}

warn "$count reads from processed\n";
warn "$qc reads did not pass mean quality cutoff of $opts{meanqual}\n"
  if $opts{meanqual};
warn
  "$threshold reads did not pass threshold quality cutoff of $opts{threshold}\n"
  if $opts{threshold};

