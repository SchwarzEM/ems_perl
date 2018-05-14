#!/usr/bin/env perl

# rawread_table.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/27/2011.
# Purpose: given three Ali files, make a single table of raw read counts for a given RNA-seq type.

use strict;
use warnings;
use Getopt::Long;

my $uniqs   = q{};
my $splices = q{};
my $multis  = q{};
my $lengths = q{};
my $help;

GetOptions ( 'uniqs=s'   => \$uniqs,
             'splices=s' => \$splices,
             'multis=s'  => \$multis,
             'lengths=s' => \$lengths,
             'help'      => \$help, );

if ($help or (! $uniqs ) or (! $splices ) or (! $multis ) ) { 
    die "Format: rawread_table.pl",
        " --uniqs|-u [unique recount]",
        " --splices|-s [splices count]",
        " --multis|-m [multis count]",
        " --lengths|-l [gene lengths]",
        " --help|-h [this message]\n";
}

my $wbgene    = q{};
my $readcount = q{};
my $length    = q{};

my @headers = qw(       all.reads 
                       uniq.reads 
                     splice.reads 
                  non_multi.reads 
                      multi.reads 
                      gene_length );
my $data_ref;

open my $UNIQS, '<', $uniqs or die "Can't open uniq readcount file $uniqs: $!";
while (my $input = <$UNIQS>) { 
    chomp $input;
    if ( $input !~ /\A .+ \s WBGene\d+ \s .* \b \d+ \s* \z /xms ) { 
        warn "Can't parse input line from uniq readcount file $uniqs: $input\n";
    }
    if ( $input =~ /\A .+ \s (WBGene\d+) \s .* \b (\d+) \s* \z/xms ) {
        $wbgene    = $1;
        $readcount = $2;
        if (exists $data_ref->{$wbgene}->{'uniq.reads'} ) { 
            die "Two different unique readcounts for the single gene $wbgene!\n";
        }
        $data_ref->{$wbgene}->{'uniq.reads'}      = $readcount;
        $data_ref->{$wbgene}->{'non_multi.reads'} = $readcount;
        $data_ref->{$wbgene}->{'all.reads'}       = $readcount;
    }
}
close $UNIQS or die "Can't close filehandle to uniq readcount file $uniqs: $!";

open my $SPLICES, '<', $splices or die "Can't open splice readcount file $splices: $!";
while (my $input = <$SPLICES>) {
    chomp $input;
    if ( $input !~ /\A .+ \s WBGene\d+ \s .* \b \d+ \s* \z/xms ) {
        warn "Can't parse input line from splice readcount file $splices: $input\n";
    }
    if ( $input =~ /\A .+ \s (WBGene\d+) \s .* \b (\d+) \s* \z/xms ) {
        $wbgene    = $1;
        $readcount = $2;
        $data_ref->{$wbgene}->{'splice.reads'}     = $readcount;
        $data_ref->{$wbgene}->{'non_multi.reads'} += $readcount;
        $data_ref->{$wbgene}->{'all.reads'}       += $readcount;
    }
}
close $SPLICES or die "Can't close filehandle to splice readcount file $splices: $!";

open my $MULTIS, '<', $multis or die "Can't open multi readcount file $multis: $!";
while (my $input = <$MULTIS>) {
    chomp $input;
    if ( $input !~ /\A .+ \s WBGene\d+ \s .* \b \d+ \s* \z/xms ) {
        warn "Can't parse input line from multis readcount file $multis: $input\n";
    }
    if ( $input =~ /\A .+ \s (WBGene\d+) \s .* \b (\d+) \s* \z/xms ) {
        $wbgene    = $1;
        $readcount = $2;
        $data_ref->{$wbgene}->{'multi.reads'}  = $readcount;
        $data_ref->{$wbgene}->{'all.reads'}   += $readcount;
    }
}
close $MULTIS or die "Can't close filehandle to multi readcount file $multis: $!";

open my $LENGTHS, '<', $lengths or die "Can't open gene length file $lengths: $!";
while (my $input = <$LENGTHS>) {
    chomp $input;
    if ( $input !~ /\A WBGene\d+ \t \d+\.\d+ \t /xms ) {
        warn "Can't parse input line from gene length file $lengths: $input\n";
    }
    if ( $input =~ /\A (WBGene\d+) \t (\d+\.\d+) \t /xms ) {
        $wbgene = $1;
        $length = $2;
        if (exists $data_ref->{$wbgene}->{'gene_length'} ) {
            die "Two different lengths for the single gene $wbgene!\n";
        }
        $data_ref->{$wbgene}->{'gene_length'} = $length;
    }
}
close $LENGTHS or die "Can't close filehandle to gene length file $lengths: $!";

my $top_line = join "\t", @headers;
$top_line    = "\t" . $top_line . "\n";
print $top_line;

foreach my $wbgene1 ( sort keys %{ $data_ref } ) { 
    my @reads = ();
    foreach my $header (@headers) { 
        push @reads, $data_ref->{$wbgene1}->{$header};
    }
    my $read_line = join "\t", @reads;
    print "$wbgene1\t$read_line\n";
}

