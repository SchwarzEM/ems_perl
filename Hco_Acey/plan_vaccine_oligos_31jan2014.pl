#!/usr/bin/env perl

# plan_vaccine_oligos_31jan2014.pl -- Erich Schwarz <ems394@cornell.edu>, 1/31/2014.
# Purpose: given CDSes going from start to stop codon, print out a report of exactly which sequences I want for which oligos.

use strict;
use warnings;

my $handle_L      = 'GAACCT';
my $kpn1_kozak    = 'GGTACCATGG';
my $start_oligo_L = $handle_L . $kpn1_kozak;

my $handle_R      = 'GACTCA';
my $xba1_stop     = 'TCTAGA';
my $start_oligo_R = $handle_R . $xba1_stop ;

my $seqname = q{};

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        $seqname = $1;
    }
    elsif ($seqname) {
        $input =~ s/\s//g;
        if ( $input =~ / [^acgtACGT] /xms ) { 
            die "From sequence $seqname, can't parse: $input\n";
        }
        $data_ref->{'seqname'}->{$seqname}->{'seqtext'} .= $input;
    }
}

my @seqnames = sort keys %{ $data_ref->{'seqname'} };
foreach my $seqname1 (@seqnames) {
    my $seqid = $seqname1;
    $seqid =~ s/\.t\d+\z//;

    my $seqtext1 = $data_ref->{'seqname'}->{$seqname1}->{'seqtext'};

    my $gene_spec_oligo_L = substr($seqtext1, 4, 20);
    my $gene_spec_oligo_R = substr($seqtext1, -24, 19);
    $gene_spec_oligo_R    = revcomp($gene_spec_oligo_R);

    my $full_oligo_L = $start_oligo_L . $gene_spec_oligo_L;
    my $full_oligo_R = $start_oligo_R . $gene_spec_oligo_R;

    print q{*** }, $seqid, q{== ***}, "\n";
    print "\n";
    print $seqname1, q{L oligo: [constant 5' end with KpnI/Kozak -- then, 20 nt: 5'- 2nt - 6x 3-nt codons -3']}, "\n";
    print "\n";
    print q{    5'-}, $start_oligo_L, q{-3' then 5'-}, $gene_spec_oligo_L, qq{-3' i.e.\n};
    print q{    5'-}, $full_oligo_L, qq{-3'\n};
    print "\n";
    print $seqname1, q{R oligo: [constant 5' end with XbaI/stop -- then, 19 nt, *antisense*: 5'- 1nt - 6x 3-nt codons -3']}, "\n";
    print "\n";
    print q{    5'-}, $start_oligo_R, q{-3' then 5'-}, $gene_spec_oligo_R, qq{-3' i.e.\n};
    print q{    5'-}, $full_oligo_R, qq{-3'\n};
    print "\n";
}

sub revcomp { 
    my $in_string = $_[0];
    $in_string =~ tr/[acgtACGT]/[tgcaTGCA]/;
    my @in_residues = split //, $in_string;
    @in_residues = reverse @in_residues;
    my $out_string = join q{}, @in_residues;
    return $out_string;
}

