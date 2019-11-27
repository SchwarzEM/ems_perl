#!/usr/bin/env perl

# best_salmon_isoforms.pl -- Erich Schwarz <ems394@cornell.edu>, 11/25/2019. 
# Purpose: given salmon transcript results () and a CDS-to-gene table, pick best isoforms (ones with most TPMs in a data set).

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles  = ();
my $cds2gene = q{};
my $help;

my $header = "Gene\tGene_TPM\tCDS_1\tCDS_1_TPM\tCDS_2\tCDS_2_TPM\tCDS_3\tCDS_3_TPM\n";

my $data_ref;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'cds2gene=s'   => \$cds2gene,
             'help'         => \$help,   );

if ( $help or (! @infiles) or (! $cds2gene) ) { 
    die "Format: best_salmon_isoforms.pl\n",
        "    --infile|-i     <input files>\n",
        "    --cds2gene|-c   [CDS-to-gene table]\n",
        "    --help|-h       [print this message]\n",
        ;
}

open my $CDS2GENE, '<', $cds2gene;
while (my $input = <$CDS2GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $cds  = $1;
        my $gene = $2;
        $data_ref->{'CDS'}->{$cds}->{'gene'} = $gene;
    }
    else {
        die "In CDS-to-gene table $cds2gene, cannot parse: $input\n";
    }
}
close $CDS2GENE;

foreach my $infile (@infiles) { 
    open my $INFILE, '<', $infile;
    while (my $input = <$INFILE>) {
        chomp $input;
        if ( ( $input !~ /\A Name /xms ) and ( $input =~ /\A (\S+) \t \S+ \t \S+ \t (\S+) \t \S+ \z/xms ) ) {
            my $cds = $1;
            my $tpm = $2;
            if ( (! looks_like_number($tpm) ) or ( $tpm < 0 ) ) {
                die "In Salmon TPM data $infile, no non-negative, numerical TPM: $input\n";
            }
            my $gene = $data_ref->{'CDS'}->{$cds}->{'gene'};

            $data_ref->{'gene'}->{$gene}->{'CDS'}->{$cds}->{'tpm'} += $tpm;
            $data_ref->{'gene'}->{$gene}->{'tpm'} += $tpm;
        } 
        elsif ( $input !~ /\A Name /xms ) {
            die "In Salmon TPM data $infile, cannot parse: $input\n";
        }
    }
    close $INFILE;
}

my @genes = ();
@genes = sort keys %{ $data_ref->{'gene'} };

foreach my $gene (@genes) {
     print $header if $header;
     $header = q{};

     my $gene_tpm = $data_ref->{'gene'}->{$gene}->{'tpm'};

     my @cdses = sort { $data_ref->{'gene'}->{$gene}->{'CDS'}->{$b}->{'tpm'} <=> 
                        $data_ref->{'gene'}->{$gene}->{'CDS'}->{$a}->{'tpm'} }
                        keys %{ $data_ref->{'gene'}->{$gene}->{'CDS'} };

     print "$gene\t$gene_tpm";

     foreach my $i (0..2) {
         my $output = q{};
         my $cds    = q{};
         my $tpm    = q{};

         if ( exists $cdses[$i]) {
             $cds = $cdses[$i];
             if ( exists $data_ref->{'gene'}->{$gene}->{'CDS'}->{$cds}->{'tpm'} ) {
                 $tpm = $data_ref->{'gene'}->{$gene}->{'CDS'}->{$cds}->{'tpm'};
             }
         }

         $output = "$cds\t$tpm";
         print "\t$output";
     }
     print "\n";
}

