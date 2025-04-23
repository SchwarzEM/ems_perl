#!/usr/bin/env perl

# clustify_ofind_tsv.pl -- Erich Schwarz <ems394@cornell.edu>, 4/23/2025.
# Purpose: read protein-based OrthoFinder TSV output (Orthogroups.tsv) into gene-based but still Clust-readable form.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $input_ofnd   = q{};
my @ofnd_groups  = ();
my $cds2gene2spp = ();
my $data_ref;
my $help;

GetOptions ( 'input_ofnd=s'   => \$input_ofnd,
             'cds2gene2spp=s' => \$cds2gene2spp,
             'help'           => \$help,
);

if ( $help or (! -r $input_ofnd ) or (! -r $cds2gene2spp ) ) {
    die "Format: prot2gene_ofind.pl\n",
        "        --input_ofnd|-i     [output file from OrthoFinder]\n",
        "        --cds2gene2spp|-c   [CDS-to-gene-to-species tab-delimited table; each CDS must be unique]\n",
        "        --help|-h           [print this message]\n",
        ;
}

open my $CDS2GENE2SPP, '<', $cds2gene2spp;
while (my $input = <$CDS2GENE2SPP>) {
    chomp $input;
    if ($input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms) {
        my $cds_id  = $1;
        my $gene    = $2;
        my $species = $3;
     
        # Enforce unique CDS names throughout all the species sets being archived.
        # This will also enforce unique mappings of CDS to gene to species.
        #
        # Note: we're mapping CDSes to genes, *not* proteins to genes.
        #     A single protein *can* be produced by 2+ (very similar) genes!
        #     Mercifully, there really is at least one CDS per gene...

        # It turns out that OrthoFinder silently changes ':', '(', and ')' in CDS/protein names to '_';
        #     so this must also be done with CDS/protein names in cds2gene2species.
        $cds_id =~ s/[:]/_/g;
        $cds_id =~ s/[(]/_/g;
        $cds_id =~ s/[)]/_/g;

        # I am dealing with wonky file format in which I must see *either* a ',' *or* a '\n' at the end of each CDS name.
        # Making matters worse, sometimes a 'final' CDS will have one or more \t between it and the final \n!

        my $cds_id_1 = "$cds_id,";
        my $cds_id_2 = "$cds_id\t";
        my $cds_id_3 = "$cds_id";

        my $gene_1   = "$gene,";
        my $gene_2   = "$gene\t";
        my $gene_3   = "$gene";

        if ( ( exists $data_ref->{'CDS'}->{$cds_id_1} ) or ( exists $data_ref->{'CDS'}->{$cds_id_2} ) ) {
            die "Redundant CDS ID: $cds_id_1 or $cds_id_2\n";
        }

        $data_ref->{'CDS'}->{$cds_id_1}->{'gene'} = $gene_1;
        $data_ref->{'CDS'}->{$cds_id_2}->{'gene'} = $gene_2;
        $data_ref->{'CDS'}->{$cds_id_3}->{'gene'} = $gene_3;
    }               
    else {
        die "From CDS-to-gene-to-species tab-delimited table $cds2gene2spp, malformatted input line: $input\n";
    }
}
close $CDS2GENE2SPP;

open my $INPUT_OFND, '<', $input_ofnd;
while (my $input = <$INPUT_OFND>) { 
    # OrthoFinder puts '\r\n' into its files for some goofy reason, so strip the '\r' characters.
    $input =~ s/\r//g;
    my $output = $input;

    # Convert the one CDS gene at the end of the line -- which may, or may not, have \t\n instead of just \n.
    if ( $output =~ /(\S+)(\t*\n)\z/ ) {
        my $cds_id_3 = $1;
        my $end      = $2;
        if ( exists $data_ref->{'CDS'}->{$cds_id_3}->{'gene'} ) {
            my $gene_3 = $data_ref->{'CDS'}->{$cds_id_3}->{'gene'};
            $output =~ s/$cds_id_3$end/$gene_3$end/;
        }
    }
    else {
        die "Cannot parse end of input line: $input";
    }

    # Go through the rest of the line and convert every 'CDS,' to a 'gene,'
    while ( $output =~ / (\S+[,]) /xmsg ) {
        my $cds_id_1 = $1;
        my $gene_1   = q{};
        if ( exists $data_ref->{'CDS'}->{$cds_id_1}->{'gene'} ) {
           $gene_1 = $data_ref->{'CDS'}->{$cds_id_1}->{'gene'};
           $output =~ s/$cds_id_1/$gene_1/g;
       }
    }

    while ( $output =~ / (\S+\t) /xmsg ) {
        my $cds_id_2 = $1;
        my $gene_2   = q{};
        if ( exists $data_ref->{'CDS'}->{$cds_id_2}->{'gene'} ) {
           $gene_2 = $data_ref->{'CDS'}->{$cds_id_2}->{'gene'};
           $output =~ s/$cds_id_2/$gene_2/g;
       }
    }

    print "$output";
}

