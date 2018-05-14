#!/usr/bin/env perl

# genpept2uniprot.pl -- Erich Schwarz <ems@emstech.org>, 7/6/2016.
# Purpose: generate FastA record that somewhat mimics UniProt header format; uses information from NCBI GenPept that NCBI .fa generally lacks.

use strict;
use warnings;
use Getopt::Long;

my $infile         = q{};

my $aa_read_toggle = 0;  # prevent text from being scanned as coding sequences

my $acc_number     = q{};
my $species        = q{};
my $gene_name      = q{};
my $product_name   = q{};

my $orf_seq_line   = q{};

my $help;

GetOptions ( 'infile=s' => \$infile,
             'help'     => \$help,   );

if ( $help or (! -r $infile) ) {
    die "Format: genpept2uniprot.pl --infile|-i [infile, GenPept] > [STDOUT] \n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A LOCUS /xms ) {
        # Set everything to zero.
        $aa_read_toggle = 0;
        $acc_number     = q{};   
        $species        = q{};
        $gene_name      = q{};
        $product_name   = q{};
    }

    # This is the invariant identifier used by NCBI.
    elsif ( $input =~ /\A ACCESSION \s+ (\S+) \s* \z/xms ) {
        $acc_number = $1;

        if ( $aa_read_toggle == 1 ) {
            die "ACCESSION out of expected order: $input\n";
        }
        $aa_read_toggle = 0;
    }

    # new format for GI numbers: VERSION     AAF57416.1  GI:7302326
    # GI numbers will soon be gone; don't use them
    # previous line authenticates data, but use this line to get *full* accession number
    elsif ( $input =~ /\A VERSION \s+ (\S+) /xms) {
        $acc_number = $1;

        if ( $aa_read_toggle == 1 ) {
            die "VERSION out of expected order: $input\n";
        }
        $aa_read_toggle = 0;
    }
    elsif ( $input =~ /\A \s+ ORGANISM \s+ (\S.+\S) \s* \z/xms ) {
        $species = $1;

        $species =~ s/\<\/a>//;         # this is a very specific kludge, designed in
        $species =~ s/\<a href.+\>//;   # 5/15/02 to remove NCBI cruft from headers...

        $species =~ s/[ ]{2,}/ /g;      # trim any extra spaces in the species name
    }

    # There seem to be two inconsistent ways that GenPept provides gene names for proteins.  Enforce consistency.
    elsif ( $input =~ /\s + \/ (?:gene|locus_tag) = (\S+) \s* \z/xms ) {
        my $gene_id = $1;
        if ($gene_name) {
            die "Redundant gene names: $gene_name versus $gene_id (in $input)\n";
        }
        $gene_name = $gene_id;
        $gene_name =~ s/["]//g;  # don't need "gene_name"
    }

    elsif ( $input =~ /\A \s+ \/ product = (.*) \z/xms ) {
        my $product_id = $1;
        $product_id =~ s/\A\s+//;
        $product_id =~ s/\s+\z//;
        $product_id =~ s/["]//g;  # don't need "protein name", either, even if it's long 

        # No attempt to enforce non-redundancy; instead, we record the first value and ignore any later values in the same GenPept record.
        # Meant to keep "serum amyloid A-3 protein precursor" from being later overrun by "Serum amyloid A-3 protein" for the same gene product.
        if (! $product_name) {
            $product_name = $product_id;
        }
    }

    elsif ($input =~ /ORIGIN/) {
        my $fasta_header = q{};

        # Enforce population of all values.
        if (! $acc_number) {
            die "Failed to identify accession ID\n";
        }
        if (! $species) {
            die "Failed to identify species name\n";
        }
        if (! $gene_name) {
            die "Failed to identify gene name\n";
        }
        if (! $product_name) {
            die "Failed to identify product name\n";
        }

        $fasta_header = ">tr|$acc_number|$gene_name $product_name OS=$species GN=$gene_name";

        print "$fasta_header\n";
        $aa_read_toggle = 1;
    }
    elsif ( $input =~ /[a-zA-Z]+/xms ) {
        if ($aa_read_toggle == 1) {
            $orf_seq_line = $input;
            $orf_seq_line =~ s/[\d]*//g;
            $orf_seq_line =~ s/[\s]*//g;
            $orf_seq_line =~ s/[\W]*//g;
            $orf_seq_line =~ tr/a-z/A-Z/;
            print "$orf_seq_line\n";
        }
    }
    elsif ( $input =~ /\A \/ \/ \s* \z/xms ) {
        if ( $aa_read_toggle != 1 ) {
            die "Unexpected end of GenPept product record, not being read for amino acids\n";
        }

        # Set everything to zero.
        $aa_read_toggle = 0;
        $acc_number     = q{};
        $species        = q{};
        $gene_name      = q{};
        $product_name   = q{};
    }
}

