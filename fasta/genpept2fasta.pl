#!/usr/bin/env perl

# genpept_to_clean_tfa.pl -- Erich Schwarz <ems@emstech.org>; first written, 9/2/2003; modernized, 11/6/2015.
# Purpose: generate FastA record that has information in headers from GenPept (that NCBI .fa generally lacks).

use strict;
use warnings;
use Getopt::Long;

my $infile         = q{};
my $format         = q{};

my $aa_read_toggle = 0;  # prevent text from being scanned as coding sequences

my $acc_number     = q{};
my $species        = q{};
my $gene_name      = q{};
my $product_name   = q{};
my $gi_number      = q{};

my $orf_seq_line   = q{};

my $help;

GetOptions ( 'infile=s' => \$infile,
             'format=s' => \$format,
             'help'     => \$help,   );

if ( $help or (! -r $infile) or ( ( $format ne 'gi' ) and ( $format ne 'acc' ) ) ) {
    die "Format: genpept2fasta.pl --infile|-i [infile, GenPept] --format|-f [format: either \"gi\" or \"acc\"] > [STDOUT] \n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A LOCUS /xms ) {
        $aa_read_toggle = 0;
    }

    # This is the invariant identifier used by NCBI.
    elsif ( $input =~ /\A ACCESSION \s+ (\S+) \s* \z/xms ) {
        $acc_number = $1;

        if ( $aa_read_toggle == 1 ) {
            die "ACCESSION out of expected order: $input\n";
        }
        $aa_read_toggle = 0;

        $species      = "no_species_name";
        $gene_name    = q{};
        $product_name = q{};
    }

    # new format for GI numbers: VERSION     AAF57416.1  GI:7302326
    # previous line authenticates data, but use this line to get *full* accession number
    elsif ( $input =~ /\A VERSION \s+ (.+) \s+ GI [:] ([0-9]+) \s* \z/xms) {
        $acc_number = $1;
        $gi_number  = $2;

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
        if ($product_name) {
            die "Redundant product names: $product_name versus $product_id (in $input)\n";
        }
        $product_name = $product_id;  # *do* keep "long wordy product name" in double quotes
    }

    elsif ($input =~ /ORIGIN/) {
        my $fasta_header = q{};

        if (! $gene_name) {
            $gene_name = 'no_gene_name';
        }
        if (! $product_name) {
            $product_name = 'no_product_name';
        }

        if ($format eq 'acc') {
            $fasta_header = '>' 
                            . "$acc_number" 
                            . "   gi|$gi_number" 
                            . "   gene:$gene_name" 
                            . "   product:$product_name" 
                            . "   [$species]"
                            ;
        }
        # $format eq 'gi':
        else {
            $fasta_header = '>'
                            . "gi|$gi_number" 
                            . "   acc|$acc_number" 
                            . "   gene:$gene_name"
                            . "   product:$product_name" 
                            . "   [$species]"     
                            ;
        }
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
        $aa_read_toggle = 0;
    }
}

