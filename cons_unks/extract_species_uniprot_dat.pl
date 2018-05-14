#!/usr/bin/env perl

# extract_species_uniprot_dat.pl -- Erich Schwarz <ems394@cornell.edu>, 8/22/2013.
# Purpose: given 1+ files with 1+ NCBI Tax ID numbers, and 1+ *.dat files with UniProt records, print each species' data to individual, safenamed output files. 
# Note documentation for .dat file format at: http://web.expasy.org/docs/userman.html
#    which gave useful explanation for 'OX' file line here, but is useful generally.
#    "OX	Taxonomy cross-reference	Once"

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my @taxfiles_input = ();
my @datfiles_input = ();
my $outfile_prefix = q{};

my $data_ref;

my $help;

GetOptions ( 'taxfiles=s{,}'    => \@taxfiles_input,
             'datfiles=s{,}'    => \@datfiles_input,
             'outfile_prefix=s' => \$outfile_prefix,
             'help'             => \$help, 
);

if (! $outfile_prefix) { 
    $outfile_prefix = 'uniprot_';
}

if ($help or (! @taxfiles_input) or (! @datfiles_input) ) { 
    die "Format: extract_species_uniprot_dat.pl\n",
        "    --taxfiles|-t        [NCBI taxon ID list file(s)]\n",
        "    --datfiles|-d        [UniProt .dat file(s) from which to extract taxon-specific data files]\n",
        "    --outfile_prefix|-o  [Prefix for outfiles; default is \"uniprot_\"]\n",
        "    --help|-h            [Print this message]\n",
        ;
}

foreach my $taxfile (@taxfiles_input) { 
    open my $TAX, '<', $taxfile or die "Can't open NCBI Tax ID number file $taxfile: $!";
    while (my $input = <$TAX>) { 
        chomp $input;
        # Note that this format allows taxids to have trailing space and comments.
        # E.g.:  559292   # OS   Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast).
        if ( $input =~ /\A (\d+) /xms ) { 
            my $tax_id = $1;
            $data_ref->{'chosen_tax_id'}->{$tax_id} = 1;
        }
        else { 
            warn "Could not parse line from NCBI Tax ID number file $taxfile: $input\n";
        }
    }
    close $TAX or die "Can't close filehandle to NCBI Tax ID number file $taxfile: $!";
}

# Require at least *one* number entered from one file:
if (! $data_ref->{'chosen_tax_id'}) {
    die "Failed to enter any NCBI Tax ID numbers from file(s): @taxfiles_input\n";
}

foreach my $datfile (@datfiles_input) { 
    open my $DAT, '<', $datfile or die "Can't open UniProt .dat file $datfile: $!";

    my @obs_ncbi_taxa = ();
    my @stored_text   = ();

    while (my $input = <$DAT>) {
        chomp $input;
        push @stored_text, $input;

        # Conveniently, each record in a UniProt .dat ends with '//'.
        if ( $input =~ /\A \/\/ /xms ) { 
            @obs_ncbi_taxa = sort @obs_ncbi_taxa;
            @obs_ncbi_taxa = uniq @obs_ncbi_taxa;
            foreach my $obs_taxon (@obs_ncbi_taxa) {
                push @{ $data_ref->{'tax_id'}->{$obs_taxon}->{'data_subset'} }, @stored_text;
            }
            @obs_ncbi_taxa = ();
            @stored_text   = ();
        }
        # Sample trigger line:
        # OX   NCBI_TaxID=654924;
        elsif ( $input =~ /\A OX \s+ NCBI_TaxID = (\d+) ; /xms ) {
            my $tax_id = $1;
            if ( $data_ref->{'chosen_tax_id'}->{$tax_id} ) { 
                push @obs_ncbi_taxa, $tax_id;
            }
        }
    }
    close $DAT or die "Can't close filehandle to UniProt .dat file $datfile: $!";
}

my @ncbi_ids = sort keys %{ $data_ref->{'chosen_tax_id'} };
foreach my $tax_id (@ncbi_ids) { 
    if ( exists $data_ref->{'tax_id'}->{$tax_id}->{'data_subset'} ) {
        my $output_file = $outfile_prefix . $tax_id . ".dat";
        $output_file    = safename($output_file);
        open my $OUT, '>', $output_file or die "Can't open output subset data file $output_file: $!";
        foreach my $text_line ( @{ $data_ref->{'tax_id'}->{$tax_id}->{'data_subset'} } ) { 
            print $OUT "$text_line\n";
        }
        close $OUT or die "Can't close filehandle to output subset data file $output_file: $!";
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

