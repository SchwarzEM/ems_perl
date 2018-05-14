#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my @uniprots    = ();
my @gene_names  = ();

my $description = q{};

my $gene_name   = q{};
my $mod_name    = q{};
my $cds_name    = q{};

# Sample input -- but note that there can be two or more different GN, etc. per protein.
# Also note that UniProt does not always agree with WormBase about how to name a gene (e.g., LGO).
#
# GN   ORFNames=SPAC2E1P5.01c, SPAPB1E7.13c;
# GN   Name=moa1; Synonyms=mug159; ORFNames=SPAC15E1.07c;
# GN   Name=mob1; ORFNames=SPBC428.13c;

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\s+\z//;

    # '//' == end of entry.
    if ( $input =~ /\A \/\/ /xms ) { 
        if (@uniprots) {
            map_print_and_clear_stored_data();
        }
    }

    # Note that there can be more than one AC line per entry; important not to overwrite first line's accessions with the 2+ lines.
    if ( $input =~ /\A AC \s+ (\S.*\S) \z/xms ) { 
        my $uniprot_text = $1;
        $uniprot_text =~ s/;//g;
        my @new_uniprots = split /\s+/, $uniprot_text;
        push @uniprots, @new_uniprots;
    }
    elsif ( @uniprots and ( $input =~ /\A DE \s+ RecName: \s+ Full= (.+) ; \z/xms ) ) {
        $description = $1;
    }

    # To try to capture everything, *hope* that UniProt consistently puts the correct CDS right before the ';'.
    # Note that having .* be nongreedy ? was crucial here!
    elsif ( @uniprots and ( $input =~ /\A GN \s+ .* ORFNames= .*? (\S+) ; /xms ) ) {
        $gene_name = $1;

        # Sometimes, but not always, this will happen:
        if ( $input =~ /\A GN \s+ Name= (\S+) ; /xms ) {
            my $human_name = $1;
            $gene_name = $gene_name . q{|} . $human_name;
        }
        push @gene_names, $gene_name;
    }
}

# Fail-safe, but shouldn't be needed.
map_print_and_clear_stored_data();

# This has been recoded to explicitly handle any possible mixture of 0 to 2+ gene names or ENSEMBL names per protein.
sub map_print_and_clear_stored_data {
    foreach my $uniprot (@uniprots) { 
        # We will only print data for UniProt proteins which have either an ENSEMBL or a human-readable gene name.
        if (@gene_names) {
            my @print_lines = ();

            # The description, at any rate, should be invariant per UniProt entry.
            $description =~ s/\A["]+//;
            $description =~ s/["]+\z//;
            if ( $description =~ /\S/xms ) {
                $description = "\"$description\"";
            }

            # If we do have gene names:
            if (@gene_names) { 
                # In this particular case, we have already made every gene name a long name!
                # Which simplifies life a great deal.
                foreach my $gene_name1 (@gene_names) {
                    my $print_text = "$uniprot\t$gene_name1\t$description";
                    push @print_lines, $print_text;
                }
            }
            foreach my $print_text (@print_lines) {
                print "$print_text\n";
            }
        }
    }

    # With lines printed, return everything to zero for the next entry.
    @uniprots    = ();   
    @gene_names  = ();
    
    $description = q{};   
    $gene_name   = q{};
    $mod_name    = q{};

    return;
}

