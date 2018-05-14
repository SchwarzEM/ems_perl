#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t ([^\t]*) \z/xms ) { 
        my $uniprot   = $1;
        my $orig_gene = $2;
        my $desc      = $3;
        $desc         =~ s/["]//g;
        $desc         =~ s/\A\s+//;
        $desc         =~ s/\s+\z//;

        my $gene_stem = $orig_gene;
        my $prefix    = q{};
        if ( $gene_stem =~ /\A (EN[A-Z]*\d+) \| (\S+) \z/xms ) { 
            $prefix    = $1;
            $gene_stem = $2;
        }

        # Both acknowledges that we've observed this, and provides a reliable uniprot to gene stem link.
        $data_ref->{'uniprot'}->{$uniprot}->{'gene_stem'}->{$gene_stem}->{'seen'} = 1;

        if ($prefix) {
            $data_ref->{'gene_stem'}->{$gene_stem}->{'prefix'}->{$prefix} = 1;
        }
        if ($desc) {
            # There may well be more than one entry from more than one version of a gene.
            # Note that we need to have primary access through the gene stem, not the uniprot ID first!
            $data_ref->{'gene_stem'}->{$gene_stem}->{'desc'}->{$desc} = 1;
        }
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };

foreach my $uniprot (@uniprots) {
    my @gene_stems = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'gene_stem'} };
    foreach my $gene_stem (@gene_stems) { 
        # Default -- to change if we have prefixes.
        my $final_name = $gene_stem;

        # If there are ENS..\d+ annotations, pick the one with the highest number.
        my @prefixes        = ();
        my $selected_prefix = q{};

        if ( exists $data_ref->{'gene_stem'}->{$gene_stem}->{'prefix'} ) { 
            # Reverse sort to get the highest index number first.
            @prefixes        = reverse sort keys %{ $data_ref->{'gene_stem'}->{$gene_stem}->{'prefix'} };
            $selected_prefix = $prefixes[0];
            $final_name      = $selected_prefix . q{|} . $final_name;
        }

        # Also pool descriptions.
        my @descs     = ();
        my $desc_text = q{};
        if ( exists  $data_ref->{'gene_stem'}->{$gene_stem}->{'desc'} ) { 
            @descs     = sort keys %{ $data_ref->{'gene_stem'}->{$gene_stem}->{'desc'} };
            $desc_text = join '; ', @descs;
            $desc_text = "\"$desc_text\"";
        }           

        # Finally!
        print "$uniprot\t$final_name\t$desc_text\n";
    }
}

