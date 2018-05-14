#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my @uniprots    = ();
my @gene_names  = ();

my $description = q{};
my $gene_name   = q{};
my $mod_name    = q{};

# Sample input -- but note that there can be two or more different GN and MOD ID per protein.
# 
# ID   143BA_DANRE             Reviewed;         244 AA.
# AC   Q5PRD0; A3KNI9;
# [...]
# DE   RecName: Full=14-3-3 protein beta/alpha-A;
# [...]
# GN   Name=ywhaba; Synonyms=ywhab1; ORFNames=wu:fb80c08;
# [...]
# DR   ZFIN; ZDB-GENE-030131-6583; ywhaba.  ** simplest thing to do is just to mine this! **

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\s+\z//;

    # '//' == end of record.
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

    # Only one description per gene.
    elsif ( @uniprots and ( $input =~ /\A DE \s+ RecName: \s+ Full= (.+) ; \z/xms ) ) {
        $description = $1;
    }

    # Sample input:
    # DR   ZFIN; ZDB-GENE-030131-6583; ywhaba.
    elsif ( @uniprots and ( $input =~ /\A DR \s+ ZFIN /xms ) ) { 
        if ( $input =~ /\A DR \s+ ZFIN; \s+ (ZDB[-]GENE[-]\d+[-]\d+) (.*) \z/xms ) {
            $mod_name  = $1;
            my $human_name = $2;

            # Default:
            $gene_name = $mod_name;

            # Trim non-name characters.
            $human_name =~ s/\W//g;  # e.g., ';' and '.'.
            $human_name =~ s/\s//g;

            # If there's anything left, append it to $gene_name.
            if ($human_name) { 
                $gene_name = $gene_name . q{|} . $human_name;
            }
            push @gene_names, $gene_name;
        } 
        else { 
            die "Can't parse input line: $input\n";
        }
    }
}

# Clear out last stored data, after the text is done running.
map_print_and_clear_stored_data();

# This has been recoded to explicitly handle any possible mixture of 0 to 2+ gene names or ENSEMBL names per protein.
sub map_print_and_clear_stored_data {
    foreach my $uniprot (@uniprots) { 
        # We only have @gene_names to keep track of, so life is simple.
        if (@gene_names) {
            # The description, at any rate, should be invariant per UniProt entry.
            $description =~ s/\A["]+//;
            $description =~ s/["]+\z//;
            if ( $description =~ /\S/xms ) {
                $description = "\"$description\"";
            }

            foreach my $gene_name1 (@gene_names) { 
                print "$uniprot\t$gene_name1\t$description\n";
            }
        }
    }

    # With lines printed, return everything to zero for the next entry.
    @uniprots    = ();   
    @gene_names  = ();
    
    $description = q{};   
    $gene_name   = q{};

    return;
}

