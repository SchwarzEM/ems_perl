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
# ID   O18367_DROME            Unreviewed;       950 AA.
# AC   O18367;
# [...]
# DE   SubName: Full=CALX;
# [...]
# GN   Name=Calx; ORFNames=CG5685;
# DR   FlyBase; FBgn0013995; Calx.

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
    # DR   FlyBase; FBgn0013995; Calx.
    # Other sample inputs:
    # DR   FlyBase; FBgn0019644; ATPsyn-b.
    # DR   FlyBase; FBgn0000409; Cyt-c-p.
    # Ergo, we have to be cautious about trimming the 'human_name' part.
    elsif ( @uniprots and ( $input =~ /\A DR \s+ FlyBase /xms ) ) { 
        if ( $input =~ /\A DR \s+ FlyBase; \s+ (FBgn\d+) (.*) \z/xms ) {
            $mod_name  = $1;
            my $human_name = $2;

            # Default:
            $gene_name = $mod_name;

            # Trim non-name characters, carefully.
            $human_name =~ s/\A;\s+//;
            $human_name =~ s/\.\z//;
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

