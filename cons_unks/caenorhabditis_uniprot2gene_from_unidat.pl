#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my @uniprots    = ();
my @gene_names  = ();

my %cds2gene    = ();

my $description = q{};

my $gene_name   = q{};
my $mod_name    = q{};
my $cds_name    = q{};

# Sample input -- but note that there can be two or more different GN, etc. per protein.
# Also note that UniProt does not always agree with WormBase about how to name a gene (e.g., LGO).
# 
# ID   14331_CAEEL             Reviewed;         248 AA.
# AC   P41932; Q21537;
# [...]
# DE   RecName: Full=14-3-3-like protein 1;
# [...]
# GN   Name=par-5; Synonyms=ftt-1; ORFNames=M117.2;  ** Note that there are *always* ORFNames= but not always Name= **
# [...]
# DR   WormBase; M117.2; CE06200; WBGene00003920; par-5.
# [...]

# More problems with this:
# GN   ORFNames=R05D3.6;
# GN   and
# GN   ORFNames=ZC262.5;
# [...]
# DR   WormBase; R05D3.6; CE00285; WBGene00019880; -.
# DR   WormBase; ZC262.5; CE00285; WBGene00022582; -.

# Yet another problem, argh.
# GN   Name=let-23; ORFNames=CELE_ZK1067.1, ZK1067.1;
# DR   WormBase; ZK1067.1c; CE42910; WBGene00002299; let-23.

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

    # For C. elegans, we will accept the UniProt version of full gene names for now, but keep in mind that there may be some inconsistencies with WS235.
    # Also, we have to build up a single monolithic name in two steps, with only CDS names as links of the IDs.
    # Note that '.*' flanks allow the OrderedLocusNames to be anywhere in the line.

    # To try to capture everything, *hope* that UniProt consistently puts the correct CDS right before the ';'.
    # Note that having .* be nongreedy ? was crucial here!
    elsif ( @uniprots and ( $input =~ /\A GN \s+ .* ORFNames= .*? (\S+) ; /xms ) ) {
        $cds_name  = $1;
        $gene_name = q{};
        # Trim any trailing [a-z] from CDS names.
        $cds_name =~ s/[a-z]\z//;

        # Sometimes, but not always, this will happen:
        if ( $input =~ /\A GN \s+ Name= (\S+) ; /xms ) {
            $gene_name = $1;
            $cds2gene{$cds_name} = $gene_name;
        }
    }
    elsif ( @uniprots and ( $input =~ /\A DR \s+ WormBase; \s+ (\S+) ; \s+ \S+ ; \s+ (WBGene\d+) ; /xms ) ) {
        $cds_name = $1;
        $mod_name = $2;
        # Trim any trailing [a-z] from CDS names.
        $cds_name =~ s/[a-z]\z//;
        $gene_name = $mod_name . q{|} . $cds_name;
        if ( exists $cds2gene{$cds_name} ) { 
            my $cgc_name = $cds2gene{$cds_name};
            $gene_name = $gene_name . q{|} . $cgc_name;
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
    %cds2gene    = ();
    
    $description = q{};   
    $gene_name   = q{};
    $mod_name    = q{};

    return;
}

