#!/usr/bin/env perl

# repartition_cdt.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/1/2010.
# Purpose: given an existing *.cdt readout of genes and two *.kgg clusterings of the same genes, produce a repartitioned *.cdt and *.kgg.

use strict;
use warnings;
use Getopt::Long;
use Tie::RefHash;
use List::MoreUtils qw(uniq);

my $input_cdt     = q{};
my @kggs          = ();
my $output_prefix = q{};
my $help;

my $kgg_num       = 0;
my @kgg_nums      = ();

my %kgg_genes     = ();
my @key_kgg_genes = ();

my @key_cdt_genes      = ();
my %seen_orig_cdt_gene = ();

# We need the magic of Tie::RefHash to let us use arrayrefs as hash keys, 
#     and we need to use it before using the hashes for anything at all.
my %grouplists = ();
tie %grouplists,  "Tie::RefHash";

my $data_ref;

GetOptions ( 'input_cdt=s'      => \$input_cdt,
             'kggs=s{,}'        => \@kggs,
             'output_prefix=s', => \$output_prefix,
             'help'             => \$help, );

if ( ($help) or (! $input_cdt) or (! @kggs) or (! $output_prefix ) ) {
    die "Format: repartition_cdt.pl\n",
        "    -i|--input_cdt [CDT to repartition]\n",
        "    -k|--kggs      [KGGs to use as guides, in order of precedence -- first KGG should match input CDT]\n",
        "    -o|--output_prefix  [prefix for repartitioned CDT and KGG output files; can specify path]\n",
        "    -h|--help\n", 
        ;
}

# Later required to enforce correct number of elements in each grouplist.
my $kgg_count = @kggs;

# Ensure that outputs are possible before bothering with anything else.
my $output_cdt = $output_prefix . '.cdt';
$output_cdt    = safename($output_cdt);
open my $OUTPUT_CDT, '>', $output_cdt or die "Can't open output CDT file $output_cdt: $!";

my $output_kgg = $output_prefix . '.kgg';
$output_kgg    = safename($output_kgg);
open my $OUTPUT_KGG, '>', $output_kgg or die "Can't open output KGG file $output_kgg: $!";

# Get a complete accounting of each gene's status in each KGG file giving clusters.
foreach my $input_kgg (@kggs) { 
    # $kgg_num goes up to 1 for first KGG file; up to N for the Nth file.
    $kgg_num++;
    my $accept_header = 1;

    # Populate @kgg_nums with (1, ... , k):
    push @kgg_nums, $kgg_num;

    my $gene_id = q{};
    my $grp_num = q{};
    open my $INPUT_KGG, '<', $input_kgg  or die "Can't open input KGG file $input_kgg: $!";
    while (my $input = <$INPUT_KGG>) { 
        chomp $input;
        if ( $input =~ /\A (\S+) \s+ (\d+) \s* \z /xms ) { 
            $gene_id = $1;
            $grp_num = $2;

            # Record list of genes from *first* (index) KGG file.  
            # This will be later checked for perfect match to input CDT genes:
            if ( $kgg_num == 1 ) { 
                push @key_kgg_genes, $gene_id;
            }

            # Nonredundantly record each gene found in *any* KGG file:
            $kgg_genes{$gene_id} = 1;

            # Die loudly if any given gene turns out to have two group annotations from any one KGG file:
            if ( exists $data_ref->{'gene'}->{$gene_id}->{'kgg'}->{$kgg_num}->{'group'} ) { 
                die "Gene $gene_id has more than one group annotation in KGG $kgg_num!\n";
            }

            # For each gene, track a unique group affiliation within each KGG.
            #     Though each gene will belong to all KGGs, 
            #     it will belong to only one group within a single KGG:
            $data_ref->{'gene'}->{$gene_id}->{'kgg'}->{$kgg_num}->{'group'} = $grp_num;

            # For each gene, also list, one by one, the group numbers from each KGG in turn.
            # This list does *not* contain KGG file numbers, just group numbers within each KGG file.
            # This will later let us do a series of sortings, from least important group to most, 
            #     which will give us the repartioning of genes we want.
            push @{ $data_ref->{'gene'}->{$gene_id}->{'grouplist'} }, $grp_num;
        }
        else { 
            # Accept only one header line, and only once:
            if ( ( $input =~ /\A Gene \s+ GROUP \s* \z/xms ) and $accept_header ) { 
                $accept_header = 0;
            }
            else { 
                die "Can't parse input line from input KGG file $input_kgg: $input\n";
            }
       }
    }
    close $INPUT_KGG or die "Can't close filehandle to input KGG file $input_kgg: $!";
}

# Before proceeding, make sure that there are no genes without some KGG/group pair and without a grouplist.
foreach my $gene_id (sort keys %kgg_genes) { 
    if (! exists $data_ref->{'gene'}->{$gene_id}->{'grouplist'} ) { 
        die "Gene $gene_id does not have a grouplist!\n";
    }
    foreach my $kgg_num1 (@kgg_nums) { 
        if (! exists $data_ref->{'gene'}->{$gene_id}->{'kgg'}->{$kgg_num1}->{'group'} ) {
            die "Gene $gene_id does not have data for any group in KGG $kgg_num1!\n";
        }
    }
}

# Make a single hash listing arrayrefs for *all* observed grouplists for all genes.
#    This will later be used to produce a properly sorted nonredundant list of group memberships,
#    which in turn will allow sorting of the key CDT's genes into our final partitioning.

foreach my $gene_id ( sort keys %{ $data_ref->{'gene'} } ) { 
    # Enforce correct number of elements in each grouplist.
    my $grouplist_len = @{ $data_ref->{'gene'}->{$gene_id}->{'grouplist'} };
    if ( $grouplist_len != $kgg_count ) { 
        die "The grouplist of gene $gene_id has $grouplist_len members, but should have $kgg_count!\n";
    }

    my @grouplist_mems = @{ $data_ref->{'gene'}->{$gene_id}->{'grouplist'} };
    my $grouplist_ref = \@grouplist_mems;

    # Make a single nonredundant list of all group lists observed for *all* genes.
    # This is possible because of Tie::RefHash; normal hashes won't take references (e.g., arrayrefs) as keys.

    $grouplists{$grouplist_ref} = 1;

    # For each individual gene, enforce uniqueness of its grouplist, then store.
    # Note that we use a plain text string, because trying to use arrayref keys gave highly 
    #    sporadic selections in the final gene list sort!!

    if ( exists $data_ref->{'gene'}->{$gene_id}->{'grouplist_txt'} ) { 
        die "Gene $gene_id has multiple grouplists.\n";
    }
    my $grouplist_txt = join q{.}, @{ $data_ref->{'gene'}->{$gene_id}->{'grouplist'} };
    $data_ref->{'gene'}->{$gene_id}->{'grouplist_txt'} = $grouplist_txt;
}

# Getting the entire input CDT into an array makes life simpler:
#    lines 0-1 (Perl numbering) are the headers and can be left untouched;
#    lines 2+  (Perl numbering) are the data, and can be sorted rationally into a new CDT.

open my $INPUT_CDT, '<', $input_cdt or die "Can't open input CDT file $input_cdt: $!";
my @input_cdt_data = <$INPUT_CDT>;

# Remember to get rid of all those pesky "\n"s.
chomp @input_cdt_data;

# Don't need this filehandle very long...
close $INPUT_CDT or die "Can't close filehandle to input CDT file $input_cdt: $!";

# Obtain, and sanity-check, Perl-style numbers for header + data lines in CDT file:
my $input_cdt_lines = @input_cdt_data;
$input_cdt_lines--;
if ( $input_cdt_lines < 2) { 
    die "No possible data lines in input CDT file $input_cdt\n";
}

# Read array of data lines into a gene-name-indexed hashref:
foreach my $i (2..$input_cdt_lines) { 
    my $data_line = $input_cdt_data[$i];
    my $gene_id = q{};
    if ( $data_line =~ /\A (\S+) \s+ \S .* \z/xms ) { 
        $gene_id = $1;

        # Enforce prior sight of gene ID:
        if (! $kgg_genes{$gene_id} ) { 
            die "Gene $gene_id in input CDT file $input_cdt does not exist in KGG data!\n";
        }

        # Record list of input CDT genes, so that 
        #     they can be later checked for an exact match to the first KGG file's gene list:
        push @key_cdt_genes, $gene_id;

        # Record each line of data for each gene:
        $data_ref->{'CDT_gene'}->{$gene_id}->{'CDT_data'} = $data_line;
    }
    else { 
        die "Can't parse data line of input CDT file $input_cdt: $data_line\n";
    }
}

# Delete the now-unneeded array of data lines: 
delete @input_cdt_data[2..$input_cdt_lines];

# Check lists of genes from input CDT and first KGG file for exact identity.
# Check for identical numbers:

my $key_kgg_genecount = @key_kgg_genes;
my $key_cdt_genecount = @key_cdt_genes;
if ( $key_kgg_genecount != $key_cdt_genecount ) {
    die "Input CDT and key KGG file have different numbers of genes!\n";
}

# Then check for identical names between the two arrays:

$key_kgg_genecount--;
foreach my $i (0..$key_kgg_genecount) { 
    if ( $key_kgg_genes[$i] ne $key_cdt_genes[$i] ) { 
        die "Gene $key_kgg_genes[$i] and CDG gene $key_cdt_genes[$i] are not the same!\n";
    }
}

# Make properly sorted list of grouplist arrayrefs:
my @grouplist_refs = keys %grouplists;
my $j = ($kgg_count - 1);
my @indices = reverse (0..$j);
foreach my $i (@indices) {   
    @grouplist_refs = sort { ${ $a }[$i] <=> ${ $b }[$i] } @grouplist_refs;
}

# Then, map this list from arrayrefs to grouplist_txt strings:
my @grouplist_txts = map { join q{.}, @{ $_ } } @grouplist_refs;
# And enforce uniqueness of each member of list:
@grouplist_txts = uniq(@grouplist_txts);

# Then use that list for a series of greps that, bit by bit, build up the desired gene list.

my @new_cdt_genelist = ();

foreach my $grouplist_txt ( @grouplist_txts ) { 
    foreach my $key_cdt_gene (@key_cdt_genes) { 
        if ( $data_ref->{'gene'}->{$key_cdt_gene}->{'grouplist_txt'} eq $grouplist_txt ) { 
            push @new_cdt_genelist, $key_cdt_gene;
            $seen_orig_cdt_gene{$key_cdt_gene}++;
            if ( $seen_orig_cdt_gene{$key_cdt_gene} >= 2 ) { 
                die "Trying to add gene $key_cdt_gene twice to new list!\n";
            }
        }
        elsif ( $data_ref->{'gene'}->{$key_cdt_gene}->{'grouplist_txt'} ne $grouplist_txt ) {
        }
        else { 
            die "Can't parse gene $key_cdt_gene for match or non-match to grouplist txt \"$grouplist_txt\"\n";
        }
    }
}

# At last!  I have the sorted list I wanted.  Print the new CDT:
print $OUTPUT_CDT "$input_cdt_data[0]\n";
print $OUTPUT_CDT "$input_cdt_data[1]\n";

foreach my $new_gene (@new_cdt_genelist) { 
    print $OUTPUT_CDT "$data_ref->{'CDT_gene'}->{$new_gene}->{'CDT_data'}\n";
}

# Done with CDT:
close $OUTPUT_CDT or die "Can't close filehandle to output CDT file $output_cdt: $!";

# Now, print the new KGG, starting with its header line:
print $OUTPUT_KGG "Gene\tGROUP\n";
# And moving on to its new multi-valued grouping:
foreach my $new_gene (@new_cdt_genelist) { 
    print $OUTPUT_KGG "$new_gene\t$data_ref->{'gene'}->{$new_gene}->{'grouplist_txt'}\n";
}
# Done!
close $OUTPUT_KGG or die "Can't close filehandle to output KGG file $output_kgg: $!";

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

