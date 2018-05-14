#!/usr/bin/env perl

# merge_two_uniprot2gene_mappings.pl -- Erich Schwarz <ems394@cornell.edu>, 8/24/2013.
# Purpose: given two UniProt to gene mappings, either create versions with consistent names, or make synonym to final gene table.

# Note: first version which actually has complete function for the 'consistent' argument.

use strict;
use warnings;

use Getopt::Long;

my $data_ref;

my $consistent;
my $nametable;

my $primary_mapping    = q{};
my @secondary_mappings = ();
my $output_suffix      = q{};
my $help;

GetOptions ( 'consistent'     => \$consistent,
             'nametable'      => \$nametable,
             'primary=s'      => \$primary_mapping,
             'secondary=s{,}' => \@secondary_mappings,
             'output=s'       => \$output_suffix,
             'help'           => \$help, 
);

if (! $output_suffix) { 
    $output_suffix = '.rev.txt';
}

if (     $help 
     or (! $primary_mapping) 
     or (! @secondary_mappings) 
     or ( $consistent and $nametable )
     or ( (! $consistent ) and (! $nametable ) )
   ) { 
    die "Format: merge_two_uniprot2gene_mappings.pl\n",
        "    --consistent|-c   [Take two or more UniProt to gene mappings, and print versions which have consistent naming]\n",
        "    --nametable|-n    [Take one or more UniProt to gene mappings, and print a table which maps synonyms to final name versions]\n",
        "    --primary|-p      [The one UniProt to gene mapping file whose names, in a full-name conflict, take precedence over others]\n",
        "    --secondary|-s    [One or more UniProt to gene mapping files whose names are to be made consistent with the primary file's names]\n",
        "    --output|-o       [For 'nametable', name of output file; for 'consistent', suffix to append to revised UniProt to gene mappings]\n",
        "    --help|-h         [Print this message]\n",
        ;
}

open my $PRIMARY, '<', $primary_mapping or die "Can't open primary mapping file $primary_mapping: $!";
while (my $input = <$PRIMARY> ) { 
    chomp $input;

    if ( $input =~ /\A (\S+) \t (\S+) \b .* \z/xms ) {
        my $uniprot   = $1;
        my $orig_name = $2;

        record_gene_versions('primary_first_word', $orig_name);
    }

    # Require that the parsing work on every single line.
    else {
        die "In primary mapping file $primary_mapping, cannot parse text line: $input;\n";
    }

    # Record successfully parsed lines, for future revision and export.
    push @{ $data_ref->{'text_lines'}->{$primary_mapping} }, $input;
}
close $PRIMARY or die "Can't close filehandle to primary mapping file $primary_mapping: $!";

foreach my $secondary_mapping (@secondary_mappings) {
    open my $SECONDARY, '<', $secondary_mapping or die "Can't open secondary mapping file $secondary_mapping: $!";
    while (my $input = <$SECONDARY>) {
        chomp $input;

        if ( $input =~ /\A (\S+) \t (\S+) \b .* \z/xms ) {
            my $uniprot   = $1;
            my $orig_name = $2;

            # Note that I am requiring 2+ secondary files to *agree* on genes.  Simpler than coding 3-way (or more!) decisions.
            record_gene_versions('secondary_first_word', $orig_name);
        }
        else {
            die "In secondary mapping file $secondary_mapping, cannot parse text line: $input;\n";
        }

        push @{ $data_ref->{'text_lines'}->{$secondary_mapping} } , $input;
    }
    close $SECONDARY or die "Can't close filehandle to secondary mapping file $secondary_mapping: $!";
}

if ($consistent) { 
    my $primary_rev_mapping = $primary_mapping . $output_suffix;
    $primary_rev_mapping  = safename($primary_rev_mapping);

    open my $REV_PRIME, '>', $primary_rev_mapping or die "Can't open revised primary mapping file $primary_rev_mapping: $!";
    foreach my $primary_map_line ( @{ $data_ref->{'text_lines'}->{$primary_mapping} } ) {
        $primary_map_line = convert_gene_names($primary_map_line);
        print $REV_PRIME "$primary_map_line\n";
    }
    close $REV_PRIME or die "Can't close filehandle to revised primary mapping file $primary_rev_mapping: $!";

    foreach my $secondary_mapping (@secondary_mappings) {
        my $secondary_rev_mapping = $secondary_mapping . $output_suffix;
        $secondary_rev_mapping    = safename($secondary_rev_mapping);

        open my $REV_SECOND, '>', $secondary_rev_mapping or die "Can't open revised secondary mapping file $secondary_rev_mapping: $!";
        foreach my $secondary_map_line ( @{ $data_ref->{'text_lines'}->{$secondary_mapping} } ) {
            $secondary_map_line = convert_gene_names($secondary_map_line);
            print $REV_SECOND "$secondary_map_line\n";
       }
       close $REV_SECOND or die "Can't close filehandle to revised secondary mapping file $secondary_rev_mapping: $!";
    }
}

# The following subroutine, 'record_gene_versions', does two things:
#
# 1. It will record varying lengths of names in a recoverable way;  
#    we will always have the longest version of a name and its shorter versions.
#
# 2. Within a given mapping file, it will require that no two names of identical 
#    word length may share a common first word, but differ in the rest of their content.
#
# These will save no end of headaches later when we need to merge names and come up with synonym tables.
# We can always ask, for a given first name, what its *longest* version in a primary or secondary table is
# We can then choose to pick the longer of two options, either within a mapping or between them.
# And we can prefer the primary mapping of length X to the secondary one of length X.

sub record_gene_versions {             
    my $_mapping_first_word = $_[0];  # E.g., 'primary_first_word'
    my $_orig_name          = $_[1];

    # There are N-1 staves ('|') splitting N words in a gene name.
    my $_word_count = ( $_orig_name =~ tr/[|]/[|]/ );
    $_word_count++;

    my $_first_word = $_orig_name;
    $_first_word    =~ s/[|]\S+\z//;

    # Within a given mapping, require that each first name map to only one version of X words in length.
    if ( ( exists $data_ref->{$_mapping_first_word}->{$_first_word}->{$_word_count} )
             and ( $_orig_name ne $data_ref->{$_mapping_first_word}->{$_first_word}->{$_word_count} )
       ) {
        die "In primary mapping file $primary_mapping, contradictory gene names:",
            " $data_ref->{$_mapping_first_word}->{$_first_word}->{$_word_count} versus $_orig_name\n",
            ;
    }
        
    # Map the name and its length to its first word.
    $data_ref->{$_mapping_first_word}->{$_first_word}->{$_word_count} = $_orig_name;
        
    # For names with 2-3 words, populate their shorter versions as well.
    while ( $_word_count > 0 ) {
            $_word_count--;
            $data_ref->{$_mapping_first_word}->{$_first_word}->{$_word_count} = $_orig_name;
    }
}

sub convert_gene_names {
    my $_input = $_[0];
    if ( $_input =~ /\A (\S+) \t (\S+) \b (.*) \z/xms ) { 
        my $_uniprot   = $1;
        my $_orig_name = $2;
        my $_comment   = $3;

        # Default: no change!
        my $_revision = $_orig_name;

        my $_first_word = $_orig_name;
        $_first_word    =~ s/[|]\S+\z//;

        # Figure out what the longest version of the name is.
        my $_prim_word_len = 0;
        my $_sec_word_len  = 0;

        if ( exists $data_ref->{'primary_first_word'}->{$_first_word} ) { 
            my @_prim_lens  = sort { $b <=> $a } 
                              keys %{ $data_ref->{'primary_first_word'}->{$_first_word} };
            $_prim_word_len = $_prim_lens[0];
        }

        if ( exists $data_ref->{'secondary_first_word'}->{$_first_word} ) {
            my @_sec_lens   = sort { $b <=> $a } 
                              keys %{ $data_ref->{'secondary_first_word'}->{$_first_word} };
            $_sec_word_len  = $_sec_lens[0];
        }

        # Assign a revision.  Break ties by preferring the primary mapping's version.
        if ( $_prim_word_len >= $_sec_word_len ) { 
            $_revision = $data_ref->{'primary_first_word'}->{$_first_word}->{$_prim_word_len};
        }
        if ( $_sec_word_len > $_prim_word_len ) { 
            $_revision = $data_ref->{'secondary_first_word'}->{$_first_word}->{$_sec_word_len};
        }

        # Finally, incorporate the (possibly) altered gene name into a revised text line.
        return "$_uniprot\t$_revision$_comment";
    }
    else {
        die "Subroutine convert_gene_names cannot parse text line: $_input;\n";
   }
}

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}


