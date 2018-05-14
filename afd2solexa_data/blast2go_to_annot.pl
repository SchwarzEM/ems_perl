#!/usr/bin/env perl

use strict;
use warnings;

my $go_annot = $ARGV[0];
my $blast2go = $ARGV[1];

my $gene          = q{};
my $go_term       = q{};
my $go_desc       = q{};

my %go_terms_seen    = ();
my @go_term_synonyms = ();

my $data_ref;

if ( (! $go_annot) or (! $blast2go ) ) { 
    die "Format: blast2go_to_annot.pl [GO term/descs OBO] [blast2go results] > [output]\n";
}

open my $GO, '<', $go_annot or die "Can't open GO annotations $go_annot\n";
while (my $input = <$GO>) { 
    chomp $input;
    if ( $input =~ /\A id[:] \s+ (GO[:]\d+) /xms ) { 
        # Capture latest.
        $go_term                 = $1;

        # Clear out backlogged data, which due to obo format I can't process until I *know* I'm out of alt IDs.
        @go_term_synonyms = sort keys %go_terms_seen;
        foreach my $go_term_synonym (@go_term_synonyms) {
            $data_ref->{'go_term'}->{$go_term_synonym} = $go_desc;
        }

        # Then zero out everything:
        %go_terms_seen           = ();
        @go_term_synonyms        = ();

        # And start the data cycle anew:
        $go_terms_seen{$go_term} = 1;
    }
    elsif ( $input =~ /\A name[:] \s (.+\S) \s* \z/xms ) {
        # Capture this datum but don't act on it until I've got all synonyms harvested:
        $go_desc = $1;
    }
    elsif ( $input =~ /\A alt_id[:] \s+ (GO[:]\d+) /xms ) {
        $go_term                 = $1;
        $go_terms_seen{$go_term} = 1;
    }
}

# One last round of data export:

@go_term_synonyms = sort keys %go_terms_seen;
foreach my $go_term_synonym (@go_term_synonyms) {
    $data_ref->{'go_term'}->{$go_term_synonym} = $go_desc;
}

close $GO or die "Can't close filehandle to GO annotations $go_annot\n";

open my $BLAST, '<', $blast2go or die "Can't open blast2go file $blast2go\n";
while (my $input = <$BLAST>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (GO[:]\d+) \s* \z/xms ) { 
        $gene = $1;
        $go_term  = $2;

        # Allow >= 2 terms/gene:
        $data_ref->{'gene'}->{$gene}->{'go_term'}->{$go_term} = 1;
    }
}
close $BLAST or die "Can't close filehandle to blast2go file $blast2go\n";

my $header = "Gene\tGO_terms\n";

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @go_terms = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'go_term'} };
    foreach my $go_term1 (@go_terms) { 
        if (! exists $data_ref->{'go_term'}->{$go_term1} ) { 
            die "Can't describe $go_term1\n";
        }
    }
    my @go_descs      = map { "$data_ref->{'go_term'}->{$_} [$_]" } @go_terms;
    my $go_descs_text = q{};
    if (@go_descs) { 
        $go_descs_text = join '; ', @go_descs;
    }
    print $header if $header;
    $header = q{};
    print "$gene1\t$go_descs_text\n";
}

