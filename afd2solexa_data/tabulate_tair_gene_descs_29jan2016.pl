#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my @genes = ();

# Unless I decide to expand this, there is no point in printing 'type';
# however the capacity might be added in future scripts.
my %useful_types = (
    protein_coding => 1,
);

LOOP:
while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+\.\d+) \t (\S+) \t (.*) /xms ) { 
        my $tx   = $1;
        my $type = $2;
        my $desc = $3;

        if (! exists $useful_types{$type} ) {
            next LOOP;
        }

        if ( $desc !~ /\A [^\t]* \t [^\t]* \t [^\t]* \z/xms ) { 
            die "Cannot parse transcript description in: $input\n";
        }

        my $gene       = $tx;
        my $short_desc = q{};
        my $summary    = q{};
        my $comp_desc  = q{};
        
        $gene =~ s/\.\d+\z//;

        if ( $desc =~ /\A ([^\t]*) \t ([^\t]*) \t ([^\t]*) \z/xms ) {
            $short_desc = $1;
            $summary    = $2;
            $comp_desc  = $3;
            if ( $comp_desc =~ /\A ([^;]+) [;] /xms ) {
                $comp_desc = $1;
            }
            $comp_desc =~ s/\A\s+//;
            $comp_desc =~ s/\s+\z//;
        }

        # Do not bother with this, as long as there is just one type!
        # $data_ref->{'gene'}->{$gene}->{'type'}->{$type} = 1;

        if ( $short_desc =~ /\S/xms ) { 
            $data_ref->{'gene'}->{$gene}->{'short_desc'}->{$short_desc} = 1;
        }

        if ( $summary =~ /\S/xms ) {
            $data_ref->{'gene'}->{$gene}->{'summary'}->{$summary} = 1;
        }   

        if ( $comp_desc =~ /\S/xms ) {
            $data_ref->{'gene'}->{$gene}->{'comp_desc'}->{$comp_desc} = 1;
        }
    }
}

if ( exists $data_ref->{'gene'} ) {
    @genes = sort keys %{ $data_ref->{'gene'} };
}

if (@genes) {
    # header for 'short description, curator summary and computational description'
    print "Gene\tShort_desc\tSummary\tComp_desc\n";

    foreach my $gene (@genes) {
        my @short_descs      = ();
        my $short_desc_text = q{};
        if ( exists $data_ref->{'gene'}->{$gene}->{'short_desc'} ) {
            @short_descs     = sort keys %{ $data_ref->{'gene'}->{$gene}->{'short_desc'} };
            $short_desc_text = join '; ', @short_descs;
            $short_desc_text =~ s/[;]\s+[;]/;/g;
        }

        my @summaries    = ();
        my $summary_text = q{};
        if ( exists $data_ref->{'gene'}->{$gene}->{'summary'} ) { 
            @summaries    = sort keys %{ $data_ref->{'gene'}->{$gene}->{'summary'} };
            $summary_text = join '; ', @summaries;
            $summary_text =~ s/[;]\s+[;]/;/g;
        }

        my @comp_descs     = ();
        my $comp_desc_text = q{};
        if ( exists $data_ref->{'gene'}->{$gene}->{'comp_desc'} ) {
            @comp_descs     = sort keys %{ $data_ref->{'gene'}->{$gene}->{'comp_desc'} };
            $comp_desc_text = join '; ', @comp_descs;
            $comp_desc_text =~ s/[;]\s+[;]/;/g;

        }
        print "$gene\t$short_desc_text\t$summary_text\t$comp_desc_text\n";
    }
}

