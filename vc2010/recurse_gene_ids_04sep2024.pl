#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

my $names = q{};
my $maps  = q{};

my @set_list = ();

my %seen     = ();
my @queries  = ();

$names = $ARGV[0] if $ARGV[0];
$maps  = $ARGV[1] if $ARGV[1];

if ( (! $names ) or (! $maps ) ) {
     die "Format: recurse_gene_ids_04sep2024.pl [gene names 1+2] [maps of 1-to-2 and 2-to-1] > [equiv groups]\n";
}

open my $MAPS, '<', $maps;
while ( my $input = <$MAPS> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) {
        my $gene = $1;
        my $txt  = $2;
        my @list = split '; ', $txt;
        @{ $data_ref->{'gene'}->{$gene}->{'list'} } = @list;
    }
    else {
        die "From mapping table $maps, cannot parse: $input\n";
    }
}
close $MAPS;

open my $NAMES, '<', $names;
while (my $gene = <$NAMES> ) {
    chomp $gene;
    if ( $gene =~ /\A \S+ \z/xms ) {
        %seen    = ();
        @queries = ();

        push @queries, $gene;

        while (@queries) {
            evaluate(@queries);
        }

        my @set_members = sort keys %seen;
        my $set_txt     = join '...', @set_members;
        push @set_list, $set_txt;
    }
    else {
        die "From gene list $names, cannot parse: $gene\n";
    }
}
close $NAMES;

@set_list = sort(@set_list);
@set_list = uniq(@set_list);

foreach my $set (@set_list) {
    # Add this transformation at the last possible minute because '...' turns out to be unparseable, though useful internally.
    $set =~ s/\.\.\./; /g;
    print "$set\n";
}

sub evaluate {
    # Import contents of global @queries, then zero it out.
    my @_queries = @_;
    @queries = ();

    # Repopulate @queries if and only if we observe *new* matches to known genes.
    foreach my $_gene1 (@_queries) {
        $seen{$_gene1} = 1;
        if ( exists $data_ref->{'gene'}->{$_gene1}->{'list'} ) {
            my @_matches = @{ $data_ref->{'gene'}->{$_gene1}->{'list'} };
            foreach my $_gene2 (@_matches) {
                if (! $seen{$_gene2} ) {
                    push @queries, $_gene2;
                }
            }
        }
    }
    return;
}

