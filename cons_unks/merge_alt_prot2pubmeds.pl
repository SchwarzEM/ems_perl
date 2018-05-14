#!/usr/bin/env perl

use strict;
use warnings;

# merge_alt_prot2pubmeds -- Erich Schwarz <ems394@cornell.edu>, 8/21/2013.
# Purpose: given two different, conflicting UniProt to PubMed mapping tables, merge them.  If a reference is 'primary', do not also let it be secondary.

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^\t]*) \t ([^\t]*) \z/xms ) { 
        my $uniprot       = $1;
        my $prim_ref_text = $2;
        my $sec_ref_text  = $3;
        map_refs_to_prot($uniprot, 'primary', $prim_ref_text);
        map_refs_to_prot($uniprot, 'secondary', $sec_ref_text);
    }
    else {
        die "Can't parse input line: $input\n";
    }
}

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };
foreach my $uniprot (@uniprots) { 
    my @prim_refs = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'primary'} };
    foreach my $prim_ref (@prim_refs) { 
        if ( exists $data_ref->{'uniprot'}->{$uniprot}->{'secondary'}->{$prim_ref} ) { 
            delete $data_ref->{'uniprot'}->{$uniprot}->{'secondary'}->{$prim_ref};
        }
    }
    my @sec_refs = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'secondary'} };
    my $prim_ref_text = join '; ', @prim_refs;
    my $sec_ref_text  = join '; ', @sec_refs;
    print "$uniprot\t$prim_ref_text\t$sec_ref_text\n";
}

sub map_refs_to_prot {
    my $_uniprot  = $_[0];
    my $_type     = $_[1];
    my $_ref_text = $_[2];
    my @_refs     = ();

    if ( $_ref_text =~ /\S/xms ) { 
        @_refs = split /; /, $_ref_text;
        foreach my $_ref (@_refs) { 
            $data_ref->{'uniprot'}->{$_uniprot}->{$_type}->{$_ref} = 1;
        }
    }

    # Do this so that I am guaranteed a full list of UniProts, even if they have no references!
    $data_ref->{'uniprot'}->{$_uniprot}->{'seen'} = 1;

    return;
}

