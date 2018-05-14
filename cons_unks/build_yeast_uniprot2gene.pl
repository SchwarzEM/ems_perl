#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my $id_table_1 = $ARGV[0];
my $id_table_2 = $ARGV[1];

open my $ID1, '<', $id_table_1 or die "Can't open yeast ID table 1 $id_table_1: $!";
while (my $input = <$ID1>) { 
    chomp $input;
    # Sample input:
    # P38903  SGD     S000005540  
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) { 
        my $uniprot = $1;
        my $db_ref  = $2;
        my $sgdid   = $3;
        if ( $db_ref eq 'SGD' ) {
            if ( $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} ) {
                my @known_sgdids = sort keys %{$data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} };
                push @known_sgdids, $sgdid;
                my $known_sgdid_text = join ' versus ', @known_sgdids;
                warn "Inconsistent mapping of UniProt $uniprot to gene IDs $known_sgdid_text\n";
            }
            if ( $sgdid !~ /\A S\d+ \z/xms ) { 
                die "Can't make sense of putative SGDID $sgdid in input from yeast ID table 1 $id_table_1: $input\n";
            }
            # Note: we have to do this, because there are instances of *identical* proteins encoded by two or more genes...
            $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'}->{$sgdid} = 1;
        }
    }
    else { 
        die "Can't parse input from yeast ID table 1 $id_table_1: $input\n";
    }
}
close $ID1 or die "Can't close filehandle to yeast ID table 1 $id_table_1: $!";

open my $ID2, '<', $id_table_2 or die "Can't open yeast ID table 2 $id_table_2: $!";
while (my $input = <$ID2>) {
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        # Sample input:
        # >YAL001C TFC3 SGDID:S000000001,
        if ( $input =~ /\A > (\S+) \s (\S+) \s SGDID:(S\d+) , /xms ) { 
            my $genome_name = $1;
            my $common_name = $2;
            my $sgdid       = $3;
            my $long_name = q{};
            if ( $genome_name ne $common_name ) { 
                $long_name = $sgdid . q{|} . $genome_name . q{|} . $common_name;
            }
            if ( $genome_name eq $common_name ) { 
                $long_name = $sgdid . q{|} . $genome_name;
            }
            if ( ( $data_ref->{'mod_gene_id'}->{$sgdid}->{'long_name'} ) and ( $data_ref->{'mod_gene_id'}->{$sgdid}->{'long_name'} ne $long_name ) ) {
                die "Inconsistent mapping of SGDID $sgdid to long gene names $data_ref->{'mod_gene_id'}->{$sgdid}->{'long_name'} and $long_name\n";
            }
            $data_ref->{'mod_gene_id'}->{$sgdid}->{'long_name'} = $long_name;
        }
        else {
            die "Can't parse input from yeast ID table 2 $id_table_2: $input\n";
        }
    }
}
close $ID2 or die "Can't close filehandle to yeast ID table 2 $id_table_2: $!";

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };
my $uniprot_count = @uniprots;
foreach my $uniprot (@uniprots) {
    my @sgdids = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} };
    my $sgdid_count = @sgdids;
    foreach my $sgdid (@sgdids) {
        if ( $data_ref->{'mod_gene_id'}->{$sgdid}->{'long_name'} ) { 
            my $long_name = $data_ref->{'mod_gene_id'}->{$sgdid}->{'long_name'};
            print "$uniprot\t$long_name\n";
        }
        else {
            warn "Can't find long name for uniprot $uniprot and sgdid $sgdid\n";
        }
    }
}

