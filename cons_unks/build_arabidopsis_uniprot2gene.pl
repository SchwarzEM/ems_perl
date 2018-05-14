#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my $id_table_1 = $ARGV[0];
my $db_name    = $ARGV[1];
my $id_table_2 = $ARGV[2];

open my $ID1, '<', $id_table_1 or die "Can't open MOD ID table 1 $id_table_1: $!";
while (my $input = <$ID1>) { 
    chomp $input;
    # Sample input:
    # Q9LPP4	TAIR	AT3G10525
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) { 
        my $uniprot = $1;
        my $db_ref  = $2;
        my $mod_id  = $3;
        if ( $db_ref eq $db_name ) {

            # IDs can be AT.G or At.g:
            if ( $mod_id !~ /\A AT \S G \d+ \z/ixms ) {
                die "Can't make sense of putative MOD ID $mod_id in input from MOD ID table 1 $id_table_1: $input\n";
            }
            # But for the FASTA headers which I will later be exploiting to get an accurate gene ID->human name map, *only* uppercase works, so:
            $mod_id =~ tr/[a-z]/[A-Z]/;

            if ( $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} ) {
                my @known_mod_ids = sort keys %{$data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} };
                push @known_mod_ids, $mod_id;
                my $known_mod_id_text = join ' vs. ', @known_mod_ids;
                warn "Multiple mapping of UniProt $uniprot to gene IDs $known_mod_id_text\n";
            }

            # Note: we have to do this, because there are instances of *identical* proteins encoded by two or more genes...
            $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'}->{$mod_id} = 1;

            # Add this, so that I have a filter later on for whether something is worth giving a full name.
            # We have to include the extra pointer to 'seen', or we crash the script later when it wants to point at 'long_name' from $mod_id instead.
            $data_ref->{'mod_gene_id'}->{$mod_id}->{'seen'} = 1;
        }
    }
    else { 
        die "Can't parse input from MOD ID table 1 $id_table_1: $input\n";
    }
}
close $ID1 or die "Can't close filehandle to MOD ID table 1 $id_table_1: $!";

open my $ID2, '<', $id_table_2 or die "Can't open MOD ID table 2 $id_table_2: $!";
while (my $input = <$ID2>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) \. (\d+) \s \| \s Symbols: (.*?) \s \| /xms ) { 
        # Sample input:
        # >AT1G73440.1 | Symbols:  | calmodulin-related | chr1:27611418-27612182 FORWARD LENGTH=254
        # >AT1G75120.1 | Symbols: RRA1 | Nucleotide-diphospho-sugar transferase family protein | chr1:28197022-28198656 REVERSE LENGTH=402
        # >AT3G10525.1 | Symbols: LGO, SMR1 | LOSS OF GIANT CELLS FROM ORGANS | chr3:3281576-3281962 REVERSE LENGTH=128

        my $mod_id      = $1;
        my $isoform_no  = $2;
        my $common_name = $3;

        # Trim off any *leading* \s, so they don't botch things up later:
        $common_name =~ s/\A\s+//g;
        my $long_name   = q{};

        # Insist on only scraping data out of the first isoform (.1) because other ones have inconsistent data!!
        if ( ( exists $data_ref->{'mod_gene_id'}->{$mod_id}->{'seen'} ) and ( $isoform_no == 1 ) ) {
            # Pick the *first* text for the canonical common name, if it's there at all -- trimming any trailing commas first.
            if ( $common_name =~ /\A (\S+) /xms ) { 
                $common_name = $1; 
                $common_name =~ s/,\z//; 
                $long_name   = $mod_id . q{|} . $common_name;
            }
            else {
                $long_name = $mod_id;
            }

            # Check this after we have dealt with dodgy regex behavior of the input.
            if (     ( exists $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} )
                 and ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ne $long_name ) ) {
                die "Inconsistent mapping of MOD ID $mod_id to long gene names $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} and $long_name\n";
            }

            $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} = $long_name;
        }
    }
    elsif ( $input =~ /\A > /xms ) {
        die "Can't parse input from MOD ID table 2 $id_table_2: $input\n";
    }
}
close $ID2 or die "Can't close filehandle to MOD ID table 2 $id_table_2: $!";

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };
my $uniprot_count = @uniprots;
foreach my $uniprot (@uniprots) {
    my @mod_ids = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} };
    foreach my $mod_id (@mod_ids) {
        if ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ) { 
            my $long_name = $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'};
            print "$uniprot\t$long_name\n";
        }
        else {
            warn "Can't find long name for uniprot $uniprot and MOD ID $mod_id\n";
            print "$uniprot\t$mod_id\n";
        }
    }
}

