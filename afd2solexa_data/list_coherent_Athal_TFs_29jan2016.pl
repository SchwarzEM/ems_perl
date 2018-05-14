#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $AtTFDB_file    = $ARGV[0];
my $PlantTFDB_file = $ARGV[1];
my $header         = "Gene\tTF";

my @genes = ();

open my $DB1, '<', $AtTFDB_file;
while (my $input = <$DB1>) {
    chomp $input;
    # sample input line:
    # ABI3VP1 At3g24650       ABI3    abscisic acid-insensitive protein 3 (ABI3)      VIEW    NA      NA      NA      No
    if ( $input =~ /\A ([^\t]+) \t (At (?: 1|2|3|4|5|c|m) g\d+) \t /xms ) {
        my $tf_type = $1;
        my $gene    = $2;
        $gene =~ tr/[a-z]/[A-Z]/;
        $data_ref->{'gene'}->{$gene}->{'tf_type'}->{$tf_type} = 1;
    }
}
close $DB1;

open my $DB2, '<', $PlantTFDB_file;
while (my $input = <$DB2>) {
    chomp $input;
    # sample input lines:
    # gene_model      plantTFDB_id    family  data_source
    # AT3G25730.1     AT3G25730.1     RAV     TAIR
    if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+)\.\d+ \t [^\t]* \t ([^\t]+) \t /xms ) {
        my $gene    = $1;
        my $tf_type = $2;
        $data_ref->{'gene'}->{$gene}->{'tf_type'}->{$tf_type} = 1;
    }
}
close $DB2;

if ( exists $data_ref->{'gene'} ) {
    @genes = sort keys %{ $data_ref->{'gene'} };
    print "$header\n";
    foreach my $gene (@genes) {
        my @tf_types     = sort keys %{ $data_ref->{'gene'}->{$gene}->{'tf_type'} };
        my $tf_type_text = join '; ', @tf_types;
        $tf_type_text    = "TF[$tf_type_text]";
        print "$gene\t$tf_type_text\n";
    }
}

