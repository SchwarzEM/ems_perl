#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $wbgene2go_file = $ARGV[0];
my $go2term_file   = $ARGV[1];

my $header = "Gene\tGO_term";

open my $WB, '<', $wbgene2go_file;
while (my $input = <$WB>) {
    chomp $input;
    if ( $input =~ /\A WB \t (WBGene\d+) \t [^\t]* \t [^\t]* \t (GO:\d+) \t /xms ) {  
        my $wbgene = $1;
        my $go_id  = $2;
        $data_ref->{'wbgene'}->{$wbgene}->{'go_id'}->{$go_id} = 1;
    }
    elsif ( $input !~ /\A [!] /xms ) { 
        die "From WBGene-to-GO file $wbgene2go_file, cannot parse: $input\n";
    }
}
close $WB;

open my $GO, '<', $go2term_file;
while (my $input = <$GO>) {
    chomp $input;
    # sample input line:
    # 27	mitochondrion inheritance	biological_process	GO:0000001	0	0	0
    if ( $input =~ /\A \d+ \t ([^\t]+) \t [^\t]+ \t (GO:\d+) \t /xms ) {
        my $go_text = $1;
        my $go_id   = $2;
        my $go_desc = "$go_text [$go_id]";
        if ( exists $data_ref->{'go_id'}->{$go_id}->{'go_desc'} ) {
            die "Redundant descriptions of GO ID $go_id: $go_desc and $data_ref->{'go_id'}->{$go_id}->{'go_desc'}\n";
        }
        $data_ref->{'go_id'}->{$go_id}->{'go_desc'} = $go_desc;
    }
}
close $GO;

my @wbgenes = sort keys %{ $data_ref->{'wbgene'} };
foreach my $wbgene (@wbgenes) {
    my @go_descs = sort
                   map { go_id_to_desc($_) } 
                   keys %{ $data_ref->{'wbgene'}->{$wbgene}->{'go_id'} };
    my $go_text  = join '; ', @go_descs;
    print "$header\n" if $header;
    $header = q{};
    print "$wbgene\t$go_text\n";
}

sub go_id_to_desc {
    my $_go_id = $_[0];
    if ( exists $data_ref->{'go_id'}->{$_go_id}->{'go_desc'} ) {
        return $data_ref->{'go_id'}->{$_go_id}->{'go_desc'};
    }
    else {
        return;
    }
}

