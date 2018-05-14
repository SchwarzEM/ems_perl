#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $infile        = q{};
my @taxa_w_ranges = ();

my $data_ref;

my $help;

GetOptions ( 'infiles=s'          => \$infile,
             'taxa_w_ranges=s{,}' => \@taxa_w_ranges,
             'help'               => \$help,   );

if ( $help or (! $infile) or (! @taxa_w_ranges) ) { 
    die "Format: count_orthofinder_taxa_15apr2017.pl\n",
        "    --infile|-i         <input stream/file>\n",
        "    --taxa_w_ranges|-t  [one or more of: non-redundant taxon; lower count of allowable occurrences; upper count o.a.o.]\n",
        "    --help|-h           [print this message]\n",
        ;
}

my $number_taxa_w_ranges = @taxa_w_ranges ;

# $number_taxa_w_ranges must be an integral multiple of 3
if ( ( $number_taxa_w_ranges % 3 ) != 0 ) {
    die "The number of arguments must be a multiple of 3: taxon, lower bound, upper bound; next taxon, etc.\n"
}

$number_taxa_w_ranges = ($number_taxa_w_ranges/3);
$number_taxa_w_ranges--;

foreach my $i (0..($number_taxa_w_ranges)) {
    my $taxon       = $taxa_w_ranges[ ($i*3) ];
    if ( exists $data_ref->{'taxon'}->{$taxon} ) {
        die "Redundant taxon \"$taxon\" in $infile\n";
    }

    my $lower_count = $taxa_w_ranges[ ( ($i*3) +1 ) ]; 
    if ($lower_count !~ /\A \d+ \z/xms ) { 
        die "Lower count does not look like number: \"$taxon $lower_count\" in @taxa_w_ranges\n";
    }

    my $upper_count = $taxa_w_ranges[ ( ($i*3) +2 ) ];
    if ($upper_count !~ /\A \d+ \z/xms ) {
        die "Upper count does not look like number: \"$taxon $lower_count $upper_count\" in @taxa_w_ranges\n";
    }

    if ( $lower_count > $upper_count ) {
        die "Upper count must be greater than or equal to lower count: \"$taxon $lower_count $upper_count\" in @taxa_w_ranges\n";
    }

    $data_ref->{'taxon'}->{$taxon}->{'lower_count'} = $lower_count;
    $data_ref->{'taxon'}->{$taxon}->{'upper_count'} = $upper_count;
}

my @taxa = sort keys %{ $data_ref->{'taxon'} };

my $INPUT_FILE;
if ($infile eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile;
}

while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    my $print_line = 1;

    foreach my $taxon (@taxa) {
        my $taxon_count = 0;

        my $lower_count = $data_ref->{'taxon'}->{$taxon}->{'lower_count'};
        my $upper_count = $data_ref->{'taxon'}->{$taxon}->{'upper_count'};

        while ( $input =~ /\( $taxon \)/xmsg ) { 
            $taxon_count++;
        }
        if ( ( $taxon_count < $lower_count ) or ( $upper_count < $taxon_count ) ) {
            $print_line = 0;
        }
    }

    if ( $print_line == 1 ) { 
        print "$input\n";
    }
}

