#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $orig_name_line = q{};
my $new_name_line  = q{};
my $header_line    = q{};

$orig_name_line = $ARGV[0] if $ARGV[0];
$new_name_line  = $ARGV[1] if $ARGV[1];
$header_line    = $ARGV[2] if $ARGV[2];

if ( (! -e $orig_name_line ) or (! -e $new_name_line ) or (! -e $header_line ) ) {
    die "Format: rename_headers_16may2021.pl [orig name line, tabbed] [corresponding new name line, tabbed] [header line, tabbed] > [new header line, tabbed]\n";
}

my @orig_names         = ();
my @new_names          = ();
my @orig_header_fields = ();
my @new_header_fields  = ();
my %orig2new           = ();

open my $ORIG, '<', $orig_name_line;
while (my $input = <$ORIG>) {
    chomp $input;
    if ( @orig_names ) {
        die "Trying to list original column names twice\n";
    }
    @orig_names = split '\t', $input;
}
close $ORIG;

open my $NEW, '<', $new_name_line;

while (my $input = <$NEW>) {
    chomp $input;
    if ( @new_names ) {
        die "Trying to list new column names twice\n";
    }
    @new_names = split '\t', $input;
    my $orig_name_count = @orig_names;
    my $new_name_count  = @new_names;
    if ( $orig_name_count != $new_name_count ) {
        die "$orig_name_count original names, but $new_name_count new names; do not match\n";
    }
}
close $NEW;

my $j = @orig_names;
$j--;

foreach my $i (0..$j) {
    my $orig_name = $orig_names[$i];
    if ( $orig_name !~ /\S/xms ) {
        die "Do not have usable value for original name: $orig_name\n";
    }
    my $new_name  = $new_names[$i];
    if ( $new_name !~ /\S/xms ) {
        die "Do not have usable value for original name: $new_name\n";
    }
    $orig2new{$orig_name} = $new_name;
}

open my $HEADER, '<', $header_line;
while (my $input = <$HEADER>) {
    chomp $input;
    if ( @orig_header_fields ) {
        die "Trying to list current column names twice\n";
    }
    @orig_header_fields = split '\t', $input;
    foreach my $header_field (@orig_header_fields) {
        if ( $header_field !~ /\S/xms ) {
            die "Do not have usable value for header field: $header_field\n";
        }
        if ( exists $orig2new{$header_field} ) {
            $header_field = $orig2new{$header_field};
        }
        push @new_header_fields, $header_field;
    }
}
close $HEADER;

my $new_header_line = join "\t", @new_header_fields;
print "$new_header_line\n";


