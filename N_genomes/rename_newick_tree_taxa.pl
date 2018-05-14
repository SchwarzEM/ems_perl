#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $names  = q{};
my $newick = q{};

my %old2new = ();

$names  = $ARGV[0] if $ARGV[0];
$newick = $ARGV[1] if $ARGV[1];

if ( (! $names) or (! $newick) ) { 
    die "Format: rename_newick_tree_taxa.pl [name conversion table] [Newick text tree file]\n";
}

open my $NAMES, '<', $names;
while (my $input = <$NAMES>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $old = $1;
        my $new = $2;
        if ( exists $old2new{$old} ) { 
            die "Old name $old maps both to $old2new{$old} and to $new\n";
        }
        $old2new{$old} = $new;
    }
    else { 
        die "In name-conversion table $names, cannot parse: $input\n";
    }
}
close $NAMES;

open my $NEWICK, '<', $newick;
while (my $input = <$NEWICK>) {
    chomp $input;
    my $output = q{};
    while ( $input =~ /\A (.*?) (,|\() ([^:,\(\s]+) [:] (.+) /xmsg ) {
        my $leader   = $1;
        my $punct    = $2;
        my $name     = $3;
        my $trailer  = $4;
        my $new_name = update_name($name);
        $output = $output . $leader . $punct . $new_name . q{:};
        $input = $trailer;
    }
    $output = $output . $input;
    print "$output\n";
}
close $NEWICK;

sub update_name {
    my $_input  = $_[0];
    my $_output = q{};

    if (exists $old2new{$_input}) {
        $_output = $old2new{$_input};
    }
    elsif (! looks_like_number($_input)) {
        die "Cannot parse name: $_input\n";
    }
    else {
       $_output = $_input;
    }
    return $_output;
}
