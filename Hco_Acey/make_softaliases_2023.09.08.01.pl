#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $orig_file_list = q{};
my $file2dir_list  = q{};

$orig_file_list = $ARGV[0] if $ARGV[0];
$file2dir_list  = $ARGV[1] if $ARGV[1];

my $data_ref;

if ( (! $orig_file_list ) or (! $file2dir_list ) ) {
    die "Format: make_softaliases_2023.09.08.01.pl [orig full file list] [file to dir list] > [softaliasing commands]\n";
}

open my $ORIG_FILES, '<', $orig_file_list;
while ( my $orig_file = <$ORIG_FILES> ) {
    chomp $orig_file;
    if ( $orig_file =~ /\A \/ \S+ \/ ([^\s\/]+) \z/xms ) {
        my $file = $1;
        if ( exists $data_ref->{'file'}->{$file}->{'orig_file'} ) {
            die "Redundant mapping of $orig_file to $file\n";
        }
        $data_ref->{'file'}->{$file}->{'orig_file'} = $orig_file;
    }
    else {
        die "From orig. file list $orig_file_list, cannot parse: $orig_file\n";
    }
}
close $ORIG_FILES;

open my $FILE2DIR, '<', $file2dir_list;
while ( my $input = <$FILE2DIR> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $file = $1;
        my $dir  = $2;
        if (! $data_ref->{'file'}->{$file}->{'orig_file'} ) {
            die "Cannot map file $file to orig. file\n";
        }
        my $orig_file = $data_ref->{'file'}->{$file}->{'orig_file'};
        print "ln -s $orig_file $dir/$file ;\n";
    }
    else {
        die "From file2dir list $file2dir_list, cannot parse: $input\n";
    }
}
close $FILE2DIR;
