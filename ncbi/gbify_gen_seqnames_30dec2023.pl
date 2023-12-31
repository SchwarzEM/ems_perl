#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $seq2gbacc = q{};
my $fasta     = q{};

$seq2gbacc = $ARGV[0] if $ARGV[0];
$fasta     = $ARGV[1] if $ARGV[1];

my %orig2gb_acc = ();

if ( (! $seq2gbacc ) or (! $fasta) ) {
    die "Format: gbify_gen_seqnames_30dec2023.pl [seq to gbacc table] [FASTA] > [GenBank-renamed FASTA] ;\n";
}

open my $SEQ2GBACC, '<', $seq2gbacc;
while (my $input = <$SEQ2GBACC>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $orig_name = $1;
        my $gb_acc    = $2;
        if (! exists $orig2gb_acc{$orig_name} ) {
            $orig2gb_acc{$orig_name} = $gb_acc;
        }
        else {
            die "In seq to gbacc table $seq2gbacc, \"$orig_name\" has redundant GenBank accessions \"$orig2gb_acc{$orig_name}\" and \"$gb_acc\"\n";
        }
    }
    else {
        die "From seq to gbacc table $seq2gbacc, cannot parse: $input\n";
    }
}
close $SEQ2GBACC;

open my $FASTA, '<', $fasta;
while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input !~ /\A [>] /xms ) {
        print "$input\n";
    }
    elsif ( $input =~ /\A [>] (\S+) (.*) /xms ) {
        my $orig_name = $1;
        my $comments  = $2;

        my $gb_acc = q{};
        if (! exists $orig2gb_acc{$orig_name} ) {
            die "In FASTA $fasta, cannot assign GenBank acc. to: $orig_name\n";
        }
        $gb_acc = $orig2gb_acc{$orig_name};

        my $acc_prefix = q{};
        if ( $gb_acc =~ /\A ([A-Z]{6}) \d+ /xms ) {
            $acc_prefix = $1;
        }
        else {
            die "Cannot extract accession prefix from \"$gb_acc\" derived from: $input\n";
        }

        my $new_name = 'gnl|WGS:' . $acc_prefix  . q{|} . $orig_name . '|gb|' . $gb_acc;
        print ">$new_name$comments\n";
    }
    else {
        die "From FASTA $fasta, cannot parse: $input\n";
    }
}
close $FASTA;

