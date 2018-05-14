#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $phobius   = q{};
my $cons_unks = q{};

my $data_ref;

if ( (! $ARGV[0]) or (! $ARGV[1]) ) {
    die "Format: extend_fly_cons_unks_12apr2016.pl [rich Phobius summary] [cons. unk. table] > [STDOUT] ;\n";
}

$phobius   = $ARGV[0];
$cons_unks = $ARGV[1];

open my $PHOBIUS, '<', $phobius;
while (my $input = <$PHOBIUS>) {
    chomp $input;
    if ( $input =~ /\A (FBgn\d+) \S* \t ([^\t]+) \z/xms ) { 
        my $gene_id = $1;
        my $phobius = $2;
        $data_ref->{'gene_id'}->{$gene_id} = $phobius;
    }
}
close $PHOBIUS;

open my $CONS_UNKS, '<', $cons_unks;
while (my $input = <$CONS_UNKS>) {
    chomp $input;

    # Edit the original text to make plain 'Phobius' more clearly read 'Worm_Phobius'.
    $input =~ s/\tPhobius/\tWorm_Phobius/g;

    if ( $input =~ /\A Gene [ ] family /xms ) { 
        print "$input\tFly_Phobius\n";
    }
    else { 
        my @gene_ids = ();
        while ( $input =~ /(FBgn\d+)/xmsg ) { 
            my $gene_id = $1;
            if ( exists $data_ref->{'gene_id'}->{$gene_id} ) { 
                push @gene_ids, $gene_id;
            }
        }
        @gene_ids = sort @gene_ids;
        @gene_ids = uniq @gene_ids;
        my $phobius_annot  = q{};
        my @phobius_annots = ();
        if (@gene_ids) {
            @phobius_annots = grep { /\S/ } map { $data_ref->{'gene_id'}->{$_} } @gene_ids ;
        }
        @phobius_annots = sort @phobius_annots;
        @phobius_annots = uniq @phobius_annots;
        if (@phobius_annots) {
            $phobius_annot = join '; ', @phobius_annots;
            $phobius_annot = "Fly_Phobius: $phobius_annot";
        }
        print "$input\t$phobius_annot\n";
    }
}
close $CONS_UNKS;

