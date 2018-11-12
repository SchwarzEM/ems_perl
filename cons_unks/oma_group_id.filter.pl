#!/usr/bin/env perl

# oma_group_id.filter.pl -- Erich Schwarz <ems394@cornell.edu>, 11/12/2018.
# Purpose: given oma-groups.txt from OMA, and given a selected set of IDs, print *only* the OMA groups that contain 2+ preferred IDs.
# Note: by default, this script does not bother printing 'groups' that contain only one ID.  That default may become overridable if there is need.
# The key value of this script is reducing the huge and copious oma-groups.txt file to something of focused use for selected taxa.

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

my $id_list = q{};
my $oma     = q{};

$id_list = $ARGV[0] if $ARGV[0];
$oma     = $ARGV[1] if $ARGV[1];

if ( (! $id_list) or (! $oma) ) {
    die "Format: oma_group_id.filter.pl [ok_id_list] [oma_groups] > [oma_groups_w_2+_ok_ids]\n";
}

open my $IDS, '<', $id_list;
while (my $oma_id = <$IDS>) {
    chomp $oma_id;
    if ( $oma_id =~ /\A \S+ \z/xms ) {
        $data_ref->{'ok_id'}->{$oma_id} = 1;
    }
    else {
        die "In OMA ID set file $id_list, cannot parse: $oma_id\n";
    }
}
close $IDS;

open my $OMA, '<', $oma;
while (my $input = <$OMA> ) {
    chomp $input;
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A (\d+ \t \S+) \t (.+) \z/xms  ) ) {
        my $header      = $1;
        my $oma_id_text = $2;
        my @oma_ids     = split /\t/, $oma_id_text;
        my @ok_oma_ids  = ();

        foreach my $oma_id (@oma_ids) {
            if ( exists $data_ref->{'ok_id'}->{$oma_id} ) {
                push @ok_oma_ids, $oma_id;
            }
        }
        @ok_oma_ids = sort @ok_oma_ids;
        @ok_oma_ids = uniq(@ok_oma_ids);
        my $ok_oma_id_no = @ok_oma_ids;
        if ( $ok_oma_id_no >= 2 ) {
             my $ok_oma_id_text = join "\t", @ok_oma_ids;
             print "$header\t$ok_oma_id_text\n";
        }
    }
    elsif ( $input !~ /\A [#] /xms ) {
        die "In OMA family file $oma, cannot parse: $input\n";
   }
}
close $OMA;

