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

my $sig_list = q{};
my $oma      = q{};

$sig_list = $ARGV[0] if $ARGV[0];
$oma      = $ARGV[1] if $ARGV[1];

if ( (! $sig_list) or (! $oma) ) {
    die "Format: oma_group_id.filter.pl [ok_signature_list] [oma_groups] > [oma_groups_w_2+_ok_sigs]\n";
}

open my $SIGS, '<', $sig_list;
while (my $oma_sig = <$SIGS>) {
    chomp $oma_sig;
    if ( $oma_sig =~ /\A \S+ \z/xms ) {
        $data_ref->{'ok_sig'}->{$oma_sig} = 1;
    }
    else {
        die "In OMA signature set file $sig_list, cannot parse: $oma_sig\n";
    }
}
close $SIGS;

open my $OMA, '<', $oma;
while (my $input = <$OMA> ) {
    chomp $input;
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A \d+ \t (\S+) \t .+ \z/xms  ) ) {
        my $oma_sig = $1;
        if ( exists $data_ref->{'ok_sig'}->{$oma_sig} ) {
            print "$input\n";
        }
    }
    elsif ( $input !~ /\A [#] /xms ) {
        die "In OMA family file $oma, cannot parse: $input\n";
   }
}
close $OMA;

