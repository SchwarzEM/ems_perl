#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $data_ref;

my $mot_file      = $ARGV[0];
my $mot_prefix    = $ARGV[1];

my $motif_no      = q{};
my $scan_pfm_vals = 0;
my $matrix_id     = q{};
my $i             = 0;

if ( $mot_prefix !~ /\A (1|2) \z/xms ) { 
    die "Motif prefix must be 1 or 2, not \"$mot_prefix\"\n";
}

open my $MOT_FILE, '<', $mot_file or die "Can't open motif file $mot_file: $!";
while (my $input = <$MOT_FILE>) { 
    chomp $input;
    if ( $input =~ /\A A \s+ (\S+) \s+ C \s+ (\S+) \s+ G \s+ (\S+) \s+ T \s+ (\S+) \s+ \z/xms ) { 
        my $a_bg = $1;
        my $c_bg = $2;
        my $g_bg = $3;
        my $t_bg = $4;

        my @bg_vals = ( $a_bg, $c_bg, $g_bg, $t_bg );
        test_for_numerical(\@bg_vals);

        if ( exists $data_ref->{'bg'} ) {
            die "Redundant observation of putative background values from: $input\n";
        }

        $data_ref->{'bg'}->{'A'} = $a_bg;
        $data_ref->{'bg'}->{'C'} = $c_bg;
        $data_ref->{'bg'}->{'G'} = $g_bg;
        $data_ref->{'bg'}->{'T'} = $t_bg;
    }
    elsif ( $input =~ /\A MOTIF \s+ (\d+) \s+ .+ \s sites \s+ [=] \s+ (\d+) .+ E-value \s+ [=] \s+ (\S+) \s* \z/xms ) {
        $motif_no      = $1;
        my $site_count = $2;
        my $e_value    = $3;
        $scan_pfm_vals = 0;

        $i++;
        $matrix_id = sprintf "%08i", $i;

        my $brief_id  = "$mot_prefix-" . $motif_no;
        my $header    = 'Position_Matrix : "WBPmat' . $matrix_id . "\"";
        $site_count   = "Sites_used  $site_count";

        $data_ref->{'matrix'}->{$matrix_id}->{'header'}     = $header;
        $data_ref->{'matrix'}->{$matrix_id}->{'brief_id'}   = $brief_id;
        $data_ref->{'matrix'}->{$matrix_id}->{'site_count'} = $site_count;
        $data_ref->{'matrix'}->{$matrix_id}->{'e_value'}    = $e_value;
    }
    elsif ( $input =~ /\A letter-probability \s+ matrix:/xms ) { 
        $scan_pfm_vals = 1;
    }
    elsif ($scan_pfm_vals) { 
        if ( $input =~ /\A \s+ (\S+) \s+ (\S+) \s+ (\S+) \s+ (\S+) \s* \z/xms ) { 
            my $a_val = $1;
            my $c_val = $2;
            my $g_val = $3;
            my $t_val = $4;

            my @site_vals = ( $a_val, $c_val, $g_val, $t_val, );
            test_for_numerical(\@site_vals);

            push @{ $data_ref->{'matrix'}->{$matrix_id}->{'vals'}->{'A'} }, $a_val;
            push @{ $data_ref->{'matrix'}->{$matrix_id}->{'vals'}->{'C'} }, $c_val;
            push @{ $data_ref->{'matrix'}->{$matrix_id}->{'vals'}->{'G'} }, $g_val;
            push @{ $data_ref->{'matrix'}->{$matrix_id}->{'vals'}->{'T'} }, $t_val;
        }
        else { 
            $scan_pfm_vals = 0;
        }
    }
}
close $MOT_FILE or die "Can't close filehandle to motif file $mot_file: $!";

print "\n";

my @matrices = sort keys %{ $data_ref->{'matrix'} };
my @residues = qw( A C G T );
foreach my $matrix (@matrices) { 
    print "$data_ref->{'matrix'}->{$matrix}->{'header'}\n";
    my $brief_id = $data_ref->{'matrix'}->{$matrix}->{'brief_id'};
    my $brief_id_text = 'Brief_id      "Mortazavi_2010.' . $brief_id . '.pfm"';
    print "$brief_id_text\n";
    print "Description   \"Motif number $brief_id predicted by Mortazavi et al. (2010).\" Paper_evidence \"WBPaper00037732\" // pmid20980554\n";
    print "Type          Frequency\n";
    print "Consensus     \"$brief_id\"\n";
    foreach my $residue (@residues) {
        print "Background_model $residue $data_ref->{'bg'}->{$residue}\n";
    }
    foreach my $residue (@residues) {
        print "Site_values $residue @{ $data_ref->{'matrix'}->{$matrix}->{'vals'}->{$residue} }\n";
    }
    print "$data_ref->{'matrix'}->{$matrix}->{'site_count'}\n";
    print "Remark      \"Non-coding DNA motif conserved between C. elegans, C. briggsae, and C. angaria; e-value $data_ref->{'matrix'}->{$matrix}->{'e_value'}.\" Paper_evidence \"WBPaper00037732\" // pmid20980554\n";
    print "\n";
}

sub test_for_numerical { 
    my @_input_values = @{ $_[0] };
    foreach my $_input_val (@_input_values) { 
        if (! looks_like_number($_input_val) ) { 
            die "Failed to get numerical background values from: @_input_values\n";
        }
    }
    return;
}

