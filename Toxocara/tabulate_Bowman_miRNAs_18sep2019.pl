#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $tag1    = q{};
my $infile1 = q{};

my $tag2    = q{};
my $infile2 = q{};

my $data_ref;

$tag1    = $ARGV[0] if $ARGV[0];
$infile1 = $ARGV[1] if $ARGV[1];
$tag2    = $ARGV[2] if $ARGV[2];
$infile2 = $ARGV[3] if $ARGV[3];

if ( (! $tag1) or (! -r $infile1) or (! $tag2) or (! -r $infile2) ) {
    die "Format: tabulate_Bowman_miRNAs_18sep2019.pl [tag1] [infile1] [tag2] [infile2]\n";
}

if ( $tag1 !~/\A \S+ \z/xms ) {
    die "tag1 is not all-text: $tag1\n";
}
if ( $tag2 !~/\A \S+ \z/xms ) {
    die	"tag2 is not all-text: $tag2\n";
}
if ( $tag1 eq $tag2 ) {
    die "tag1 and tag2 are the same: $tag1, $tag2\n";
}

open my $INFILE1, '<', $infile1;
while (my $input = <$INFILE1>) {
    # Sample lines:
    # Name    Length  EffectiveLength TPM     NumReads
    # male.t0001910   23	3.000   3.081405        2.632
    # male.t0002066   24	3.000   0.000000        0.000
    # male.t0004361   21	3.000   1438.767171     1229.153
    chomp $input;
    if ( ( $input !~ /\AName/xms ) and ( $input =~ /\A (\S+) \s+ \d+ \s+ \S+ \s+ (\S+) \s+ (\S+) \z/xms ) ) {
        my $tx    = $1;
        my $tpm   = $2;
        my $reads = $3;
        $data_ref->{'tx'}->{$tx}->{'tag'}->{$tag1}->{'TPM'}   = $tpm;
       	$data_ref->{'tx'}->{$tx}->{'tag'}->{$tag1}->{'reads'} = $reads;
    }
    elsif ( $input !~ /\AName/xms ) {
        die "In input 1 file $infile1, cannot format: $input\n";
    }
}
close $INFILE1;

open my $INFILE2, '<', $infile2;
while (my $input = <$INFILE2>) {
    chomp $input;
    if ( ( $input !~ /\AName/xms ) and ( $input =~ /\A (\S+) \s+ \d+ \s+ \S+ \s+ (\S+) \s+ (\S+) \z/xms ) ) {
        my $tx    = $1;
        my $tpm   = $2;
        my $reads = $3;
       	$data_ref->{'tx'}->{$tx}->{'tag'}->{$tag2}->{'TPM'}   = $tpm;
        $data_ref->{'tx'}->{$tx}->{'tag'}->{$tag2}->{'reads'} = $reads;
    }
    elsif ( $input !~ /\AName/xms ) {
	die "In input 1 file $infile2, cannot format: $input\n";
    }
}
close $INFILE2;

my $header = "Transcript\tTPM.$tag1\tTPM.$tag2\treads.$tag1\treads.$tag2\n";

my @txs = sort keys %{ $data_ref->{'tx'} };
foreach my $tx1 (@txs) {
    print $header if $header;
    $header = q{};
    my $tpm_tag1 = $data_ref->{'tx'}->{$tx1}->{'tag'}->{$tag1}->{'TPM'};
    my $tpm_tag2 = $data_ref->{'tx'}->{$tx1}->{'tag'}->{$tag2}->{'TPM'};
    my $rds_tag1 = $data_ref->{'tx'}->{$tx1}->{'tag'}->{$tag1}->{'reads'};
    my $rds_tag2 = $data_ref->{'tx'}->{$tx1}->{'tag'}->{$tag2}->{'reads'};

    if ( ( $tpm_tag1 > 0 ) or ( $tpm_tag2 > 0 ) ) {
        print "$tx1\t$tpm_tag1\t$tpm_tag2\t$rds_tag1\t$rds_tag2\n";
    }

    if ( ( $tpm_tag1 > 0 ) and ( $tpm_tag2 > 0 ) ) {
        $data_ref->{'both_txs'}++;
    }
    elsif ( $tpm_tag1 > 0 ) {
        $data_ref->{'tag1_txs'}++;
    }
    elsif ( $tpm_tag2 > 0 ) {
       	$data_ref->{'tag2_txs'}++;
    }
}

$data_ref->{'both_txs'} = commify($data_ref->{'both_txs'});
$data_ref->{'tag1_txs'} = commify($data_ref->{'tag1_txs'});
$data_ref->{'tag2_txs'} = commify($data_ref->{'tag2_txs'});

# Print to STDERR so that this doesn't go into the main TSV output.
warn "\n";
warn "Total $tag1/$tag2 txs: $data_ref->{'both_txs'}\n";
warn "Total $tag1 only txs: $data_ref->{'tag1_txs'}\n";
warn "Total $tag2 only txs: $data_ref->{'tag2_txs'}\n";
warn "\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}


