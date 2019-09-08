#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $format   = q{};
my $proteome = q{};
my $gene     = q{};
my $tx       = q{};

my $data_ref;
my $max;

my %ok_formats = (
    'augustus'    => 'aug',
    'aug'         => 'aug',
    'ensembl'     => 'ens',
    'ens'         => 'ens',
    'ra'          => 'ra',
    'arabidopsis' => 'tair',
    'tair'        => 'tair',
    'wormbase'    => 'wb',
    'wb'          => 'wb',
);

my $help;

GetOptions ( 'format=s'   => \$format,
             'proteome=s' => \$proteome, 
             'max'        => \$max,
             'help'       => \$help, );

if ($help or (! $proteome) or (! exists $ok_formats{$format}) ) { 
    die "Format: tabulate_proteome_sizes.pl\n",
        "    --format|-f   [proteome header format -- acceptable format names include:\n",
        "                  'ensembl' or 'ens'\n",
       	"                  'wormbase' or 'wb'\n",
       	"                  'augustus' or 'aug'\n",
       	"                  'arabidopsis' or 'tair'\n", 
       	"               or 'ra' for -RA-suffixed]\n",
        "    --proteome|-p [protein file; gets gene names, transcript names, and protein lengths]\n",
        "    --max|-m      [instead of human-readable range, give number of aa for largest isoform only]\n",
        "    --help|-h     [prints this message]\n",
        ;
}

# Map all format names to smaller subset:
$format = $ok_formats{$format};

open my $PROT, '<', $proteome or die "Can't open proteome: $proteome\n";
while (my $input = <$PROT>) { 
    if ( $input =~ /\A > /xms ) { 
        if ( ( $format eq 'aug' ) and ( $input =~ /\A > ((\S+\.g\d+)\.t\d+) /xms ) ) { 
            $tx   = $1;
            $gene = $2;
        }
	elsif ( ( $format eq 'ens' ) and ( $input =~ /\A > (\S+) \s .* gene:(\S+) /xms ) ) {
            $tx   = $1;
            $gene = $2;
        }
	elsif ( ( $format eq 'ra' ) and ( $input =~ /\A > ((\S+)\-R[A-Z]) \s /xms ) ) {
            $tx   = $1;
            $gene = $2;
        }
	elsif ( ( $format eq 'tair' ) and ( $input =~ /\A > ((\S+)\.\d+)\b /xms ) ) {
            $tx   = $1;
            $gene = $2;
        }
        elsif ( ( $format eq 'wb' ) and ( $input =~ /\A > (\S+) \s .* (WBGene\d+) /xms ) ) {
            $tx   = $1;
            $gene = $2;
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }
    elsif ( $input =~ /\S/xms ) { 
        $input =~ s/\s//;
        $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx}->{'aa_length'} += length($input);
    }
}
close $PROT or die "Can't close filehandle to proteome: $proteome\n";

my $header = "Gene\tProt_size\n";
if ($max) {
    $header = "Gene\tMax_prot_size\n";
}

my @genes = grep { /\S/ } sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @txs     = grep { /\S/ } sort keys %{ $data_ref->{'gene'}->{$gene1}->{'tx'} };
    my @lengths = ();
    foreach my $tx1 (@txs) { 
        my $length = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'aa_length'};
        push @lengths, $length;
    }
    @lengths = sort { $a <=> $b } @lengths;
    my $shortest = $lengths[0];
    my $longest  = $lengths[-1];
    my $length_text = q{};
    if (! $max) { 
        if ( $shortest == $longest ) { 
            $length_text = "\"$shortest aa\"";
        }
        if ( $longest != $shortest ) { 
            $length_text = "\"$shortest-$longest aa\"";
        }
    }
    if ($max) { 
        $length_text = $longest;
    }
    print $header if $header;
    $header = q{};
    print "$gene1\t$length_text\n";
}

