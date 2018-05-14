#!/usr/bin/env perl

# link_Reinke_CDS_WS220.pl -- Erich Schwarz <emsch@caltech.edu>, 6/21/2012.
# Purpose: given existing files in wormpep220, link Reinke CDSes from approximately WS110 to the present.

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $data_ref;

my $wb220_hist = q{};
my $wpep110    = q{};
my $wpep220    = q{};
my $reinke     = q{};

my $cds        = q{};
my $cepep      = q{};
my $wbgene     = q{};
my $descriptor = q{};

($wb220_hist, $wpep110, $wpep220, $reinke) = @ARGV;

# Identify which CEpeps link to which CDS gene names, and which of those gene names are current:

open my $WB220_HIST, '<', $wb220_hist or die "Can't open wormpep220 history ($wb220_hist): $!";
while (my $input = <$WB220_HIST>) {
    chomp $input;

    # Sample input:
    # AH6.12  CE01453 8
    # AH6.13  CE01454 8       12
    # ZK593.6	CE06627	11	220
    # ZK593.6a	CE06627	220
    # ZK593.6b	CE45416	220

    # CDS sequence forms get converted to CDS gene names.

    # If a CEpep has a starting number but no ending number, it is both linked to a CDS and is *live*:
    if ( $input =~ /\A (\S+) \t (\S+) \t \d+ \s* \z/xms ) {
        $cds   = $1;
        $cepep = $2;
        $cds   =~ s/[a-z]+\z//;
        if ( $cds =~ /\d\z/xms ) { 
            $data_ref->{'cepep'}->{$cepep}->{'cds'}->{$cds} = 1;
            $data_ref->{'live_cds'}->{$cds} = 1;
        }
    }
    # If it has both starting and ending numbers, it is still linked to a CDS, but it is *archaic*.
    # I am assuming that this does not have outright name-changes of CDSes, because I can't cope with them!
    elsif ( $input =~ /\A (\S+) \t (\S+) \t \d+ \t \d+ \z/xms ) { 
        $cds   = $1;
        $cepep = $2;
        $cds   =~ s/[a-z]+\z//;
        if ( $cds =~ /\d\z/xms ) {
            $data_ref->{'cepep'}->{$cepep}->{'cds'}->{$cds} = 1;
        }
    }
    else { 
        die "Can't parse input line from wormpep220 history ($wb220_hist): $input\n";
    }
}
close $WB220_HIST or die "Can't close filehandle to wormpep220 history ($wb220_hist): $!";

# Identify which WS110 CDS gene names link to which WS110 CEpep names:

open my $WP110, '<', $wpep110 or die "Can't open wormpep110 ($wpep110): $!";
while (my $input = <$WP110>) { 
    chomp $input;

    # Sample input:
    # >2L52.1 CE32090   Zinc finger, C2H2 type status:Predicted TR:Q9XWB3 protein_id:CAA21776.2
    # >2RSSE.1 CE32785    status:Predicted TR:Q8I133 protein_id:CAD59137.1
    # >cTel54X.1 CE20619    status:Partially_confirmed TR:Q9XXA5 protein_id:CAA20224.1
    # >cTel7X.1 CE20618    status:Predicted TR:Q9XXA0 protein_id:CAA20336.1

    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (\S+) \s+ (CE\d+) /xms ) { 
            $cds   = $1;
            $cepep = $2;
            $cds   =~ s/[a-z]+\z//;
            if ( $cds =~ /\d\z/xms ) {
                $data_ref->{'ws110_cds'}->{$cds}->{'cepep'}->{$cepep} = 1;
            }
        }
        else { 
            die "From wormpep110 ($wpep110), can't parse header: $input\n";
        }
    }
}
close $WP110 or die "Can't close filehandle to wormpep110 ($wpep110): $!";

# Identify which WS220 CDS gene names link to which WS220 WBGene names:

open my $WP220, '<', $wpep220 or die "Can't open wormpep220 ($wpep220): $!";
while (my $input = <$WP220>) { 
    chomp $input;

    # Sample input:
    # >2RSSE.2        CE38260 WBGene00044165  status:Predicted        UniProt:A4F338  protein_id:ABO33280.1
    # >3R5.1  CE24758 WBGene00007065  locus:pot-3     status:Partially_confirmed      UniProt:Q9XWB2  protein_id:CAB76729.1
    # >4R79.1a        CE35820 WBGene00003525  locus:nas-6     Zinc-binding metalloprotease domain     status:Partially_confirmed      UniProt:Q9U3S9  protein_id:CAB63429.2
    # >4R79.1b        CE39659 WBGene00003525  locus:nas-6     status:Partially_confirmed      UniProt:Q2HQL9  protein_id:CAJ76926.1

    if ( $input =~ /\A > /xms ) {
        if ( $input =~ /\A > (\S+) \s+ CE\d+ \s+ (WBGene\d+) /xms ) {
            $cds = $1;
            $wbgene = $2;
            $cds   =~ s/[a-z]+\z//;
            if ( $cds =~ /\d\z/xms ) {
                $data_ref->{'ws220_cds'}->{$cds}->{'wbgene'} = $wbgene;
            }
        }
        else {
            die "From wormpep220 ($wpep220), can't parse header: $input\n";
        } 
    }
}
close $WP220 or die "Can't close filehandle to wormpep220 ($wpep220): $!";

# Get CDS gene names and descriptors for Reinke genes.
# Map as follows: WS110 CDS gene name -> WS110 CEpep -> WS220 CDS gene name from 'history' -> *live* WS220 CDS gene names -> WS220 WBGene name.
# Where the mapping succeeds, print out: WS220 WBGene name, descriptor.

open my $REINKE, '<', $reinke or die "Can't open Reinke sperm CDS/descriptor file ($reinke): $!";
while (my $input = <$REINKE>) {
    chomp $input;

    # Sample input:
    # Y71G12B.2	male sperm
    # Y73B6A.1	shared sperm
    # Y73B6A.2	herm sperm

    if ( $input =~ /\A (\S+) \t ([^\t]+) \s* \z/xms ) {
        my $ws110_cds = $1;
        $descriptor = $2;
        $ws110_cds =~ s/[a-z]+\z//;
        if ( exists $data_ref->{'ws110_cds'}->{$ws110_cds}->{'cepep'} ) { 
            my @cepeps = sort keys %{ $data_ref->{'ws110_cds'}->{$ws110_cds}->{'cepep'} };
            my @ws220_cdses = ();
            my @wbgenes     = ();
            foreach my $ws110_cepep (@cepeps) { 
                if ( exists $data_ref->{'cepep'}->{$ws110_cepep}->{'cds'} ) { 
                    my @tmp_cdses = grep { ( exists $data_ref->{'live_cds'}->{$_} ) } 
                                    keys %{ $data_ref->{'cepep'}->{$ws110_cepep}->{'cds'} }; 
                    push @ws220_cdses, @tmp_cdses;
                }
            }

            @ws220_cdses = sort @ws220_cdses;
            @ws220_cdses = uniq @ws220_cdses;

            if (@ws220_cdses) { 
                # Occam's razor: if there is just one CDS that matches the WS110 name, then stick with it;
                #    this has the pleasing effect of suppressing many cases where one cepep maps to 4-10 genes...
                my @subset_ws220_cdses = grep { $_ eq $ws110_cds } @ws220_cdses;
                if (@subset_ws220_cdses) { 
                    @ws220_cdses = @subset_ws220_cdses;
                }

                foreach my $ws220_cds (@ws220_cdses) { 
                    if ( exists $data_ref->{'ws220_cds'}->{$ws220_cds}->{'wbgene'} ) { 
                        my $tmp_wbgene = $data_ref->{'ws220_cds'}->{$ws220_cds}->{'wbgene'};
                        push @wbgenes, $tmp_wbgene;
                    }
                }
                if (@wbgenes) { 
                    @wbgenes = sort @wbgenes;
                    @wbgenes = uniq @wbgenes;
                    foreach my $wbgene (@wbgenes) { 
                        print "$wbgene\t$descriptor\n";
                    }
                }
            }
            if (! @ws220_cdses) {
                warn "From Reinke file $reinke, failed to detect any WS220 CDSes for ws110_cds $ws110_cds and descriptor \"$descriptor\"\n";
            }
        }
        else { 
            warn "From Reinke file $reinke, can't map ws110_cds $ws110_cds with descriptor \"$descriptor\" to any cepep\n";
        }
    }
    else { 
        die "Can't parse input from Reinke file $reinke: $input\n"
    }
}
close $REINKE or die "Can't close filehandle to Reinke sperm CDS/descriptor file ($reinke): $!";

