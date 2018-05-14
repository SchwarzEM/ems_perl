#!/usr/bin/env perl

# differential_pfam_interpro.pl -- Erich Schwarz <ems394@cornell.edu>, 3/30/2013.
# Purpose: given tabular Hmmer 3.0/Pfam-A or raw InterProScan hits, identify motifs which are present in required, absent from banned, and possibly present in optional taxa.

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $format_input = q{};
my @required     = ();
my @banned       = ();
my @optional     = ();
my $thresh_req   = 1;
my $thresh_ban   = 1;
my $file2species = q{};

my $data_ref;

my $help;

GetOptions ( 'format=s'      => \$format_input,
             'required=s{,}' => \@required,
             'banned=s{,}'   => \@banned,
             'optional=s{,}' => \@optional,
             'thresh_req=f'  => \$thresh_req,
             'thresh_ban=f'  => \$thresh_ban,
             'species=s'     => \$file2species,
             'help'          => \$help, );

if ( $help or (! $format_input) or (! @required) ) { 
    die "Format: differential_pfam_interpro\n",
        "            --format|-f    [input format: 'pfam' or 'interpro|ipr']\n",
        "            --required|-r  [1+ tables in which motif must be always present]\n",
        "            --banned|-b    [1+ tables in which motif must never be present (not mandatory)]\n",
        "            --optional|-o  [1+ tables in which motif may be present (not mandatory)]\n",
        "            --thresh_req   [E-value maximum acceptable threshold for required or optional hits; default is 1]\n",
        "            --thresh_ban   [E-value maximum acceptable threshold for banned hits; default is 1]\n",
        "            --species|-s   [table linking file names to species names (not mandatory)]\n",
        "            --help|-h      [print this message]\n",
        ;
}

# Enforce and parse $format_input.
if ( ( $format_input ne 'pfam' ) and ( $format_input ne 'interpro' ) and ( $format_input ne 'ipr' ) ) {
    die "Format input --format|-f must be 'pfam' or 'interpro|ipr', not $format_input\n";
}
# Use full name from here on:
if ( $format_input eq 'ipr' ) {
    $format_input = 'interpro';
}

# Enforce sane values for $thresh_req and $thresh_ban.
my %thresholds = ( 'thres_req'  => $thresh_req,
                   'thresh_ban' => $thresh_ban, );
my @thres_keys = sort keys %thresholds;
foreach my $thresh_key (@thres_keys) {
    if ( (! looks_like_number($thresholds{$thresh_key}) ) or ( $thresholds{$thresh_key} < 0 ) ) { 
        die "The threshold value $thresh_key must be a non-negative number, not $thresholds{$thresh_key}\n";
    }
    if ( $thresholds{$thresh_key} > 1 ) {
        warn "The threshold value $thresh_key is dubiously large value of $thresholds{$thresh_key}\n";
    }
}

# Optionally give files short species names; otherwise, just let each file's 'species' be its base filename:
if ($file2species) { 
    open my $SPP, '<', $file2species or die "Can't open file-to-species table $file2species: $!";
    while (my $input = <$SPP>) { 
        chomp $input;
        if ( $input =~ /\A ([^\t]+) \t ([^\t]+) \z/xms ) { 
            my $file    = $1;
            my $species = $2;
            # Enforce unique naming.
            if ( exists $data_ref->{'species'}->{$species}) { 
                die "Two or more files both assigned to species $species\n";
            }
            $data_ref->{'species'}->{$species} = 1;
            $data_ref->{'file'}->{$file}->{'species'} = $species;
        }
        else { 
            die "From file-to-species table $file2species, can't parse input: $input\n";
        }
    }
    close $SPP or die "Can't close filehandle to file-to-species table: $!";
}
my @total_files = ();
push @total_files, @required, @banned, @optional;
foreach my $file (@total_files) { 
    my $basename = basename $file;
    if (! $data_ref->{'file'}->{$file}->{'species'}) {
        # Again, enforce unique naming.
        if ( exists $data_ref->{'species'}->{$basename} ) { 
            die "Two or more files both assigned to faux-species $basename\n";
        }
        $data_ref->{'species'}->{$basename} = 1; 
        $data_ref->{'file'}->{$file}->{'species'} = $basename;
    }
}

# Finally, start the real work.  Define the search pattern (which is dependent on what format we'll accept).
my $search_pattern = q{};

# Sample partial PFAM input:
#
# #                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
# #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
# BUX.s01109.415       -          1-cysPrx_C           PF10417.4    2.3e-13   49.8   0.1   4.2e-13   48.9   0.1   1.5   1   0   0   1   1   1   1 BUX.s01109.415  BUX.gene.s01109.415


# Note that we are *not* trying to quote this in such a way that it works well in regex searches. Why not?
# Source: Perl Cookbook, 2cd. ed., ch. 6, part 10.
# $search_pattern = qr/$search_pattern/;

# Note that this is not yet fully coded -- I need to put in the real case for when I'm scanning InterPro!
if ( $format_input eq 'pfam' ) { 
    $search_pattern = '\A (\S+) \s+ \S+ \s+ (\S+) \s+ (\S+) \s+ (\S+) \s';
                    # start of line      motif_name  accession  E-value
}

# Sample partial interpro RAW format:
#
# Gene3D  G3DSA:1.10.167.10       no description  0.001
# Gene3D  G3DSA:3.30.200.20       no description  7.2e-11
# Seg     seg     seg     NA
# PatternScan     PS00108 PROTEIN_KINASE_ST       NA
# HMMPfam PF01484 Col_cuticle_N   2e-16
# superfamily     SSF48097        Regulator of G-protein signaling, RGS   1.2e-13
# [etc.]

# Note that InterPro has at least some components whose 'descriptions' are empty, so allow that in this pattern:
if ( $format_input eq 'interpro' ) {
    $search_pattern = '\A ([^\t]+) \t ([^\t]+) \t ([^\t]*) \t ([^\t]+) \z';
        # Start of line ^  classif.   access.   description  p-/E-value ^ end of line
}

# Find banned motifs.
foreach my $banned_file (@banned) {
    open my $BAN, '<', $banned_file or die "Can't open banned file $banned_file: $!";
    while (my $input = <$BAN>) { 
        chomp $input;
        if ( $input =~ /$search_pattern/xms ) {
            if ( $format_input eq 'pfam' ) { 
                my $start      = $1;
                my $accession  = $3;
                my $e_value    = $4;

                if ( ( $start !~ /\A [#] /xms ) and ( $e_value <= $thresh_ban ) ) { 
                    $data_ref->{'accession'}->{$accession}->{'banned'} = 1;
                }
            }
            elsif ( $format_input eq 'interpro' ) {
                my $classification = $1;
                my $accession      = $2;
                my $description    = $3;
                my $e_value        = $4;

                # InterProScan can put out numbers like '1E-100', so throughout the script, silently correct that to '1e-100'.
                $e_value =~ s/E/e/g;
                my $full_accession = $classification . q{:} . $accession;

                if ( ( looks_like_number($e_value) ) and ( $e_value <= $thresh_ban ) ) { 
                    $data_ref->{'accession'}->{$full_accession}->{'banned'} = 1;
                }
            }
        }
        else { 
            die "In banned file $banned_file, can't parse input line: $input\n";
        }
    }
    close $BAN or die "Can't open banned file $banned_file: $!";
}

# Find required motifs.
foreach my $req_file (@required) {
    open my $REQ, '<', $req_file or die "Can't open required file $req_file: $!";
    while (my $input = <$REQ>) {
        chomp $input;
        if ( $input =~ /$search_pattern/xms ) {
            if ( $format_input eq 'pfam' ) {
                my $start      = $1;
                my $motif_name = $2;
                my $accession  = $3;
                my $e_value    = $4;

                if (     ( $start !~ /\A [#]/xms                                      ) 
                     and ( $e_value <= $thresh_req                                    ) 
                     and (! exists $data_ref->{'accession'}->{$accession}->{'banned'} ) ) { 
                    $data_ref->{'accession'}->{$accession}->{'req_file'}->{$req_file} = 1;
                    $data_ref->{'accession'}->{$accession}->{'motif_name'} = $motif_name;
                }
            }
            elsif ( $format_input eq 'interpro' ) {
                my $classification = $1;
                my $accession      = $2;
                my $description    = $3;
                my $e_value        = $4;

                $e_value =~ s/E/e/g;
                my $full_accession = $classification . q{:} . $accession;

                if (     ( looks_like_number($e_value)                                     )
                     and ( $e_value <= $thresh_req                                         ) 
                     and (! exists $data_ref->{'accession'}->{$full_accession}->{'banned'} ) ) {
                    $data_ref->{'accession'}->{$full_accession}->{'req_file'}->{$req_file} = 1;
                    $data_ref->{'accession'}->{$full_accession}->{'motif_name'} = $description;
                }
            }
        }
        else {
            die "In required file $req_file, can't parse input line: $input\n";
        }
    }
    close $REQ or die "Can't close filehandle to required file $req_file: $!";
}   

# Given all the possible hits, keep only those which were found in *each* required input file; ban all others.
my $req_count       = @required;
my @possible_motifs = grep { exists( $data_ref->{'accession'}->{$_}->{'req_file'} ) } keys %{ $data_ref->{'accession'} };
foreach my $possible_motif (@possible_motifs) { 
    if ( exists $data_ref->{'accession'}->{$possible_motif}->{'req_file'} ) { 
        my @files_for_poss_motif = keys %{ $data_ref->{'accession'}->{$possible_motif}->{'req_file'} };
        my $count_for_poss_motif = @files_for_poss_motif;
        if ( $count_for_poss_motif < $req_count ) { 
            $data_ref->{'accession'}->{$possible_motif}->{'banned'} = 1;
            delete $data_ref->{'accession'}->{$possible_motif}->{'req_file'};
        }
    }
    else { 
        die "This shouldn't happen!  Obtuse Perl syntax wins again.\n";
    }
}

# Find the optional positive hits.
foreach my $opt_file (@optional) {
    open my $OPT, '<', $opt_file or die "Can't open optional file $opt_file: $!";
    while (my $input = <$OPT>) {
        chomp $input;
        if ( $input =~ /$search_pattern/xms ) {
            if ( $format_input eq 'pfam' ) {
                my $start      = $1;
                my $accession  = $3;
                my $e_value    = $4;
                if ( ( $start !~ /\A [#]/xms ) and ( $e_value <= $thresh_req ) and (! exists $data_ref->{'accession'}->{$accession}->{'banned'} ) ) {
                    $data_ref->{'accession'}->{$accession}->{'opt_file'}->{$opt_file} = 1;
                }
            }
            elsif ( $format_input eq 'interpro' ) {
                my $classification = $1;
                my $accession      = $2;
                my $description    = $3;
                my $e_value        = $4;

                $e_value =~ s/E/e/g;
                my $full_accession = $classification . q{:} . $accession;

                if ( ( looks_like_number($e_value) ) and ( $e_value <= $thresh_req ) ) {
                    $data_ref->{'accession'}->{$full_accession}->{'opt_file'}->{$opt_file} = 1;
                }
            }
        }
        else {
            die "In optional file $opt_file, can't parse input line: $input\n";
        }
    }
    close $OPT or die "Can't close filehandle to optional file $opt_file: $!";
}

# For each motif which somehow survived all this, get a table listing name, accession, required, banned, and optional-positive file/species.
my @final_motifs = grep { exists( $data_ref->{'accession'}->{$_}->{'req_file'} ) } keys %{ $data_ref->{'accession'} };
foreach my $final_motif (@final_motifs) {
    if ( exists $data_ref->{'accession'}->{$final_motif}->{'req_file'} ) {
        my $motif_name = $data_ref->{'accession'}->{$final_motif}->{'motif_name'};

        my @req_spp_for_final_motif = map { $data_ref->{'file'}->{$_}->{'species'} } @required;
        my $req_spp_txt             = join '; ', @req_spp_for_final_motif;

        my @banned_spp_for_final_motif = map { $data_ref->{'file'}->{$_}->{'species'} } @banned;
        my $banned_spp_txt             = join '; ', @banned_spp_for_final_motif;

        my @opt_spp_for_final_motif = keys %{ $data_ref->{'accession'}->{$final_motif}->{'opt_file'} };
        @opt_spp_for_final_motif    = map { $data_ref->{'file'}->{$_}->{'species'} } @opt_spp_for_final_motif;
        my $opt_spp_txt             = join '; ', @opt_spp_for_final_motif;

        print "$final_motif\t$motif_name\t$req_spp_txt\t$banned_spp_txt\t$opt_spp_txt\n";
    }
    else {
        die "This *really* shouldn't happen.  Ach, Perl!\n";
    }
}   

