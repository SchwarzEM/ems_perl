#!/usr/bin/env perl

# get_cds_from_rpkm.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/21/2008.
# Purpose: get +/- CDSes sorted list from +/- (filtered) rpkm.

use strict;
use warnings;

my $DEBUG       = 0;
my $require_all = 0;

my %wanted_CDS = ();

my $input_args = join " ", @ARGV;
if ($input_args !~ /\A --(?:plus|all) \s .+ \s --minus \s .+  \z/xms ) { 
    die "Format: ./get_cds_from_rpkm.pl --plus|all [1+ rpkms]",
        " --minus [1+ rpkms]\n",
        ;
}
if ($input_args =~ /\A --(plus|all) \s (.+) \s --minus \s (.+)  \z/xms ) { 
    my $condition  = $1;
    my $plus_args  = $2;
    my $minus_args = $3;
    if ($condition eq 'all') { 
        $require_all = 1;
    }
    if ( ($condition ne 'all') and ($condition ne 'plus') ) { 
        die "Incorrect argument (\"$condition\")\n";
    }
    my @plus_files  = split /\s+/, $plus_args;
    my @minus_files = split /\s+/, $minus_args;
    foreach my $putative_file (@plus_files, @minus_files) { 
        if (! -e $putative_file ) { 
            die "Can't find putative file $putative_file!\n";
        }
        if ( ($DEBUG) and (-e $putative_file) ) {
            warn "Will try to read putative file $putative_file!\n";
        }
    }
    record_cdses(\@plus_files, \%wanted_CDS, $require_all);
    delete_cdses(\@minus_files, \%wanted_CDS);
    foreach my $specific_cds (sort keys %wanted_CDS) { 
        print "$specific_cds\n";
    }
}

sub record_cdses { 
    my $plus_files_ref   = $_[0];
    my $hash_to_fill_ref = $_[1];
    my $enforce_in_all   = $_[2];
    my $number_of_files  = @{ $plus_files_ref };
    my $seen_ref;
    foreach my $input_file ( @{ $plus_files_ref } ) { 
        open my $INFILE, '<', $input_file
            or die "Can't open input file $input_file: $!";
        while (my $input = <$INFILE>) { 
            if ($input =~ /\A (\S+) \s+ \S+ \s+ (\S+) /xms) { 
                my $cds  = $1;
                my $rpkm = $2;
                if ( ($rpkm ne '0.00') 
                     and (! $seen_ref->{$input_file}->{$cds} ) ) { 
                    $hash_to_fill_ref->{$cds} += 1;
                    $seen_ref->{$input_file}->{$cds} = 1;
                }
            }
        }
        close $INFILE or die "Can't close filehandle for $input_file: $!";
    }
    if ( $enforce_in_all == 1 ) {
        foreach my $seen_cds ( sort keys %{ $hash_to_fill_ref } ) { 
            if ( $hash_to_fill_ref->{$seen_cds} != $number_of_files ) { 
                delete $hash_to_fill_ref->{$seen_cds};
            }
        }
    }
    return;
}

sub delete_cdses { 
    my $minus_files_ref  = $_[0];
    my $hash_to_weed_ref = $_[1];
    foreach my $input_file ( @{ $minus_files_ref } ) {
        open my $INFILE, '<', $input_file
            or die "Can't open input file $input_file: $!";
        while (my $input = <$INFILE>) {
            if ($input =~ /\A (\S+) \s+ \S+ \s+ (\S+) /xms) {
                my $cds  = $1;
                my $rpkm = $2;
                if ( ($rpkm ne '0.00') 
                     and (exists $hash_to_weed_ref->{$cds} ) ) {
                    delete $hash_to_weed_ref->{$cds};   
                }
            }
        }  
        close $INFILE or die "Can't close filehandle for $input_file: $!";
    }
    return;
}

