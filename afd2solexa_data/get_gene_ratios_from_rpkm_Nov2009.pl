#!/usr/bin/env perl

# get_gene_ratios_from_rpkm_Nov2009.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/2008 w/ renaming 12/2/2009.
# Purpose: get list of [+]WBGenes with ratios of RPKM means to [-]WBGenes, from +/- (filtered) rpkm.
# NOTE: this is a permanent archiving of an older version used for preliminary work in 9/2008 and 12/2009.
#       I reviewed the code; it looks reliable but a bit klunkier than it needs to be.

use strict;
use warnings;

my $allow_zeroed_plusses  = 1;
my $allow_nonzero_minuses = 1;
my %wanted_CDS = ();

my $input_args = join " ", @ARGV;
if ($input_args !~ /\A --(?:plus|all) \s .+ \s --(?:minus|none) \s .+  \z/xms ) { 
    die "Format: ./get_cds_ratios_from_rpkm.pl --plus|all [1+ rpkms]",
        " --minus|none [1+ rpkms]\n",
        ;
}
if ($input_args =~ /\A --(plus|all) \s (.+) \s --(minus|none) \s (.+)  \z/xms ) { 
    my $condition1  = $1;
    my $plus_args   = $2;
    my $condition2  = $3;
    my $minus_args  = $4;

    if ($condition1 eq 'all') { 
        $allow_zeroed_plusses = 0;
    }
    if ($condition2 eq 'none') {   
        $allow_nonzero_minuses = 0;
    }
    if (    ( ($condition1 ne 'plus')  and ($condition1 ne 'all') )
         or ( ($condition2 ne 'minus') and ($condition2 ne 'none') ) ) { 
        die "Incorrect argument (\"$condition1\")\n";
    }

    if ( ($condition1 ne 'plus') and ($condition1 ne 'all') ) {
        die "Argument \"$condition1\" should be \"plus\" or \"all\"\n";
    }
    if ( ($condition2 ne 'minus') and ($condition2 ne 'none') ) {
        die "Argument \"$condition2\" should be \"minus\" or \"none\"\n";
    }

    my @plus_files  = split /\s+/, $plus_args;
    my @minus_files = split /\s+/, $minus_args;
    my %seen_plus_file = ();

    foreach my $plus_query (@plus_files) { 
        if (! -e $plus_query ) { 
            die "Can't find putative file $plus_query!\n";
        }
        $seen_plus_file{$plus_query} = 1;
    }

    foreach my $minus_query (@minus_files) { 
        if (! -e $minus_query ) {
            die "Can't find putative file $minus_query!\n";
        }
        if ( exists $seen_plus_file{$minus_query} ) { 
            die "File $minus_query is shared by positive",
                " and negative file sets!\n",
                ;
        }
    }

    record_positive_rpkms(\@plus_files, \%wanted_CDS, $allow_zeroed_plusses);
    generate_pos_neg_rpkm_ratios(\@minus_files, \%wanted_CDS, $allow_nonzero_minuses);

    my @specific_cdses =
        map $_->[0],
        sort { $b->[1] <=> $a->[1] } 
        map [ $_, $wanted_CDS{$_}->{'uniq_rpkm_ratio'} ],
        # W/o next grep, '-all' or '-none' filters destroy the mapping:
        grep { (! $wanted_CDS{$_}->{'delete'} ) }
        keys %wanted_CDS;

    foreach my $specific_cds (@specific_cdses) { 
        if ( $wanted_CDS{$specific_cds}->{'rpkm_ratio'} > 0 ) { 
            my $a = sprintf "%.2f", $wanted_CDS{$specific_cds}->{'uniq_rpkm_ratio'};
            my $b = sprintf "%g", $wanted_CDS{$specific_cds}->{'plus_uniq_rpkm'};
            my $c = sprintf "%.2f", $wanted_CDS{$specific_cds}->{'rpkm_ratio'};
            my $d = sprintf "%g", $wanted_CDS{$specific_cds}->{'plus_rpkm'};
            print $specific_cds, "\t", $a, "\t", $b, "\t", $c, "\t", $d, "\n", ;           
        }
    }
}

sub record_positive_rpkms { 
    my $plus_files_ref   = $_[0];
    my $hash_to_fill_ref = $_[1];
    my $permit_not_all   = $_[2];
    my $number_of_files  = @{ $plus_files_ref };
    my $seen_ref;

    # Read in data from 'plus'/'all' files:
    foreach my $input_file ( @{ $plus_files_ref } ) { 
        open my $INFILE, '<', $input_file
            or die "Can't open input file $input_file: $!";
        while (my $input = <$INFILE>) { 
            if ($input =~ /\A (\S+) \s+ \S+ \s+ (\S+) \s+ (\S+) /xms) { 
                my $cds    = $1;
                my $rpkm   = $2;
                my $multis = $3;
                my $uniqs  = (1.00 - $multis);

                # Mark genes with even one zero-score for oblivion, if 'all' required:
                if ( ($rpkm < 0.01) and ($permit_not_all == 0) ) { 
                    $hash_to_fill_ref->{$cds}->{'delete'} = 1;
                }

                if (! $seen_ref->{$input_file}->{$cds} ) { 
                    $hash_to_fill_ref->{$cds}->{'positive'}  += 1;
                    $hash_to_fill_ref->{$cds}->{'plus_rpkm'} 
                        += ( $rpkm / $number_of_files );
                    $hash_to_fill_ref->{$cds}->{'plus_uniq_rpkm'}
                        += ( ( $rpkm * $uniqs ) / $number_of_files );
                    $seen_ref->{$input_file}->{$cds} = 1;
                }
            }
        }
        close $INFILE or die "Can't close filehandle for $input_file: $!";
    }

    if ( $permit_not_all == 0 ) {
        foreach my $seen_cds ( sort keys %{ $hash_to_fill_ref } ) { 
            if ( $hash_to_fill_ref->{$seen_cds}->{'delete'} ) { 
                # Absolutely, positively remove all other data.
                delete $hash_to_fill_ref->{$seen_cds};
                # Then (re-)mark this unambivalently with a do-not-touch flag.
                $hash_to_fill_ref->{$seen_cds}->{'delete'} = 1;
            }
        }
    }
    return;
}

sub generate_pos_neg_rpkm_ratios {
    my $minus_files_ref   = $_[0];
    my $hash_to_fill_ref  = $_[1];
    my $permit_nonzeros   = $_[2];
    my $number_of_files   = @{ $minus_files_ref };
    my $seen_ref;

    # Read in data from 'minus'/'none' files:
    foreach my $input_file ( @{ $minus_files_ref } ) {
        open my $INFILE, '<', $input_file
            or die "Can't open input file $input_file: $!";
        while (my $input = <$INFILE>) {

            if ($input =~ /\A (\S+) \s+ \S+ \s+ (\S+) \s+ (\S+) /xms) {
                my $cds    = $1;
                my $rpkm   = $2; 
                my $multis = $3;
                my $uniqs  = (1.00 - $multis);

                # Mark genes with even one non-zero score for oblivion, if 'none' required:
                if ( ($rpkm >= 0.01) and ($permit_nonzeros == 0) ) {
                    $hash_to_fill_ref->{$cds}->{'delete'} = 1;
                }

                # Skip any gene banned by the 'all filter earlier.
                if ( (! exists $hash_to_fill_ref->{$cds}->{'delete'} ) 
                     and (! $seen_ref->{$input_file}->{$cds} ) ) {
                    $hash_to_fill_ref->{$cds}->{'minus_rpkm'}
                        += ( $rpkm / $number_of_files );
                    $hash_to_fill_ref->{$cds}->{'minus_uniq_rpkm'}
                        += ( ( $rpkm * $uniqs ) / $number_of_files );
                    $seen_ref->{$input_file}->{$cds} = 1;
                }
            }
        }
        close $INFILE or die "Can't close filehandle for $input_file: $!";
    }

    # Code copied straight from previous subroutine:
    if ( $permit_nonzeros == 0 ) {
        foreach my $seen_cds ( sort keys %{ $hash_to_fill_ref } ) {
            if ( $hash_to_fill_ref->{$seen_cds}->{'delete'} ) {
                # Absolutely, positively remove all other data.
                delete $hash_to_fill_ref->{$seen_cds};
                # Then (re-)mark this unambivalently with a do-not-touch flag.
                $hash_to_fill_ref->{$seen_cds}->{'delete'} = 1;
            }
        }
    }

    # Proceed to rest of calculations:
    foreach my $comp_cds (sort keys %{ $hash_to_fill_ref } ) { 
        # Basic comparison of RPKMs: ratio of avg. positive to avg. negative.
        if ( exists $hash_to_fill_ref->{$comp_cds}->{'plus_rpkm'} ) { 
            my $plus_rpkms  = $hash_to_fill_ref->{$comp_cds}->{'plus_rpkm'};
            my $minus_rpkms = $hash_to_fill_ref->{$comp_cds}->{'minus_rpkm'};
            # To avoid division by zero, enforce nonzero minimum RPKM:
            if ( $minus_rpkms < 0.01) { 
                $minus_rpkms = 0.01;
            }
            $hash_to_fill_ref->{$comp_cds}->{'rpkm_ratio'} 
                = ( $plus_rpkms / $minus_rpkms );
        }

        # Exactly same comparison but with only non-multi reads.
        if ( exists $hash_to_fill_ref->{$comp_cds}->{'plus_uniq_rpkm'} ) {
            my $plus_uniq_rpkms  = $hash_to_fill_ref->{$comp_cds}->{'plus_uniq_rpkm'};
            my $minus_uniq_rpkms = $hash_to_fill_ref->{$comp_cds}->{'minus_uniq_rpkm'};
            # To avoid division by zero, enforce nonzero minimum RPKM:
            if ( $minus_uniq_rpkms < 0.01) {
                $minus_uniq_rpkms = 0.01;
            }
            $hash_to_fill_ref->{$comp_cds}->{'uniq_rpkm_ratio'}
                = ( $plus_uniq_rpkms / $minus_uniq_rpkms );
        }
    }
    return;
}

