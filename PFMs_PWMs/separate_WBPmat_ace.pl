#!/usr/bin/env perl

# separate_WBPmat_ace.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/8/2009.
# Purpose: split .ace matrices into: PFMs; PFM-derived PWMs; free PWMs; anomalies.

use strict;
use warnings;

my $infile     = 'infile';
my $type       = q{};
my $object     = q{};
my @temp_lines = ();

my $object2info_ref;

while (my $input = <>) { 
    if ($ARGV) { 
        $infile = $ARGV;
    }
    if ( $input =~ /\A \s* \z/xms ) { 
        @{ $object2info_ref->{$type}->{$object} } = @temp_lines;
        $type       = q{};
        $object     = q{};
        @temp_lines = ();
    }

    # Order of next three 'if's matters!

    if ( ($type eq 'Weight') and ($input =~ / \A Derived_from_matrix \s /xms) ) { 
        $type = 'Weight.Derived';
    }
    if ( ($object) and ( $input =~ /\S/xms ) ) {
        if ( $input =~ /\A Type \s+ (\S+) /xms ) {
            $type = $1;
            if ( ($type ne 'Frequency') and ($type ne 'Weight') ) { 
                $type = 'Anomaly';
            }
        }
        push @temp_lines, $input;
    }
    if ( (! $object) and ( $input =~ /\S/xms ) ) { 
        if ( $input =~ /\A Position_Matrix \s+ : \s+ \" (WBPmat\d+) \"  /xms ) { 
            $object = $1;
            push @temp_lines, $input;
        }
    }
}

# Clear stored data at EOF:
@{ $object2info_ref->{$type}->{$object} } = @temp_lines;
$type       = q{};
$object     = q{};
@temp_lines = ();

# Output everything to sorted files:
foreach my $matrix_type ( grep { /\S/ } sort keys %{ $object2info_ref } ) { 
    open my $PFMS, '>', "$infile.$matrix_type" 
        or die "Can't open output file for $matrix_type matrices, $infile.$matrix_type: $!";
    foreach my $matrix ( sort keys %{ $object2info_ref->{$matrix_type} } ) { 
        print $PFMS "\n";
        foreach my $line ( @{ $object2info_ref->{$matrix_type}->{$matrix} } ) { 
            print $PFMS $line;
        }
        print $PFMS "\n";
    }
    close $PFMS or die "Can't close filehandle for $infile.$matrix_type: $!";
}

