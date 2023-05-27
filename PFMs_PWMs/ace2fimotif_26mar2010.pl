#!/usr/bin/env perl

# ace2fimotif.pl -- Erich Schwarz <ems394@cornell.edu>, 3/26/2010 (legacy version).
# Purpose: generate a minimal MEME motif usable by FIMO.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);  # Perl Cookbook 2.1.

my $name         = q{};  # Optionally for *one* motif; can't list 2+; otherwise, get names from .ace.
my $cg_freq      = q{};  # Optionally user-provided; otherwise defaults to C. elegans value.
my $sites_used   = q{};  # Optionally user-provided (defaults to 1 -- to hand-edit! -- or no. in .ace).
my $width        = q{};
my $files_output = q{};

my @res_list   = qw( A C G T );
my @raw_values = ();
my @sum_values = ();
my $max_val_ref;
my %scanned    = ();

my $fimotifs_ref;
my $help;

my $record_data  = 0;
my $residue      = q{};
my $value_string = q{};

my $matrix_type = 'letter-probability matrix';  # log-odds matrix' not allowed, AFAIK.

# Use this variable to allow rounding with user-specified precision:
my $DIGITS = 'no';

GetOptions ( 'name:s'          => \$name,
             'cg_freq:s'       => \$cg_freq,
             'sites_used:s'    => \$sites_used,
             'rounding=i{0,1}' => \$DIGITS,
             'files_output'    => \$files_output,
             'help'            => \$help,   );

if ( ($help) or (! @ARGV) ) { 
    die "Format: ace2fimotif.pl",
        " --name|-n [if not in input stream/file]",
        " --cg_freq|-c [0.00-1.00]",
        " --sites_used|-s [REQ.: 1+ sites used]",
        " --rounding|-r [digits rounding; default 4]",
        " --files_output|-f [output 1 file/motif]",
        " <input stream/files>\n",
        ;
}

if ( $cg_freq and ( ( $cg_freq < 0 ) or ( $cg_freq > 1 ) ) ) { 
    die "CG frequency must be in 0.00-1.00 range, not be: $cg_freq!\n";
}

$cg_freq  ||= 0.354;     # WS210: precise value is ( 35_537_706 / 100_272_210 );
my $at_freq = 1 - $cg_freq;
my $a_freq  = $at_freq/2;
my $t_freq  = $a_freq;
my $c_freq  = $cg_freq/2;
my $g_freq  = $c_freq;

$sites_used ||= 1;

# If user just says '-r ', then default to four digits.
$DIGITS or ($DIGITS = 4);
# If it's still 'no', *then* set to q{}:
if ($DIGITS eq 'no') { 
    $DIGITS = q{};
}

while (my $input = <>) { 
    chomp $input;

    # Default is getting names from .ace; this allows >= 2 names/matrices.
    if ( $input =~ /\A Position_Matrix \s+ : \s+ \" ( [^\"\s]+ ) \" /xms ) { 
        $name         = $1;
        zero_several_values();
        # Default, but can be overridden by data in .ace:
        $fimotifs_ref->{$name}->{'sites_used'} = $sites_used;
    }

    # Specifically require that this be seen for data to be recorded:
    if ( $input =~ /\A Type \s+ Frequency \b /xms ) { 
        $record_data = 1;
    }

    # Since fimo doesn't seem to allow PWM inputs, don't use them.
    # Stop reading, zero out recently read data, delete the motif's data by name, and *then* delete its name.
    # Being any less thorough leads to unhappy error messages later, when the program tries to print half-entered motifs.
    if ( $input =~ /\A Type \s+ Weight \b /xms ) {
        zero_several_values();
        delete ${ $fimotifs_ref }{$name};
        $name = q{};
    }

    # If the .ace actually gives a site count, use that and not the default:
    if ( $record_data and ( $input =~ /\A Sites_used \s+ (\d+) \b /xms ) ) { 
        $fimotifs_ref->{$name}->{'sites_used'} = $1;
    }

    # Read residue values from the PFM:
    if ( $record_data and ( $input =~ /\A Site_values \s+ ([ACGT]) \s+ ( \S+ (?: \s+ \S+)* ) /xms ) ) { 
        $residue      = $1;
        $value_string = $2;

        # Enforce a name:
        if (! $name) { 
            die "What is the name of the matrix?\n";
        }

        # Remove comment, and then any trailing whitespace:
        $value_string =~ s/[#].*\z//;
        $value_string =~ s/\s+\z//;

        # Get a list of raw values; enforce their number-ness.
        @raw_values = split /\s+/, $value_string;
        foreach my $value (@raw_values) { 
            if (! looks_like_number($value) ) { 
                die "Non-numerical value ($value) in input data: $input\n";
            }            
        }

        # Get the width and enforce uniform widths.
        $width = @raw_values;
        if ( ( exists $fimotifs_ref->{$name}->{'width'} )
              and ( $width != $fimotifs_ref->{$name}->{'width'} ) ) { 
            die "Matrix $name has inconsistent widths: $fimotifs_ref->{$name}->{'width'} versus $width!\n";
        }
        if (! exists $fimotifs_ref->{$name}->{'width'} ) { 
            $fimotifs_ref->{$name}->{'width'} = $width;
        }

        # Store the raw values:
        @{ $fimotifs_ref->{$name}->{$residue}->{'raw'} } = @raw_values;

        # Mark off yet another residue's values as scanned:
        $scanned{$residue} = 1;

        # If we've worked through A, C, G, and T, then normalize the data, and polish off to 1.0 sum:
        if ( $scanned{'T'} and $scanned{'G'} and $scanned{'C'} and $scanned{'A'} ) {

            # Remember: Perl counts from 0 to $i-1!
            foreach my $i (0 .. ($width-1)) { 
                $sum_values[$i] = 0;
                foreach my $residue2 (@res_list) { 
                    $sum_values[$i] += $fimotifs_ref->{$name}->{$residue2}->{'raw'}->[$i];
                }
                if ( $sum_values[$i] == 0 ) { 
                    warn "Can't process: $input\n";
                    die "Can't normalize values to zero-sum total!\n";
                } 

                # Get ready to polish to 1.0 sum:
                $max_val_ref->{'maxval'} = 0;

                foreach my $residue3 (@res_list) { 
                    # This works, but tends to give absurdly precise decimals:
                    my $normalized 
                        = ( $fimotifs_ref->{$name}->{$residue3}->{'raw'}->[$i] / $sum_values[$i] );

                    # So (as a conscious choice, not a default) round them off:
                    if ($DIGITS) { 
                        ( $normalized = sprintf("%.${DIGITS}f", $normalized) ) 
                            or die "Can't round off $normalized!\n";
                    }

                    # Whether rounded or not, store the data:
                    $fimotifs_ref->{$name}->{$residue3}->{'norm'}->[$i] = $normalized;

                    # Keep track of which value was largest:
                    if ( $normalized > $max_val_ref->{'maxval'} ) { 
                        $max_val_ref->{'maxval'} = $normalized;
                        $max_val_ref->{'maxres'} = $residue3;
                    }
                }
                my $sum_non_maxvals = 0;
                foreach my $residue4 ( grep { $_ ne $max_val_ref->{'maxres'} } @res_list) { 
                    $sum_non_maxvals += $fimotifs_ref->{$name}->{$residue4}->{'norm'}->[$i];
                }
                my $polished_maxval = 1 - $sum_non_maxvals;
                my $maxres = $max_val_ref->{'maxres'};

                # Again, round off to amount needed/wanted:
                if ($DIGITS) {
                    ( $polished_maxval = sprintf("%.${DIGITS}f", $polished_maxval) )
                            or die "Can't round off $polished_maxval!\n";
                }
                $fimotifs_ref->{$name}->{$maxres}->{'norm'}->[$i] = $polished_maxval;
            }
        }
    } 
}

# Output a minimal motif file per MEME:

if (! $files_output ) { 
    print_fimotif_header();
    print_indiv_fimotifs($fimotifs_ref);
}

if ($files_output) { 
    print_indiv_fimotif_files($fimotifs_ref);
}

sub zero_several_values {
    $record_data  = 0;  
    $value_string = q{};
    @raw_values   = ();
    @sum_values   = ();
    undef $max_val_ref;
    %scanned      = ();
}

sub print_fimotif_header { 
    # The documentation says '3.0'; I assume accuracy makes more sense.
    print "MEME version 4.3.0\n",
          "\n",
          "ALPHABET= ", @res_list, "\n",
          "strands: + -\n",
          "\n",
          "Background letter frequencies (from .ace file via ace2fimotif.pl):\n",
          "A $a_freq C $c_freq G $g_freq T $t_freq\n",
          "\n",
          ;
}

sub print_indiv_fimotifs { 
    my $_fimotifs_ref = $_[0];
    foreach my $_fimotif (sort keys %{ $_fimotifs_ref } ) { 
        print "MOTIF $_fimotif\n";
        my $_width = $fimotifs_ref->{$_fimotif}->{'width'};
        my $_sites_used = $_fimotifs_ref->{$_fimotif}->{'sites_used'};
        # Does this next line need >0 numbers?
        print "BL   MOTIF $_fimotif width=0 seqs=0\n";
        print "$matrix_type: alength= 4 w= $_width nsites= $_sites_used E= 0\n";
        # $_width, not $width:
        foreach my $i (0 .. ($_width-1)) { 
            my @out_vals = ();
            foreach my $_residue (@res_list) { 
                push @out_vals, $_fimotifs_ref->{$_fimotif}->{$_residue}->{'norm'}->[$i];
            }
            my $out_vals_line = join "\t", @out_vals;
            print "$out_vals_line\n";
        }
        print "\n";
    }
}

sub print_indiv_fimotif_files { 
    my $_fimotifs_ref = $_[0];
    foreach my $_fimotif (sort keys %{ $_fimotifs_ref } ) { 
        my $_outfile = $_fimotif;
        $_outfile = failsafe_name($_outfile);
        open my $_OUTFILE, '>', $_outfile 
            or die "Can't open output file $_outfile: $!";
        print {$_OUTFILE} "MEME version 3.0\n",
                          "\n",
                          "ALPHABET= ", @res_list, "\n",
                          "strands: + -\n",
                          "\n",
                          "Background letter frequencies (from .ace file via ace2fimotif.pl):\n",
                          "A $a_freq C $c_freq G $g_freq T $t_freq\n",
                          "\n",
                          ;
        print {$_OUTFILE} "MOTIF $_fimotif\n";
        my $_width = $fimotifs_ref->{$_fimotif}->{'width'};
        my $_sites_used = $_fimotifs_ref->{$_fimotif}->{'sites_used'};
        # Does this next line need >0 numbers?
        print {$_OUTFILE} "BL   MOTIF $_fimotif width=0 seqs=0\n";
        print {$_OUTFILE} "$matrix_type: alength= 4 w= $_width nsites= $_sites_used E= 0\n";
        # $_width, not $width:
        foreach my $i (0 .. ($_width-1)) {
            my @out_vals = ();
            foreach my $_residue (@res_list) {
                push @out_vals, $_fimotifs_ref->{$_fimotif}->{$_residue}->{'norm'}->[$i]; 
            }
            my $out_vals_line = join "\t", @out_vals;
            print {$_OUTFILE} "$out_vals_line\n";
        }
        print {$_OUTFILE} "\n";
        close $_OUTFILE or die "Can't close output file $_outfile: $!";
    }
}

sub failsafe_name {
    my $filename = $_[0];
    if (-e $filename) {
        my $suffix = 0;
        while (-e $filename) {
            $suffix++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix";
        }
    }
    return $filename;
}

