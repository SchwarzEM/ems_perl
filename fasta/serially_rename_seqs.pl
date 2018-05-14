#!/usr/bin/env perl

# serially_rename_seqs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/26/2015.
# Purpose: rename scaffolds with zero-padded serial numbers; optionally, keep original names but add nt counts to header data.

use strict;
use warnings;
use Getopt::Long;

my $i            = 0;
my $prefix       = q{};
my $infile       = q{};
my $header       = q{};
my $orig_seqname = q{};
my $residues     = q{};

my $data_ref;
my $dna;
my $suppress;
my $keep;
my $badnames;
my $help;

GetOptions ( 'infile=s'   => \$infile,
             'prefix:s'   => \$prefix,
             'residues:s' => \$residues,
             'dna'        => \$dna,
             'suppress'   => \$suppress,
             'keep'       => \$keep, 
             'badnames'   => \$badnames,
             'help'       => \$help,     );

# Have human-readable default value:
$prefix   ||= 'Scaffold_';
# Default to counting nt, not aa:
$residues ||= 'nt';

if ( $help or (! $infile ) ) { 
    die "Format: serially_rename_seqs.pl\n",
        "    --infile|-i    <input stream/files>\n",
        "    --prefix|-p    [optional prefix; default is \"Scaffold_\"]\n",
        "    --residues|-r  [optional description of residues being counted; default is \"nt\"]\n",
        "    --dna|-d       [only accept sequence letters ACGTN/acgtn]\n",
        "    --suppress|-s  [do not save the previous text of the header in the post-name section]\n",
        "    --keep|-k      [keep original scaffold names, rather than using new serialized names]\n",
        "    --badnames|-b  [tolerate bad, redundant sequence names; use with caution and only for messed-up FASTA]\n",
        "    --help|-h\n",
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FILE;
if ($infile eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
}

while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    if ( $input =~ /\A > ((\S+) .*) \z/xms ) { 
        $header       = $1;
        $orig_seqname = $2;
        if ( exists $data_ref->{'orig_seqname'}->{$orig_seqname} ) { 
            # standard default, correct for almost all circumstances
            if (! $badnames ) {
                die "Redundant sequence name: $orig_seqname\n";
            }

            # use this if, and only if, it is necessary to deal with a truly boneheadly-named FASTA file
            else { 
                my $i = 1;
                my $renamed_seqname = "$orig_seqname.renamed.$i";
                if ( exists $data_ref->{'orig_seqname'}->{$renamed_seqname} ) { 
                    while ( exists $data_ref->{'orig_seqname'}->{$renamed_seqname} ) { 
                        $i++;
                        $renamed_seqname = "$orig_seqname.renamed.$i";
                    }
                }
                warn "Renaming one sequence from redundant \"$orig_seqname\" to nonredundant \"$renamed_seqname\"\n";
                $orig_seqname = $renamed_seqname;
                if ( exists $data_ref->{'orig_seqname'}->{$orig_seqname} ) {
                        die "This should not have happened, but redundant sequence name: $orig_seqname\n";
                }
            }
        }
        $data_ref->{'orig_seqname'}->{$orig_seqname}->{'header'} = $header;
    }
    elsif ( $input =~ / \A \s* [A-Za-z] /xms ) { 
        $input =~ s/\s//g;
        if ($dna and ( $input =~ / [^ACGTNacgtn] /xms ) ) { 
            die "Can't parse: $input\n";
        }
        $data_ref->{'orig_seqname'}->{$orig_seqname}->{'sequence'} .= $input;
    }
    else { 
        if ( $input !~ /\A \s* \z/xms ) { 
            die "Can't parse: $input\n";
        }
    }
}

foreach my $orig_seq1 ( keys %{ $data_ref->{'orig_seqname'} } ) { 
    $data_ref->{'orig_seqname'}->{$orig_seq1}->{'length'} 
        = length( $data_ref->{'orig_seqname'}->{$orig_seq1}->{'sequence'} );
}

my $seq_count = keys %{ $data_ref->{'orig_seqname'} };
my $DIGITS = length($seq_count);

# $b <=> $a gives *descending* scaffold sizes, which is what we want:

my @size_sorted_orig_seqnames = sort {     $data_ref->{'orig_seqname'}->{$b}->{'length'} 
                                       <=> $data_ref->{'orig_seqname'}->{$a}->{'length'} } 
                                       keys %{ $data_ref->{'orig_seqname'} };

foreach my $orig_seq2 (@size_sorted_orig_seqnames) { 
    $i++;
    my $serial_no = $i;
    my $sf_format = '%0' . $DIGITS . 'u';
    $serial_no = sprintf($sf_format, $serial_no) or die "Can't zero-pad serial number $serial_no\n";

    # Either serialized names or, optionally, original names.
    my $new_name = q{};
    if ($keep) {
        $new_name = $orig_seq2;
    }
    else {
        $new_name = $prefix . $serial_no;
    }    

    my $humread_len = commify($data_ref->{'orig_seqname'}->{$orig_seq2}->{'length'});
    print '>',
          $new_name,
          q{  },
          "$humread_len ",
          "$residues",
          ;
    if (! $suppress) { 
        print q{  }, $data_ref->{'orig_seqname'}->{$orig_seq2}->{'header'};
    }
    print "\n";

    my @output_lines 
        = unpack("a60" 
                 x ($data_ref->{'orig_seqname'}->{$orig_seq2}->{'length'}/60 + 1), 
                 $data_ref->{'orig_seqname'}->{$orig_seq2}->{'sequence'});
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
    }
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

