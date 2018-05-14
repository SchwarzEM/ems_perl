#!/usr/bin/env perl

# count_fasta_residues.pl -- Erich Schwarz <ems394@cornell.edu>, 2/8/2013. 
# Purpose: get basic statistics from FASTA files, either 'dna' or 'protein'; for DNA, get both scaffolds and true contigs; counts residues directly.

use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;
use List::MoreUtils qw(uniq);

my @infiles         = ();
my $sequence_type   = q{};
my $scaffold_name   = q{};
my $scaffold_abbrev = q{};
my $residue_type    = q{};
my $z               = q{};  # Text spacer to make output lines come out even; '$z' is easier to scan in code for printed lines.

my $basename;
my $concise;
my $extra;
my $redundant;
my $help;

my $nt_count     = 0;
my $N_count      = 0;
my $GC_count     = 0;
my $soft_count   = 0;

my $scaf_name    = q{};
my @scaf_sizes   = ();
my @contig_sizes = ();
my %scaf_seqs    = ();

my $opening_line = "\n";   # First empty line to make reports more readable; for 2+ reports, make it q{} after the first report.

GetOptions ( 'infiles=s{,}' => \@infiles,
             'basename'     => \$basename,
             'type=s'       => \$sequence_type,
             'concise'      => \$concise,
             'extra'        => \$extra,
             'redundant'    => \$redundant,
             'help'         => \$help,   );

if ( $help or (! @infiles) or ( $concise and $extra ) ) { 
    die "Format: count_fasta_residues.pl\n",
        "    --infile|-i     <input stream/files>\n",
        "    --basename|-b   [only give input file basenames]\n",
        "    --type|-t       [sequence type: either 'dna' or 'protein|prot'; default is 'dna']\n",
        "    --concise|-c    [only compute number of sequences and their total length; mutually incompatible with --extra]\n",
        "    --extra|-e      [give mean, median, SDs of contigs/scaffolds]\n",
        "    --redundant|-r  [tolerate redundant sequence names; use with caution and only for messed-up FASTA]\n",
        "    --help|-h       [print this message]\n",
        ;
}

if ($sequence_type) { 
    if ( ( $sequence_type ne 'dna' ) and ( $sequence_type ne 'protein' ) and ( $sequence_type ne 'prot' ) ) { 
        die "Sequence type must be 'dna' or 'protein|prot', not \"$sequence_type\"\n";
    }
    if ( $sequence_type eq 'prot' ) {
        $sequence_type = 'protein';
    }
    if ( $sequence_type eq 'protein' ) { 
        $scaffold_name   = 'Protein';
        $scaffold_abbrev = 'Prot.';
        $residue_type    = 'aa';
        $z               = ' ';
    }
}

if (! $sequence_type) {
    $sequence_type = 'dna';
}

if ( $sequence_type eq 'dna' ) {  
    $scaffold_name   = 'Scaffold';
    $scaffold_abbrev = 'Scaf.';
    $residue_type    = 'nt';
}

foreach my $infile (@infiles) { 
    # For each new file, rezero all of these!  Otherwise, one gets bad carryover errors.
    $nt_count     = 0;
    $N_count      = 0;
    $GC_count     = 0;
    $soft_count   = 0;

    $scaf_name    = q{};
    @scaf_sizes   = ();
    @contig_sizes = (); 
    %scaf_seqs    = (); 

    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
        $infile = basename $infile if $basename;
    }
    while (my $input = <$INPUT_FILE>) { 
        chomp $input;
        if ( $input =~ / \A > (\S+) /xms ) { 
            $scaf_name = $1;
            if ( exists $scaf_seqs{$scaf_name} ) {
                # standard default, correct for almost all circumstances
                if (! $redundant ) {
                    die "Redundant scaffold name: $scaf_name\n";
                }
                # use this if, and only if, it is necessary to count FASTA residues in a truly boneheadly-named file
                else { 
                    my $i = 1;
                    my $new_scaf_name = "$scaf_name.$i";
                    if ( exists $scaf_seqs{$new_scaf_name} ) {
                        while ( exists $scaf_seqs{$new_scaf_name} ) { 
                            $i++;
                            $new_scaf_name = "$scaf_name.$i";
                        }
                    } 
                    warn "Renaming one sequence from redundant name \"$scaf_name\" to non-redundant name \"$new_scaf_name\"\n";
                    $scaf_name = $new_scaf_name;
                    if ( exists $scaf_seqs{$scaf_name} ) {
                        die "This should not have happened, but redundant scaffold name: $scaf_name\n";
                    }
                }
            }
            if ($nt_count) { 
                push @scaf_sizes, $nt_count;
            }
            $nt_count = 0;
        }
        if ( ( $input !~ / \A > /xms) 
            and ($input =~ /\A \S+ \z/xms ) ) { 
            $input =~ s/\s//g;
            $input =~ s/[^a-zA-Z]//g;
            $scaf_seqs{$scaf_name} .= $input;

            $nt_count += length($input);

            if ( $sequence_type eq 'dna' ) {
                $N_count  += ( $input =~ tr/n/n/ );
                $N_count  += ( $input =~ tr/N/N/ );

                $GC_count += ( $input =~ tr/c/c/ );
                $GC_count += ( $input =~ tr/C/C/ );
                $GC_count += ( $input =~ tr/g/g/ );
                $GC_count += ( $input =~ tr/G/G/ );

                $soft_count += ( $input =~ tr/a/a/ );
                $soft_count += ( $input =~ tr/c/c/ );
                $soft_count += ( $input =~ tr/g/g/ );
                $soft_count += ( $input =~ tr/t/t/ );
            }
        }
    }
    # Finish off data at end of file.
    if ($nt_count) {
        push @scaf_sizes, $nt_count;
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";

    if ( ( $sequence_type eq 'dna' ) and (! $concise ) ) {
        foreach my $scaf_name2 (sort keys %scaf_seqs) { 
            $scaf_seqs{$scaf_name2} =~ tr/n/N/;
            # Relentlessly weed out any 'contig' that is not what it is supposed to be -- 1+ nt of pure ACGTacgt.
            my @contig_seqs = grep { $_ =~ /\A [ACGTacgt]+ \z/xms } split /[N]+/, $scaf_seqs{$scaf_name2};
            foreach my $contig_seq (@contig_seqs) { 
                my $contig_len = 0;
                if ( $contig_seq and ( length($contig_seq) >= 1 ) ) {
                    $contig_len = length($contig_seq);
                }
                if ($contig_len) { 
                    push @contig_sizes, $contig_len;
                }
            }
        }
    }

    if (! $concise ) { 
        # Sort in ascending numerical order, so that homebrewed N50 subroutine can work.
        @scaf_sizes   = sort { $a <=> $b } @scaf_sizes;
        @contig_sizes = sort { $a <=> $b } @contig_sizes if ( $sequence_type eq 'dna' );
    }

    my $stat1 = Statistics::Descriptive::Full->new();
    $stat1->add_data(@scaf_sizes);

    my $stat2;
    if ( ( $sequence_type eq 'dna' ) and (! $concise ) ) { 
        $stat2 = Statistics::Descriptive::Full->new(); 
        $stat2->add_data(@contig_sizes);
    }

    # Total nt of sequence in scaffolds, including N residues:
    my $scaf_sum          = $stat1->sum();
    my $raw_scaf_sum      = $scaf_sum;

    my $raw_nonN_scaf_sum = q{};
    my $perc_nonN         = q{};
    my $perc_GC           = q{};
    my $perc_soft1        = q{};
    my $perc_soft2        = q{};

    if ( ( $sequence_type eq 'dna' ) and (! $concise ) ) {
        $raw_nonN_scaf_sum = $raw_scaf_sum - $N_count;
        $perc_nonN         = ($raw_nonN_scaf_sum / $raw_scaf_sum) * 100;
        $perc_GC           = ($GC_count / $raw_nonN_scaf_sum )    * 100;
        $perc_soft1        = ($soft_count / $raw_scaf_sum)        * 100;
        $perc_soft2        = ($soft_count / $raw_nonN_scaf_sum )  * 100;
    }

    $scaf_sum = commify($scaf_sum);

    my $nonN_scaf_sum = q{};
    if ( ( $sequence_type eq 'dna' ) and (! $concise ) ) {
        $nonN_scaf_sum = commify($raw_nonN_scaf_sum);
        $N_count       = commify($N_count);
        $perc_nonN     = sprintf("%.2f", $perc_nonN);
        $perc_GC       = sprintf("%.2f", $perc_GC);
        $perc_soft1    = sprintf("%.2f", $perc_soft1);
        $perc_soft2    = sprintf("%.2f", $perc_soft2);
    }

    # Total nt of sequence in contigs, *not* including N residues joining them in scaffolds:
    my $raw_contig_sum = $stat2->sum() if ( ( $sequence_type eq 'dna' ) and (! $concise ) ); 

    my $scaffolds = $stat1->count();  # Total number of scaffolds.
    $scaffolds    = commify($scaffolds);

    my $scaf_min = q{};
    my $scaf_max = q{};
    my $scaf_n50 = q{};
    my $scaf_n90 = q{};

    if (! $concise ) { 
        $scaf_min     = $stat1->min();
        $scaf_min        = commify($scaf_min);

        $scaf_max     = $stat1->max();
        $scaf_max        = commify($scaf_max);

        $scaf_n50 = get_nXX( \@scaf_sizes, $raw_scaf_sum, 50 );
        $scaf_n50    = sprintf("%.1f", $scaf_n50);
        $scaf_n50    = commify($scaf_n50);

        $scaf_n90 = get_nXX( \@scaf_sizes, $raw_scaf_sum, 90 );
        $scaf_n90    = sprintf("%.1f", $scaf_n90);
        $scaf_n90    = commify($scaf_n90);
    }

    my $contigs    = q{};  # Total number of true contigs.
    my $contig_min = q{};
    my $contig_max = q{};
    my $contig_n50 = q{};
    my $contig_n90 = q{};

    if ( ( $sequence_type eq 'dna' ) and (! $concise ) ) {
        $contigs    = $stat2->count();  # Total number of true contigs.
        $contigs    = commify($contigs);

        $contig_min = $stat2->min();
        $contig_min = commify($contig_min);

        $contig_max = $stat2->max();
        $contig_max = commify($contig_max);

        $contig_n50 = get_nXX( \@contig_sizes, $raw_contig_sum, 50 ); 
        $contig_n50 = sprintf("%.1f", $contig_n50);
        $contig_n50 = commify($contig_n50);
    
        $contig_n90 = get_nXX( \@contig_sizes, $raw_contig_sum, 90 );
        $contig_n90 = sprintf("%.1f", $contig_n90);
        $contig_n90 = commify($contig_n90);
    }

    my $scaf_mean    = q{};
    my $scaf_std_dev = q{};
    my $scaf_median = q{};

    my $contig_mean    = q{};
    my $contig_std_dev = q{};
    my $contig_median = q{};

    if ($extra) { 
        $scaf_mean   = $stat1->mean();
        $scaf_mean   = sprintf("%.1f", $scaf_mean);
        $scaf_mean   = commify($scaf_mean);

        $scaf_std_dev   = $stat1->standard_deviation();
        $scaf_std_dev   = sprintf("%.1f", $scaf_std_dev);
        $scaf_std_dev   = commify($scaf_std_dev);

        $scaf_median   = $stat1->median();
        $scaf_median   = sprintf("%.1f", $scaf_median);
        $scaf_median   = commify($scaf_median);

        if ( $sequence_type eq 'dna' ) {
            $contig_mean = $stat2->mean();
            $contig_mean = sprintf("%.1f", $contig_mean);
            $contig_mean = commify($contig_mean);

            $contig_std_dev = $stat2->standard_deviation();
            $contig_std_dev = sprintf("%.1f", $contig_std_dev);
            $contig_std_dev = commify($contig_std_dev);

            $contig_median = $stat2->median();
            $contig_median = sprintf("%.1f", $contig_median);
            $contig_median = commify($contig_median);
        }
    }

    print $opening_line if $opening_line;  # "\n", but only the first time
    $opening_line = q{};

    print "Sequence: $infile\n";

    print "\n";
    print "Total $residue_type:           $scaf_sum\n";
    print $scaffold_name, "s:$z          $scaffolds\n";
    print "Contigs:            $contigs\n" if ( ( $sequence_type eq 'dna' ) and (! $concise ) );
    print "\n";

    if ( ( $sequence_type eq 'dna' ) and (! $concise ) ) {
        print "ACGT $residue_type:            $nonN_scaf_sum\n";
        print "N-res. $residue_type:          $N_count\n";
        print "% non-N:            $perc_nonN\n";
        print "% GC:               $perc_GC\n";
        print "% softmasked [1]:   $perc_soft1   ([1]: versus entire assembly size)\n";
        print "% softmasked [2]:   $perc_soft2   ([2]: versus non-N, ACGT only, assembly size)\n";
        print "\n";
    }

    if (! $concise ) {
        print "$scaffold_name N50 $residue_type:$z    $scaf_n50\n";
        print "$scaffold_name N90 $residue_type:$z    $scaf_n90\n";
        print "\n" if ( $sequence_type eq 'protein' );
        print "$scaffold_abbrev max. $residue_type:      $scaf_max\n";
        print "$scaffold_abbrev min. $residue_type:      $scaf_min\n";
        print "\n";
    }

    if ($extra) { 
        print "$scaffold_abbrev mean $residue_type:      $scaf_mean\n";
        print "$scaffold_abbrev median $residue_type:    $scaf_median\n";
        print "$scaffold_abbrev s. dev. $residue_type:   $scaf_std_dev\n";
        print "\n";
    }

    if ( ( $sequence_type eq 'dna' ) and (! $concise ) ) { 
        print "Contig N50 nt:      $contig_n50\n";
        print "Contig N90 nt:      $contig_n90\n";
        print "Contig max. nt:     $contig_max\n";
        print "Contig min. nt:     $contig_min\n";
        print "\n";

        if ($extra) {
            print "Contig mean nt:     $contig_mean\n";
            print "Cont. median nt:    $contig_median\n";
            print "Cont. s. dev. nt:   $contig_std_dev\n";
            print "\n";
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

# Assumes array sorted in ascending values, with precomputed total value.
# Intended to enable getting N50, N90, etc. as needed.
sub get_nXX {
    my ($sizes_ref, $sum, $nval) = @_;
    $nval = 100 - $nval;
    if ( ( $nval <= 0 ) or ( $nval >= 100 ) ) { 
        die "Subroutine get_nXX only accepts arguments for Nvalue greater than 0 and less than 100\n";
    }
    my $midsum = ($sum * ( $nval / 100 ) );
    my $tally = 0;
    my @_sizes = @{ $sizes_ref };
    foreach my $value (@_sizes) {
        $tally += $value;
        if ($tally > $midsum) {
            return $value;
        }
    }
    return;      
}                

