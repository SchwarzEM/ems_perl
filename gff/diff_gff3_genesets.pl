#!/usr/bin/env perl

# diff_gff3_genesets.pl -- Erich Schwarz <ems394@cornell.edu>, 5/29/2018.
# Purpose: list non-/overlap genes in GFF3s; optional diskfile, prog. rep.; requires 2018 format (GFF3-like must be reformatted).

use strict;
use warnings;
use autodie;

use Getopt::Long;
use DBM::Deep;
use Tie::File;
use File::Basename;

my $tiefile;
my $prog_rept;

# We may have 2+ sets of genes with identical names.
# To keep this from being a problem, use internal names ("gene_set.$set.$gene"), 
# but print their external names ("$gene").
my $set = 0;

# Compare coordinates of genes listed in 2 (or more) input files.
my @files = ();

my $data_ref;

my $help;

GetOptions ( 'infiles=s{,}' => \@files,
             'tiefile'      => \$tiefile, 
             'progress'     => \$prog_rept,
             'help'         => \$help,      ); 

if ( $help or (! @files ) ) {
     die "Format: diff_genesets.pl\n",
         "    --infile|-i     [input files: compare coordinates of genes listed in 2 (or more) input GFF3 files]\n",
         "    --tiefile|-t    [optionally, to save on RAM use, tie a big hashref to a uniquely named DB file]\n",
         "    --progress|-p   [optionally, track progress as the (potentially long) program runs]\n",
         "    --help|-h       [print this message]\n",
         ;
}

my %tx2gene = ();

# Link transcripts/CDSes to genes.
#
# Instance of orthodox GFF3, should be parsable:
# 
# chrI_pilon	AUGUSTUS	transcript	6054	11981	0.98	-	.	ID=chrI_pilon.g1.t1;Parent=chrI_pilon.g1
# 
# Instance of *un*northodox GFF3, will *not* be parsable unless reformatted first:
# 
# chrIV_pilon	WormBase	transcript	6198	20429	.	+	.	gene_id "WBGene00021406"; transcript_id "Y38C1AB.4"; \
# gene_source "WormBase"; gene_biotype "protein_coding"; transcript_source "WormBase"; transcript_biotype "protein_coding";

my $text_pattern_1
    =   ' \A \S+ '
      . ' \t [^\t]* \t transcript \t '
      . ' \d+ \t \d+ '
      . ' \t [^\t]* \t '
      . ' [+|-] \t [^\t]* '
      . ' \t ID=([^;\s]+);Parent=([^;\s]+) '  # transcript; gene 
      ;

# Second, parse lines that link transcript/CDSes to exonic coordinates:
# 
# Instance of orthodox GFF3, should be parsable:
#
# chrI_pilon	AUGUSTUS	CDS	6054	6191	1	-	0	ID=chrI_pilon.g1.t1.cds;Parent=chrI_pilon.g1.t1
#
# Instance of *un*northodox GFF3, will *not* be parsable unless reformatted first:
# 
# chrIV_pilon	WormBase	CDS	6198	6568	.	+	0	gene_id "WBGene00021406"; transcript_id "Y38C1AB.4"; \
# exon_number "1"; gene_source "WormBase"; gene_biotype "protein_coding"; transcript_source "WormBase"; transcript_biotype "protein_coding"; protein_id "Y38C1AB.4";

my $text_pattern_2
    =   ' \A (\S+) '                        # chr./contig
      . ' \t [^\t]* \t CDS \t '
      . ' (\d+) \t (\d+) '                  # start/end nt
      . ' \t [^\t]* \t ' 
      . ' ([+|-]) \t (\d+) '                # strand ori.; phase
      . ' \t ID=[^;\s]+;Parent=([^;\s]+) '  # tx. name
      ;

# By default, data stored in hashrefs in RAM.
my $gene_exists_ref;
my $gene2equivs_ref;
my $coord2gene_ref;
my $date = q{};

# Optionally, tie a big hashref to a uniquely named DB file:
if ($tiefile) { 
    $date = join('.', &get_local_date());
    my $db_name1       = 'coord2gene.'  . $date . '.db';
    $db_name1          = failsafe_name($db_name1);
    $coord2gene_ref = DBM::Deep->new( $db_name1 );
}

# Optionally, track progress as the (potentially long) program runs.
my $format_string1 = 'Reading: gene %s of %s genes in file %s';
my $format_string2 = 'Reading: %s:%s-%s from gene %s in file %s';
my @progress       = ( q{}, q{}, q{} );
my $prog_rep       = 'diff_genesets.' . $date . '.PROGRESS.txt';
$prog_rep          = failsafe_name($prog_rep);
if ( $prog_rept ) { 
    tie(@progress, "Tie::File", $prog_rep) or die "Can't manipulate progress report $prog_rep: $!";
}

# Read data from (usually two) GFF3 files.
foreach my $file (@files) { 
    # For internal gene names ("gene_set.$set.$gene") versus external ones ("$gene").
    $set++;

    my $no_scanned = 0;

    # Do first file-reading just to link transcripts to genes:
    open my $FILE, '<', $file;
    while (my $input = <$FILE>) { 
        chomp $input;

        # Link transcripts or CDSes to *unique* genes:
        if ( $input =~ / $text_pattern_1 /xms ) { 
            my $tx       = $1;
            my $gene     = $2;
            my $int_gene = "gene_set.$set.$gene";

            # We care about uniqueness of transcript-to-gene mapping, but only *within* one gene/transcript set.
            if ( exists $data_ref->{'tx'}->{$tx}->{'input_file'}->{$file}->{'int_gene'} ) { 
                die "Transcript/CDS $tx already mapped in $file to (internal) gene",
                    " $data_ref->{'tx'}->{$tx}->{'input_file'}->{$file}->{'int_gene'}",
                    " cannot also map to (internal) gene $int_gene!\n",
                    ;
            }
            $data_ref->{'tx'}->{$tx}->{'int_gene'} = $int_gene;
            $data_ref->{'int_gene'}->{$int_gene}->{'gene'} = $gene;
            $data_ref->{'tx'}->{$tx}->{'input_file'}->{$file}->{'int_gene'} = $int_gene;

            if ( $prog_rept ) { 
                $no_scanned++;
                my $commas_no_scnd = commify($no_scanned);
                $progress[1] = sprintf( "$format_string1\n", 
                                         $gene, 
                                         $commas_no_scnd, 
                                         $file, );
            }
        }
    }
    close $FILE;

    # Do second file-reading to link genes (via transcripts/CDSes) to exons:
    open $FILE, '<', $file;
    while (my $input = <$FILE>) {
        chomp $input;
        if ( $input =~ / $text_pattern_2 /xms ) {
            my $chr = $1;
            my $nt1 = $2;
            my $nt2 = $3;
            my $ori = $4;
            my $phs = $5;
            my $tx  = $6;

            my $int_gene = q{};
            my $gene     = q{};

            # Report errors, but don't kill whole analysis over one error:
            if (! exists $data_ref->{'tx'}->{$tx}->{'int_gene'} ) { 
                warn "Failed to link transcript $tx to gene, so cannot parse: $input\n";
            }
            else {
                $int_gene = $data_ref->{'tx'}->{$tx}->{'int_gene'};
                if (! exists $data_ref->{'int_gene'}->{$int_gene}->{'gene'} ) {
                    warn "Failed to link internal gene name $int_gene to external gene name, so cannot parse: $input\n";
                }
                else {
                    $gene = $data_ref->{'int_gene'}->{$int_gene}->{'gene'};
                }
            }
            # But do *not* list $gene in $gene_exists_ref until it gets defined coordinates.

            if ( $prog_rept ) { 
                 $progress[1] = sprintf( "$format_string2\n", $chr, $nt1, $nt2, $gene, $file, );
            }

            # For general work, to save wasted effort, do not automatically 
            #     cross-check every nucleotide for overlapping genes:
            my $check_indiv_nt = 0;

            # N.B., data being checked is operationally defined by this step:
            # 
            # $coord2gene_ref->{$chr}->{$nt}->{$ori}->{$phs}->{$int_gene} = 1;
            # 
            # Imprecise queries will end up trying to match HASH(0x8631928) to '+' or '-'!

            # Do laborious cross-checking only if several conditions exist:
            if (      ( exists $coord2gene_ref->{$chr}           ) 
                  and ( @{ $coord2gene_ref->{$chr} }{$nt1..$nt2} ) ) {  
                my @prior_oris = ();
                foreach my $nt ($nt1..$nt2) { 
                    my @new_oris = grep { $_ eq $ori } 
                                   keys %{ $coord2gene_ref->{$chr}->{$nt} };
                    push @prior_oris, @new_oris;
                }
                if (@prior_oris) { 
                    $check_indiv_nt = 1;
                }
            }

            # Having made the decision whether to cross-check coords.:
            foreach my $nt ($nt1..$nt2) { 
                my @overlaps = ();

                # Optionally, check for overlaps nucleotide by nucleotide:
                if (     ( $check_indiv_nt                               ) 
                     and ( exists $coord2gene_ref->{$chr}->{$nt}->{$ori}->{$phs} ) ) { 
                    # Keys are gene names; values are all '1'.
                    # Do not let multiple CDS/txs make a gene 'match' itself!
                    @overlaps = grep { $_ ne $int_gene } 
                                keys %{ $coord2gene_ref->{$chr}->{$nt}->{$ori}->{$phs} };
                    if (@overlaps) { 
                        foreach my $overlapping (@overlaps) { 
                            # Record link in both directions!
                            $gene2equivs_ref->{$int_gene}->{$overlapping} = 1;
                            $gene2equivs_ref->{$overlapping}->{$int_gene} = 1;
                        }
                    }
                }

                # And then archive the gene's own data.
                $coord2gene_ref->{$chr}->{$nt}->{$ori}->{$phs}->{$int_gene} = 1;

                # Once we know that a gene has defined coordinates,
                #    *then* we list it for later recollection.
                #    We are not interested in genes w/o coordinates!
                # Also, it is useful to know what genes were from what files.
                $gene_exists_ref->{$int_gene} = $file;
            }
        }
    }
    close $FILE;
}

# At long last!!  A nice, simple printout:
foreach my $int_gene ( sort keys %{ $gene_exists_ref } ) { 
    my $gene = q{}; 
    $gene = $data_ref->{'int_gene'}->{$int_gene}->{'gene'};
    print "$gene";

    # I want the *basenames* of source files, not full-length names.
    my $file_source = $gene_exists_ref->{$int_gene};
    $file_source = basename($file_source);
    print "\t$file_source";
    $file_source = q{};   # zero out previous values to prevent carry-over

    foreach my $overlap (sort keys %{ $gene2equivs_ref->{$int_gene} } ) { 
        $gene = $data_ref->{'int_gene'}->{$overlap}->{'gene'};
        print "\t$gene";

        # And, get basename of source file for $overlap:
        $file_source = $gene_exists_ref->{$overlap};
        $file_source = basename($file_source);
        print "\t$file_source";
        $file_source = q{};   # zero out
    }
    print "\n";
}


### Subroutines. ###

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

sub get_local_date {     
    my @ltime = localtime;
    my @ldate = ( (sprintf ("%04u", ($ltime[5] + 1900)) ),     # $year
                  (sprintf ("%02u", ($ltime[4] + 1))    ),     # $mon
                  (sprintf ("%02u", ($ltime[3] + 0))    ),     # $mday
                  (sprintf ("%02u", ($ltime[2] + 0))    ),     # $hour
                  (sprintf ("%02u", ($ltime[1] + 0))    ),     # $min
                  (sprintf ("%02u", ($ltime[0] + 0))    ), );  # $sec
    return @ldate;   
}

sub commify { 
    my $text = reverse $_[0];
    $text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $text;
}

