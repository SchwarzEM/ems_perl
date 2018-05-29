#!/usr/bin/env perl

# diff_genesets.pl -- Erich Schwarz <ems394@cornell.edu>, 5/29/2018.
# Purpose: list non-/overlap genes in GFFs; optional diskfile, prog. rep.

use strict;
use warnings;

use Getopt::Long;
use DBM::Deep;
use Tie::File;

my $tiefile;
my $prog_rept;

# Compare coordinates of genes listed in 2 (or more) input files.
my @files = ();

my $help;

GetOptions ( 'infiles=s{,}' => \@files,
             'tiefile'      => \$tiefile, 
             'progress'     => \$prog_rept,
             'help'         => \$help,      ); 

if ( $help or (! @files ) ) {
     die "Format: diff_genesets.pl\n",
         "    --infile|-i     [input files: compare coordinates of genes listed in 2 (or more) input files]\n",
         "    --tiefile|-t    [optionally, to save on RAM use, tie a big hashref to a uniquely named DB file]\n",
         "    --progress|-p   [optionally, track progress as the (potentially long) program runs]\n",
         "    --help|-h       [print this message]\n",
         ;
}

# Link transcripts/CDSes to genes;
my %tx2gene = ();

# GFF2/3 use different syntaxes to encode data.  We will be reading 
#     willy-nilly from either sort of file, and need to 
#     flexibly pattern-match either syntax as needed.

# First, parse lines that link transcript/CDSes to genes:
my $text_pattern_1  = q{}; 

# Sample inputs, GFF2 vs. GFF3:
# 
# GFF2: 
# II	curated	CDS	9198260	9206852	.	+	.	CDS "ZK1067.1" ; 
#     Note "tyrosine-protein kinase (Epidermal growth factor receptor subfamily)" ; 
#     WormPep "WP:CE03840" ; Note "let-23" ;  Status "Confirmed" ;  Gene "WBGene00002299" ; 
# 
# GFF3:
# 
# chrI	mGENE	mRNA	24696	27804	.	+	.	
#     ID=Transcript:mGene_pred_1.1;Parent=Gene:mGene_pred_1;prediction_status=predicted

my $text_pattern_1_GFF2 
    = '\s CDS \s+ \" ([^\"]+) \" .+ \s+ Gene \s+ \" (WBGene\d+) \"';

# The only syntax that seems to work is using nongreedy '(\S+?)', not (?:^(\s;)+)':
my $text_pattern_1_GFF3 = 'ID=Transcript:(\S+?);Parent=Gene:(\S+?);';

# Second, parse lines that link transcript/CDSes to exonic coordinates:
my $text_pattern_2  = q{};

# Sample inputs, GFF2 vs. GFF3:
#
# GFF2:    
# II	Coding_transcript	coding_exon	9198260	9198374	.	+	0	
#     Transcript "ZK1067.1" ; CDS "ZK1067.1"
# 
# GFF3:
# chrI	mGENE	CDS	24703	24787	.	+	0	
#     ID=CDS:mGene_pred_1.1;Parent=Transcript:mGene_pred_1.1

my $text_pattern_2_GFF2 
    =   ' \A (\S+) '                   # chr./contig
      . ' \t \S+ \t coding_exon \t '
      . ' (\d+) \t (\d+) '             # start/end nt
      . ' \t \. \t '
      . ' ([+|-]) \t (\d+) '           # strand ori.; phase
      . ' \t .+ \s+ CDS \s+ '
      . ' \" ([^\"]+) \" '             # CDS name
      ;

my $text_pattern_2_GFF3 
    =   ' \A (\S+) '                         # chr./contig
      . ' \t \S+ \t CDS \t '
      . ' (\d+) \t (\d+) '                   # start/end nt
      . ' \t \. \t ' 
      . ' ([+|-]) \t (\d+) '                 # strand ori.; phase
      . ' \t ID=\S+Parent=Transcript:(\S+)'  # tx. name
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
    tie(@progress, "Tie::File", $prog_rep) 
        or die "Can't manipulate progress report $prog_rep: $!";
}

# Read data from (usually two) GFF2/3 files.
foreach my $file (@files) { 
    my $no_scanned = 0;

    # Choose data-selection patterns for GFF2 or GFF3 file:
    if ( $file =~ / \.gff \z /xms ) {
        $text_pattern_1 = $text_pattern_1_GFF2;
        $text_pattern_2 = $text_pattern_2_GFF2;
    }
    if ( $file =~ / \.gff3 \z /xms ) { 
        $text_pattern_1 = $text_pattern_1_GFF3;
        $text_pattern_2 = $text_pattern_2_GFF3;
    }
    if ( $file !~ / \.gff[3]{0,1} \z /xms ) { 
        die "Can't interpret file format (unfamiliar suffix)!\n";
    }

    # Do first file-reading just to link transcripts to genes:
    open my $FILE, '<', $file or die "Can't open $file!\n";
    while (my $input = <$FILE>) { 
        chomp $input;

        # Link transcripts or CDSes to *unique* genes:
        if ( $input =~ / $text_pattern_1 /xms ) { 
            my $tx = $1;
            my $gene = $2;
            if ($tx2gene{$tx}) { 
                die "Transcript/CDS $tx already mapped to gene $tx2gene{$tx};",
                    " cannot also map to gene $gene!\n",
                    ;
            }
            $tx2gene{$tx} = $gene;

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
    close $FILE or die "Can't close $file filehandle!\n";

    # Do second file-reading to link genes (via transcripts/CDSes) to exons:
    open $FILE, '<', $file or die "Can't open $file!\n";
    while (my $input = <$FILE>) {
        chomp $input;
        if ( $input =~ / $text_pattern_2 /xms ) {
            my $chr = $1;
            my $nt1 = $2;
            my $nt2 = $3;
            my $ori = $4;
            my $phs = $5;
            my $tx  = $6;
            my $gene  = q{};

            # Report errors, but don't kill whole analysis over one error:
            if (! $tx2gene{$tx}) { 
                warn "Failed to previously link transcript $tx to gene!\n";
                warn "   and so cannot parse: $input\n";
            }
            if ( $tx2gene{$tx} ) { 
                $gene = $tx2gene{$tx};
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
            # $coord2gene_ref->{$chr}->{$nt}->{$ori}->{$phs}->{$gene} = 1;
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
                    @overlaps = grep { $_ ne $gene } 
                                keys %{ $coord2gene_ref->{$chr}->{$nt}->{$ori}->{$phs} };
                    if (@overlaps) { 
                        foreach my $overlapping (@overlaps) { 
                            # Record link in both directions!
                            $gene2equivs_ref->{$gene}->{$overlapping} = 1;
                            $gene2equivs_ref->{$overlapping}->{$gene} = 1;
                        }
                    }
                }

                # And then archive the gene's own data.
                $coord2gene_ref->{$chr}->{$nt}->{$ori}->{$phs}->{$gene} = 1;

                # Once we know that a gene has defined coordinates,
                #    *then* we list it for later recollection.
                #    We are not interested in genes w/o coordinates!
                # Also, it is useful to know what genes were from what files.
                $gene_exists_ref->{$gene} = $file;
            }
        }
    }
    close $FILE or die "Can't close $file filehandle!\n";
}

# At long last!!  A nice, simple printout:
foreach my $gene ( sort keys %{ $gene_exists_ref } ) { 
    print "$gene";
    print "\t$gene_exists_ref->{$gene}";
    foreach my $overlap (sort keys %{ $gene2equivs_ref->{$gene} } ) { 
        print "\t$overlap";
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

