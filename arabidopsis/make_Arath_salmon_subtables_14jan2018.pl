#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my %bigname2tag   = ();

my $name_type     = q{};
my $bigname_table = q{};
my @salmon_files  = ();

my $help;

GetOptions ( 
    'name_type=s'       => \$name_type,
    'bigname_table=s'   => \$bigname_table,
    'salmon_files=s{,}' => \@salmon_files,
    'help'              => \$help, 
);

if ( $help 
     or ( ( $name_type ne 'Gene' ) and ( $name_type ne 'Transcript' ) )   
     or (! -r $bigname_table ) 
     or (! @salmon_files ) 
   ) {
    die "make_Arath_salmon_subtables_14jan2018.pl\n",
        "    --name_type|-n      [either \"Gene\" or \"Transcript\"]\n",
        "    --bigname_table|-b  [big file names to simple data tags, in TSV table]\n",
        "    --salmon_files|-s   [input Salmon result files]\n",
        "    --help|-h           [print this message]\n",
        ;
}

open my $NAMES, '<', $bigname_table;
while (my $input = <$NAMES>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $big_name = $1;
        my $data_tag = $2;
        $bigname2tag{$big_name} = $data_tag;
    }
    else {
        die "In big-name table $bigname_table, cannot parse: $input\n";
    }
}
close $NAMES;

foreach my $salmon_file (@salmon_files) {
    # Sample input lines:
    # /home/bioinformatics/sepal_rnaseq_dec2016/salmon_CDS/results/Col_WT_rep1.trim_exact_3nt.salmon.cds.13jan2018/quant.genes.sf
    # /home/bioinformatics/sepal_rnaseq_dec2016/salmon_cDNA_w_UTRs/results/PDF1__FLAG_ATML1[...].13jan2018/quant.sf

    if ( $salmon_file =~ /\A \S* \/ ([^\s\/]+) \/ (quant (?:\.genes){0,1} \.sf) \z/xms) {
        my $bigname     = $1;
        my $salmon_type = $2;

        if (! exists $bigname2tag{$bigname} ) {
            die "Cannot parse big-name to tag in: $salmon_file\n";
        }
        my $tag = $bigname2tag{$bigname};

        my $out_tpm_file  = "$tag.$salmon_type.TPMs.txt";
        $out_tpm_file     = safename($out_tpm_file);

        my $out_reads_file = "$tag.$salmon_type.reads.txt";
        $out_reads_file    = safename($out_reads_file);

        open my $INFILE,  '<', $salmon_file;
        open my $OUT_TPMS, '>', $out_tpm_file;
        open my $OUT_READS, '>', $out_reads_file;

        while ( my $input = <$INFILE> ) {
            chomp $input;
            # Header line:
            # Name	Length	EffectiveLength	TPM	NumReads
            if ( $input =~ /\A (\S+) (?: \t \S+){2} \t (\S+) \t (\S+) \z/xms ) {
                my $gene_id = $1;  # 'Name'
                my $tpm     = $2;  # 'TPM'
                my $reads   = $3;  # 'NumReads'

                # Change the header line so that each data column is uniquely tagged, 
                # and thus can be merged into a single large table.
                if ( $gene_id eq 'Name' ) {
                    $gene_id = $name_type;
                    $tpm     = $tag . '_TPM'; 
                    $reads   = $tag . '_reads';
                }
                elsif ( looks_like_number($reads) ) {
                    # before rounding off decimals with 'int', add 0.5 to decimal $reads, 
                    # so that values of x.5 get rounded up to x+1
                    $reads = ($reads + 0.5);
                    $reads = int($reads);
                }
                # Print the TPMs and reads to separate data-slice files.
                print $OUT_TPMS  "$gene_id\t$tpm\n";
                print $OUT_READS "$gene_id\t$reads\n";
            }
            else { 
                die "From infile ($salmon_file), can't parse input: $input\n";
            }
        }

        close $INFILE;
        close $OUT_TPMS;
        close $OUT_READS;
    }
    else {
        die "Can't parse Salmon infile: $salmon_file\n";
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

