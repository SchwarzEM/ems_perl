#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $data_ref;

my $annots   = q{};
my $table    = q{};
my @deg_sets = ();
my $help;

my %allowed_oddity = (
    'Gene' => 1,
    'Athal_18S_RRNA' => 1,
    'Athal_25S_RRNA' => 1,
    'Athal_5.8S_RRNA' => 1,
    'GFP_coding_seq' => 1,
    'PE_1.0_adaptor' => 1,
    'PE_2.0_adaptor' => 1,
);

GetOptions ( 'annots=s'      => \$annots,
             'table=s'       => \$table,
             'deg_sets=s{,}' => \@deg_sets,
             'help'          => \$help,     );

if ($help or (! $annots) or (! $table) or (! @deg_sets) ) { 
    die "Format: split_sepal_go_terms_24feb2016.pl\n",
        "            -a|--annots    [complex gene annotation table]\n",
        "            -t|--table     [GO table to expand with genes]\n",
        "            -d|--deg_sets  [one or more lists of genes sig. assoc. w/ changes]\n",
        "            -h|--help      [print this message]\n",
        "            [print output to <STDOUT>]\n",
        ;
}

# Get genes that we denoted as *significantly* changed in an up or down condition.
foreach my $deg_set (@deg_sets) {
    # Get a plain-text description of what the DEG set is about.
    my $condition = basename($deg_set);

    # Revise the plain-text desc.
    if ( $condition =~ /\A (\S+reg) \.gene_list\.txt \z/xms ) { 
        $condition = $1;
        $condition =~ s/ATML1__LGO/LGOoe/g; 
        $condition =~ s/_/ /g;

        # Mostly, the above is fine, but one term needs to *keep* '_':
        $condition =~ s/Col WT/Col_WT/g;

        $condition =~ s/\.vs\./ vs. /g; 
        $condition =~ s/\.upreg/, upreg./g; 
        $condition =~ s/\.downreg/, downreg./g; 
    }
    else {
        die "Cannot extract DEG type from name of $deg_set\n";
    }

    open my $DEG, '<', $deg_set;
    while (my $input = <$DEG>) {
        chomp $input;
        if ( $input =~ / \A (AT (?: 1|2|3|4|5|C|M) G\d+) \z/xms ) {
            my $tair_gene = $1;
            $data_ref->{'condition'}->{$condition}->{'tair_gene'}->{$tair_gene} = 1;
        }
        else {
            die "From DEG set list $deg_set, cannot parse gene name $input\n";
        }
    }
    close $DEG;
}

open my $ANNOTS, '<', $annots;
while (my $input = <$ANNOTS>) {
    chomp $input;
    if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+) \t ([^\t]*) (?: \t [^\t]* ){5} \t ([^\t]*) ( \t [^\t]* ){3} \t /xms ) { 
        my $tair_gene = $1;
        my $aliases   = $2;
        my $go_text   = $3;
        my $other_go  = $4;

        $aliases =~ s/ \[ [^\[\]]+ \] //gxms;
        $aliases =~ s/[;]/|/g;
        $aliases =~ s/\s//g;
        my $name = $tair_gene;
        if ( $aliases =~ /\S/xms ) {
            $name = "$name|$aliases";
        }

        # Count *all* genes in the genome that are annotated for *any* GO terms.
        if ( ( $go_text =~ / GO:\d+ /xms ) or ( $other_go =~ / GO:\d+ /xms ) ) {
            $data_ref->{'gene_w_go_term'}->{$name} = 1;
        }

        if ( $go_text =~ /GO:\d+/ ) {
            my @go_annots = split /; /, $go_text;
            foreach my $go_annot (@go_annots) {
                if ( $go_annot =~ /\A (\S.* [ ] \[GO:\d+\]) [ ] \( (.+) \) \z/xms ) {
                    my $go_term   = $1;
                    my $cond_text = $2;

                    # Get *total* number of genes annotated to a given GO term, regardless of conditions.
                    $data_ref->{'go_term'}->{$go_term}->{'total_annot'}->{'gene'}->{$name} = 1;

                    if ( $cond_text ne 'not sig. gene' ) { 
                        # mask ',' in '[condition], upreg.' and '[condition], downreg.';
                        # then do split; then put ',' for 'upreg.' and 'downreg.' back.
                        $cond_text =~ s/, upreg\./XXCOMMAXX upreg./g;
                        $cond_text =~ s/, downreg\./XXCOMMAXX downreg./g;
                        my @conditions = grep { $_ !~ /batch/ } split /, /, $cond_text;
                        foreach my $condition (@conditions) {
                            $condition =~ s/XXCOMMAXX upreg\./, upreg./g;
                            $condition =~ s/XXCOMMAXX downreg\./, downreg./g;
                            if (! exists $data_ref->{'condition'}->{$condition}->{'tair_gene'} ) {
                                die "From annots file $annots, cannot recognize nominal condition \"$condition\" in: $input\n";
                            }
                            if ( $condition !~ /batch/xms ) { 
                                $data_ref->{'go_term'}->{$go_term}->{'condition'}->{$condition}->{'gene'}->{$name} = 1;
                                $data_ref->{'obs_condition'}->{$condition} = 1;
                            }
                        }
                    }
                }
                else {
                    die "From annots file $annots, cannot parse GO annot \"$go_annot\" in: $input\n";
                }
            }
        }
    }
    elsif ( ( $input !~ /\A \S+ \t /xms ) or ( ( $input =~ /\A (\S+) \t /xms ) and (! exists $allowed_oddity{$1} ) ) ) {
        die "From annots file $annots cannot parse: $input\n";
    }
}
close $ANNOTS;

my @go_annot_genes      = sort keys %{ $data_ref->{'gene_w_go_term'} };
my $go_annot_gene_count = @go_annot_genes;

open my $TABLE, '<', $table;
while (my $input = <$TABLE>) {
    chomp $input;
    if ( $input eq "GO term\tp-value [condition]" ) {
        print "GO term\tp-value [condition]\tAll GO annots.\tCond. GO annots.\tGene count\tGenes\n";
    }
    elsif ( $input =~ /\A ([^\t]+) \t ([^\t]+) \z/xms ) {
        my $go_term  = $1;
        my $p_values = $2;

        my @all_genes      = ();
        my $all_gene_count = 0;

        my @conditions     = ();
        my @cond_all_genes = ();
        my @cond_go_genes  = ();

        # For some reason, there are GO terms showing up in downstream analyses from FUNC that are not in TAIR GO!
        # So, make the following steps use 'warn', not 'die'.
        if (! exists $data_ref->{'go_term'}->{$go_term}->{'total_annot'}->{'gene'} ) {
            warn "From GO table $table, cannot map GO term \"$go_term\" to genes\n";
        }
        if (! exists $data_ref->{'go_term'}->{$go_term}->{'condition'} ) {
            warn "From GO table $table, cannot map GO term \"$go_term\" to conditions\n";
        }   

        if (     ( exists $data_ref->{'go_term'}->{$go_term}->{'total_annot'}->{'gene'} ) 
             and ( exists $data_ref->{'go_term'}->{$go_term}->{'condition'}             ) ) {

            @all_genes      = sort keys %{ $data_ref->{'go_term'}->{$go_term}->{'total_annot'}->{'gene'} };
            $all_gene_count = @all_genes;

            @conditions     = sort keys %{ $data_ref->{'go_term'}->{$go_term}->{'condition'} };
            @cond_all_genes = ();
            @cond_go_genes  = ();
        }

        foreach my $condition (@conditions) {
            my @new_cond_all_genes = sort keys %{ $data_ref->{'condition'}->{$condition}->{'tair_gene'} };
            my @new_cond_go_genes  = sort keys %{ $data_ref->{'go_term'}->{$go_term}->{'condition'}->{$condition}->{'gene'} };
            push @cond_all_genes, @new_cond_all_genes;
            push @cond_go_genes, @new_cond_go_genes ;
        }

        @cond_all_genes = sort @cond_all_genes;
        @cond_all_genes = uniq @cond_all_genes;

        @cond_go_genes = sort @cond_go_genes;
        @cond_go_genes = uniq @cond_go_genes;

        my $cond_all_gene_count = @cond_all_genes;

        my $cond_go_gene_count = @cond_go_genes ;
        my $cond_go_gene_text  = join '; ', @cond_go_genes ;

        my $all_annot_ratio = 'n/a';
        my $sig_annot_ratio = 'n/a';
        if ( $go_annot_gene_count >= 1 ) {
            $all_annot_ratio = ( ( $all_gene_count / $go_annot_gene_count ) );
        }
        if ( $cond_all_gene_count >= 1 ) {
            $sig_annot_ratio = ( ( $cond_go_gene_count / $cond_all_gene_count ) );
        }

        my $go_annot_gene_no = commify($go_annot_gene_count);
        my $all_gene_no      = commify($all_gene_count) ;
        my $cond_all_gene_no = commify($cond_all_gene_count);
        my $cond_go_gene_no  = commify($cond_go_gene_count);

        my $all_summary  = "$all_gene_no/$go_annot_gene_no [$all_annot_ratio]";
        my $cond_summary = "$cond_go_gene_no/$cond_all_gene_no [$sig_annot_ratio]";

        print "$go_term\t$p_values\t$all_summary\t$cond_summary\t$cond_go_gene_count\t$cond_go_gene_text\n";
    }
    else {
        die "From table $table, cannot parse line: $input\n";
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

