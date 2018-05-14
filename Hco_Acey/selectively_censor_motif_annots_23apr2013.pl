#!/usr/bin/env perl

use strict;
use warnings;

my $header = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.*) \z/xms ) { 
        my $gene       = $1;
        my $annot_text = $2;

        # Treat the very first line as the header, and print the header as-is once.
        if (! $header) { 
            $header = $input;
            print "$header\n";
        }

        # Treat all subsequent non-header lines as annotation text to be censored if a single Pfam-A or InterPro motif occurs.
        else { 
            my @unfiltered_annots = split /\t/, $annot_text;

            # Note that there should be one less \t than the annots which the \t-s are separating, so this index need not be changed to work for an array.
            my $unfiltered_annot_index = $annot_text =~ tr/\t/\t/;

            my $censor = 0;
            # A single instance of Pfam-A or InterPro disqualifies *all* non-Pfam-A/InterPro annots for a given gene.
            if ( ( $annot_text =~ /PF\d+\.\d+/xms ) or ( $annot_text =~ /IPR\d+/xms ) ) {
                $censor = 1;
            }

            # Build, as appropriate, either an uncensored or a censored annotation line; every gene gets either its Pfam-A or InterPro motifs,
            #    or annotations which consist solely of non-Pfam-A/InterPro motifs *if* that is truly all that gene has.
            # Subsequent Wilcoxon ranking will only 'see' non-Pfam-A/InterPro motifs which do not co-occur with Pfam-A/InterPro.
            # The systematic error, here, is that highly ubiquitous and uninformative Pfam-A/InterPro motifs can mask legitimately interesting non-Pfam-A/InterPro motifs.
            # But anything which does make it through this process and ends up getting a significant Wilcoxon ranking really *will* be a novelty.

            my $revised_annot_text = q{};

            # I append text, rather than building text from a \t-joined array, because getting Perl to parse arrays with q{} elements is annoyingly crashy.
            foreach my $unfiltered_annot_val (0..$unfiltered_annot_index) {
                # Default, which Perl will let me initialize *if* I do so explicitly:
                my $unfiltered_annot = q{};

                # However, if there is some actual text there, use that instead of q{}:
                if ( $unfiltered_annots[$unfiltered_annot_val] ) { 
                    $unfiltered_annot = $unfiltered_annots[$unfiltered_annot_val];
                }

                # But wait!  If censorship is in force, only allow Pfam-A and InterPro annotations to exist; reduce all else to q{} again.
                if ( $censor and ( ( $unfiltered_annot !~ /\[PF\d+\.\d+\]/xms ) and ( $unfiltered_annot !~ /\[IPR\d+\]/xms ) ) ) {
                    $unfiltered_annot = q{};
                }

                # Finally, append whatever is left.
                $revised_annot_text .= "\t$unfiltered_annot";
            }
            # No need to specify "\t" between the two scalars, since $revised_annot_text starts with a "\t":
            print "$gene$revised_annot_text\n";
        }
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

