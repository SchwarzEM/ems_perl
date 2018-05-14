#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @scripts = qw(
    job_rsem_sepals_atml1_3_rep1.trim_exact_3nt_2016.01.21.01
    job_rsem_sepals_atml1_3_rep2.trim_exact_3nt_2016.01.21.10
    job_rsem_sepals_atml1_3_rep3.trim_exact_3nt_2016.01.21.17
    job_rsem_sepals_atml1_4_rep1.trim_exact_3nt_2016.01.21.02
    job_rsem_sepals_ATML1__LGO_atml1_3_rep1.trim_exact_3nt_2016.01.21.03
    job_rsem_sepals_ATML1__LGO_atml1_3_rep2.trim_exact_3nt_2016.01.21.11
    job_rsem_sepals_ATML1__LGO_atml1_3_rep3.trim_exact_3nt_2016.01.21.18
    job_rsem_sepals_ATML1__LGO_atml1_4_rep1.trim_exact_3nt_2016.01.21.04
    job_rsem_sepals_ATML1__LGO_rep1.trim_exact_3nt_2016.01.21.05
    job_rsem_sepals_ATML1__LGO_rep2.trim_exact_3nt_2016.01.21.12
    job_rsem_sepals_ATML1__LGO_rep3.trim_exact_3nt_2016.01.21.19
    job_rsem_sepals_Col_WT_rep1.trim_exact_3nt_2016.01.21.06
    job_rsem_sepals_Col_WT_rep2.trim_exact_3nt_2016.01.21.13
    job_rsem_sepals_Col_WT_rep3.trim_exact_3nt_2016.01.21.20
    job_rsem_sepals_lgo_2_rep1.trim_exact_3nt_2016.01.21.07
    job_rsem_sepals_lgo_2_rep2.trim_exact_3nt_2016.01.21.14
    job_rsem_sepals_lgo_2_rep3.trim_exact_3nt_2016.01.21.21
    job_rsem_sepals_PDF1__FLAG_ATML1_lgo_2_rep1.trim_exact_3nt_2016.01.21.08
    job_rsem_sepals_PDF1__FLAG_ATML1_lgo_2_rep2.trim_exact_3nt_2016.01.21.15
    job_rsem_sepals_PDF1__FLAG_ATML1_lgo_2_rep3.trim_exact_3nt_2016.01.21.22
    job_rsem_sepals_PDF1__FLAG_ATML1_rep1.trim_exact_3nt_2016.01.21.09
    job_rsem_sepals_PDF1__FLAG_ATML1_rep2.trim_exact_3nt_2016.01.21.16
    job_rsem_sepals_PDF1__FLAG_ATML1_rep3.trim_exact_3nt_2016.01.21.23
);

my @system_ids = qw(
    29522317
    29522318
    29522319
    29522320
    29522321
    29522322
    29522323
    29522324
    29522325
    29522329
    29522330
    29522333
    29522337
    29522345
    29522351
    29522362
    29522365
    29522367
    29522370
    29522373
    29522374
    29522375
    29522376
);

my $count = @system_ids;
$count--;

foreach my $i (0..$count) {
    my $script = $scripts[$i];
    my $sys_id = $system_ids[$i];
    my $report = "/mnt/home/emsch/work/2015/adrienne/misc_docs/$script.e$sys_id";
    if (! -r $report ) {
        die "Unreadable putative report: $report\n";
    }
    print "$report\n";
}

