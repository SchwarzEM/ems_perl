#!/usr/bin/perl
#
# (c) Liisa Holm
# EMBL-EBI September 1997
#
# reference: Removing near-neighbour redundancy from large protein
#	sequence collections.  L. Holm & C. Sander, Bioinformatics 14:423-429 (1998)
#
# nrdb90.pl generates a representative set at 90 % identity threshold
#
# input: 
#	1st argument [trivial.rdb] - empty file or result from previous run
#	2nd argument [courant_nrdb] - fasta (pearson) format sequence database
#		(it is recommended to remove duplicates beforehand using NCBI's nrdb program)
# output: 
#	clean_nrdb - fasta format sequence database 
#	clusters - representatives + groupies
#	trivial.rdb - relational database
#
# work files:
#	new.rdb - sequences with undefined rep
#	old.rdb - sequences with defined rep
#	deleted.rdb - sequences not in courant_nrdb
#
$[=1;
$|=1;

# parameters should satisfy
# 	idecutoff>=1-1/ktup1
# 	fcutoff<=1-(1-idecutoff)*ktup2
#
# For example (ktup1=10, ktup2=5):
# - nrdb90 (default, 90 % identity threshold): idecutoff=0.90 fcutoff=0.50
# - nrdb95 (95 % identity): idecutoff=0.95 fcutoff=0.75
# - nrdb98 (98 % identity): idecutoff=0.98 fcutoff=0.90
# - nrdb99 (99 % identity): idecutoff=0.99 fcutoff=0.95
#

$ktup1=10;
$ktup2=5;
$fcutoff=0.5;
$idecutoff=0.9;
$gapfrac=0.05; # percentage of gaps allowed for diagonal hashing
$mintupfrac=0.02; # required tuple matches to check diagonal
$maxnkey=1500000; # memory limit
$maxislesize=2000; # island size limit
$rdbheader="\# Perl-RDB\n\# $0 $rdbfile $courant_nrdb \nseqno\trepno\tseqlen\tname\tseq\n8N\t8N\t8N\t40S\t40S\n";
$quiet=1;

chop($jobstart=`date`); 

if($#ARGV!=2) { die "Usage: $0 trivial.rdb '/data/research/bin/nrdb -l10 /data/research/db/nrdb \"//\" |' \n"; }
($rdbfile,$courant_nrdb)=@ARGV;

# in: rdbfile+fastafile; out: old.rdb new.rdb deleted.rdb [renumbered sequences]
&diffnrdb; 

$nsw=0; $nal=0;

open(IN,'<new.rdb');
while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip header
while(<IN>) { 
	chop; 
	($iprot,$dummy1,$dummy2,$dummy3,$seq)=split(/\t/); 
	$seq{$iprot}=$seq; 
	$new{$iprot}=1; 
	push(@outpool,$iprot); 
}
close(IN);
print "$#outpool new sequences from new.rdb\n";
open(IN,'<old.rdb');
while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip header
while(<IN>) { 
	($iprot,$rep)=split(/\t/); 
	next if($iprot!=$rep);
	push(@old,$iprot); 
	$rev{$iprot}=$iprot; 
}
close(IN);
print "$#old reps from old.rdb\n"; $n=0;
open(IN,'<old.rdb');
while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip header
while(<IN>) { 
	chop; 
	($iprot,$rep,$dummy2,$dummy3,$seq)=split(/\t/); 
	$seq{$iprot}=$seq; 
	if($iprot!=$rep) { $n++; $rev{$rep}.=":$iprot"; } 
}
close(IN);
print "$n groupies from old.rdb\n";

@x=@outpool; @outpool=sort { $a <=> $b } @x; undef (@x);
while($outpool[$[]) {
	$iprot=shift(@outpool);
	$seq=$seq{$iprot}; $lenseq=length($seq); 
	&cluster($iprot,$lenseq,$seq);
	if($nkey>$maxnkey||$maxsize>$maxislesize||!$outpool[$[]) {
		print "$#outpool sequences left in outpool; $#old old representatives\n";
		# map oldreps to clusters
		$n=0; undef(%oldmates);
		foreach $jprot (@old) {
			$n++; if($n=~/000$/) { print "oldmates ... $n\n" unless $quiet;}
			next if(!$rev{$jprot}); ## splice replaced oldreps from @old!
			undef(%mate);
			(@x)=split(//,$seq{$jprot}); (@w)=splice(@x,1,$ktup1-1);
			while($x[$[]) {
				push(@w,shift(@x)); $key=join('',@w); shift(@w);
				$x=$dekahash{$key};
				if($x) { $mate{$current{$x}}++; } 
			}
			foreach(keys(%mate)) { $oldmates{$_}.="$jprot "; }
			undef(%mate);
		}
		undef(%dekahash); undef(%current); undef(%tnerruc); undef(%cluster);
		undef(%size);
		foreach $r (keys(%protclus)) { 
			print "protclus $r -> $protclus{$r} oldmates -> $oldmates{$r}\n" unless $quiet;
			# pentahash oldmates
			foreach $q (split(/\s+/,$oldmates{$q})) {
				&load_m($seq{$q},$ktup2,$q,*hash);
			}
			(@pool)=split(/\s+/,$protclus{$r}); 
			undef(@rep2); undef(@newreps);
			while($newreps[$[]||$pool[$[]) {
				if($pool[$[]) { 
					(@newreps)=@rep2; undef(@rep2);
					&pool_loop; 
				}
				@x=@newreps; @newreps=sort { $a <=> $b } @x; undef(@x);
				while($newreps[$[]) { 
					$rep2=shift(@newreps);
					&mates_m($rep2);
					if(!$skip) { push(@rep2,$rep2); }
				}
			}
			undef(%hash);
			foreach $q (@rep2) { delete $new{$q}; }
			push(@old,@rep2); 
			# print "cleared reps: @rep2\n";
		}
		undef(%oldmates); undef(%protclus);
		$iclus=0; $nkey=0; $maxsize=0;
		@x=@outpool; @outpool= sort { $a <=> $b } @x; undef(@x);
	}
}

open(OUT,">clusters"); open(OUT1,">clean_nrdb");
&rdb2clusters; 
close(OUT); close(OUT1);

print "$nsw calls to diagonal; $nal calls to align\n";

print "writing to $rdbfile and clean_nrdb\n";
open(OUT,">$rdbfile");
print OUT $rdbheader;
foreach $rep (keys(%rev)) { foreach(split(/\:/,$rev{$rep})) { $rep{$_}=$rep; } }
foreach $iprot (sort { $a <=> $b } keys(%seq)) {
	$l=length($seq{$iprot}); 
	print OUT "$iprot\t$rep{$iprot}\t$l\t$name{$iprot}\t$seq{$iprot}\n";
}
close(OUT);

chop($jobend=`date`);
print "timestamp: $jobstart ... $jobend\n";

exit(1);

###############################################################################

sub getused {
	local ($seq,$ktup,*used)=@_;
	# hash seq
	(@x)=split(//,$seq);
	(@w)=splice(@x,1,$ktup-1);
	while($x[$[]) { 
		push(@w,shift(@x)); $key=join('',@w); shift(@w);
		if(!defined($used{$key})) { 
			$used{$key}=1; 
		} else {
			$used{$key}++; 
		}
	}
}

sub load_m {
	local ($seq,$ktup,$label,*hash)=@_;
	# print "this is load_m $label $ktup ",length($seq),"\n";
	undef(%used);
	&getused($seq,$ktup,*used);
	foreach $key (keys(%used)) {
		$hash{$key}.="$label:$used{$key} ";
	}
}

##############################################################################

sub mates_m {
	local($iprot)=@_; # returns ilen minmatch window minwinmatch %mate $seq
	$seq=$seq{$iprot};
	$ilen=length($seq); 
	$minmatch=$ilen*$fcutoff; # require 50 % match of tetrapeptides -> 90%ide
	undef(%mate); 
	undef(%used); 
	&getused($seq,$ktup2,*used);
	foreach $key (keys(%used)) {
		$n0=$used{$key};
		# print "$key ($n0) -> $hash{$key}\n";
		foreach(split(/\s+/,$hash{$key})) {
			($jprot,$n)=split(/\:/);
			if($n0<$n) { $n=$n0; }
			$mate{$jprot}+=$n; 
		}
	}
	# only keep mates that are above threshold! compare sorted by matecount
	undef(@array);
	foreach $jprot (keys(%mate)) {
		$n=$mate{$jprot};
		if($iprot>$jprot&&$n>=$minmatch) {
			push(@array,$jprot); 
		} elsif ($iprot<$jprot) { # use 0.5*jlen-4 exact formula
			$x=$fcutoff*length($seq{$jprot})-$ktup2+1;
			if($n>=$x) { push(@array,$jprot); }
		}
	}
	$skip=0;
	foreach $jprot (sort { $mate{$b} <=> $mate{$a} } @array) {
		# print "testing $iprot jprot $jprot [$rev{$jprot}] ($new{$jprot})\n";
		next if(!$rev{$jprot}); # has been merged before
		$skipthis=0;
		$domflag=0;
		$n=$mate{$jprot};
		$jlen=length($seq{$jprot}); 
		if($ilen<=$jlen) { ## default: query is shorter than match
			$domflag=1;
			$skipthis=&diagonal($seq{$iprot},$seq{$jprot},$ktup2,$gapfrac,$fcutoff,$mintupfrac,$idecutoff); 
		} else { ## compare longer query to shorter match
			$domflag=2;
			$skipthis=&diagonal($seq{$jprot},$seq{$iprot},$ktup2,$gapfrac,$fcutoff,$mintupfrac,$idecutoff); 
		}
		# print "matched $iprot=$rev{$iprot} $jprot=$rev{$jprot} $n $ilen $jlen $skipthis\n";
		if($skipthis) { # means merger
			if($domflag==1) {
				&dissolve($iprot,$jprot,$jlen);		
				$skip=1;
			} else {
				if(!defined($rev{$iprot})) { $rev{$iprot}="$iprot"; } # replacement inside iterated pool_loop
				&dissolve($jprot,$iprot,$ilen);
			}
		}
		last if($skip);
	}
}

sub dissolve {
	local($iprot,$jprot,$jlen)=@_;
	# merge rep2==iprot with rep1==jprot
	$rev{$jprot}.=":$iprot";
	(@groupies)=split(/\:/,$rev{$iprot}); shift(@groupies);
	print "merge $iprot ($new{$iprot}) -> $jprot ($new{$jprot}) $jlen $#groupies groupies (domflag $domflag)\n" unless $quiet;
	delete $rev{$iprot};
	# dissolve groupies of iprot -> follow jprot else ->@pool
	foreach $groupie (@groupies) {
		# if($groupie<$jprot) { ## never possible since iprot is shorter than jprot
		#	$skipthis=0; 
		#} else {
			$skipthis=&diagonal($seq{$groupie},$seq{$jprot},$ktup2,$gapfrac,$fcutoff,$mintupfrac,$idecutoff);
		#}
		if($skipthis) { 
			$rev{$jprot}.=":$groupie"; 
			print "groupie $groupie followed $iprot -> $jprot\n" unless $quiet;
		} else { 
			if($new{$groupie}) { 
				push(@pool,$groupie); 
				print "groupie $groupie of $iprot pooled ($#pool)\n" unless $quiet;
			} else { 
				push(@outpool,$groupie);
				print "groupie $groupie of $iprot outpooled ($#outpool)\n" unless $quiet;
				$new{$groupie}=1;
			} 
		}
	}
}

##############################################################################
#              diagonal identity check                                       #
##############################################################################

sub diagonal {
	local ($seq1,$seq2,$ktup,$gapfrac,$fcutoff,$mintupfrac,$idecutoff)=@_;
	local ($ilen,$jlen,$minsumtup,$minsumres,$ires,@x,@w,%hash,%diag,$d,
		@d,@diags,$maxdiag,%mat,%best,$score);
	$ilen=length($seq1);
	$jlen=length($seq2);
	$mintup=$mintupfrac*$ilen; # diagonal must match at least 2 tups per 100 aa
	$minsumtup=$fcutoff*$ilen; # seq1 is shorter one!
	$minsumres=$idecutoff*$ilen; # seq1 is shorter one!
	# print "this is subroutine diagonal($ktup,$gapfrac,$fcutoff,$mintupfrac)  $ilen $jlen $minsumtup $minsumres\n";
	$nsw++;
	# step 1: find best gapfrac diagonals by ktup hashing
	(@x)=split(//,$seq1); (@w)=splice(@x,1,$ktup-1);
	$ires=0;
	while($x[$[]) {
		$ires++;
		push(@w,shift(@x)); $key=join('',@w); shift(@w);
#		ambiguous letters always mismatch!
		next if($key=~/[BXZ]/);
		if(!defined($hash{$key})) { $hash{$key}=$ires; } 
		else { $hash{$key}.=":$ires"; }
	}
	
	(@x)=split(//,$seq2); (@w)=splice(@x,1,$ktup-1);
	$jres=0;
	while($x[$[]) {
		$jres++;
		push(@w,shift(@x)); $key=join('',@w); shift(@w); 
#		ambiguous letters always mismatch!
		next if($key=~/[BXZ]/);
		foreach(split(/\:/,$hash{$key})) {
			$d=$_-$jres;
			$diag{$d}+=1;
		}
	}
	
	(@d)=sort { $diag{$b} <=> $diag{$a} } keys(%diag);
	$maxdiag=1+$gapfrac*$jlen;
	splice(@d,$maxdiag+1);
	$s=0; $above=0; undef(@diags); 
	foreach(@d) { 
		last if($diag{$_}<$mintup);
		if($diag{$_}>$minsumres) { # single diagonal
			# print "single diagonal $diag{$_} > $minsumres\n";
			return(1); 
		} 
		$s+=$diag{$_}; $f=$s/$minsumtup;
		# print "best diagonals: $_\t$diag{$_}\t$s\t$f\t$minsumtup\n"; 
		# build up alignment matrix
		push(@diags,$_);
	}
	if($#diags>0) { 
		($above)=&align_sw; 
	}
	# print "return $above from diagonal $#diags\n";
	return($above);
}	

sub align_sw {
	local($i,$ires,$jres,$xres);
	local(@mat,@prev,@best);
	# align few best diagonals -> return score>minsum
	$nal++;
	# print "aligned diagonals @diags\n";
	(@seq1)=split(//,$seq1);
	(@seq2)=split(//,$seq2);
	foreach $ires (1..$ilen) {
		$xres=$seq1[$ires];
		foreach $shift (1..$#diags) { 
			$jres=$ires-$diags[$shift];
			next if($jres<1 || $jres>$jlen);
			if($xres ne $seq2[$jres]) {
				$mat[$ires][$shift]=-0.15; # mismatch
			} else {
				$mat[$ires][$shift]=1; # identity
			}
			# print "$mat[$ires][$shift]\t$seq2[$jres]\t";
		}
		# print "$xres\n";
	}
	$globscore=0; $globcell="0:0";
	foreach $shift (1..$#diags) { 
		$prev[1][$shift]="0:0"; 
		$best[1][$shift]=$mat[1][$shift]; 
	}
	foreach $ires (2..$ilen) {
		foreach $shift (1..$#diags) {
			$s=$mat[$ires][$shift];
			$curdiag=$diags[$shift];
			$jres=$ires-$curdiag;
			# determine best place to come from
			# 	(ires-1,[1..jres-1])  \
			# 	(ires-1,jres-1)	       ---> (ires,jres)
			# 	([1..ires-1],jres-1)  /
			$bestscore=0;
			$bestcell="0:0";
			foreach $diag (1..$#diags) {
				if($diag==$shift) { # continue match
					$xres=$ires-1;
					$score=$s+$best[$xres][$diag];
				} elsif($diags[$diag]>$curdiag) { # deletion in query
					$xres=$ires-1;
					$score=$s+$best[$xres][$diag]-1;
				} else { # insert in query
					$xres=$jres+$diags[$diag];
					next if($xres<1 || $xres>$ilen);
					$score=$s+$best[$xres][$diag]-($ires-$xres);
				}
				if($score>$bestscore) {
					$bestcell="$xres:$diag"; 
					$bestscore=$score;
				}
			}
			$prev[$ires][$shift]=$bestcell;
			$best[$ires][$shift]=$bestscore;
			if($bestscore>$globscore) {
				$globscore=$bestscore;
				$globcell="$ires:$shift";
			}
			# print "$diags[$shift]:$best[$ires][$shift]\t";
		}
		# print "$ires\t$bestcell\n";
	}
	# traceback to calculate nide
	# print "globcell: $globcell\n";
	$nide=0;
	($ires,$shift)=split(/\:/,$globcell);
	while($ires>0) {
		if($mat[$ires][$shift]>0) { $nide++; }
		# print "$ires:$shift\t$nide\n";
		($ires,$shift)=split(/\:/,$prev[$ires][$shift]);
	}

	if($nide>$minsumres) { $x=1; } else { $x=0; }
	# print "score = $nide / $minsumres\treturn $x from align\n";
	return ($x);
}

##############################################################################
##############################################################################

sub cluster {
	local($iprot,$lenseq,$seq)=@_;
	#
	# dekapeptide dekahash
	#
        # ktups of query sequence
        (@x)=split(//,$seq); (@w)=splice(@x,1,$ktup1-1);
	undef(%merge);
	$iclus++;
        while($x[$[]) {
                push(@w,shift(@x)); ($key)=join('',@w); shift(@w);
                next if($key=~/[BXZ]/); # ignore ambiguity codes
		# lookup islands linked by query protein
		$xclus=$dekahash{$key};
		if(!$xclus) {
			$dekahash{$key}=$iclus;
			$nkey++;
		} else {
			next if($xclus==$iclus);
			$merge{$current{$xclus}}=1;
		}
        }
	(@merge)=sort {$a <=> $b} keys(%merge);
	# print "$nkey words; to merge clusters @merge with $iclus\n";
	$kclus=$merge[$[];
	if($kclus) { 
		# merge linked clusters
		foreach $jclus (@merge[$[+1..$#merge]) {
			# print "* * * jclus->current: $jclus, $current{$jclus}, $tnerruc{$jclus} kclus: $kclus <- $tnerruc{$kclus}\n";
			foreach(split(/\s+/,$protclus{$jclus})) { $cluster{$_}=$kclus; }
			$protclus{$kclus}.=$protclus{$jclus};
			$size{$kclus}+=$size{$jclus};
			delete $protclus{$jclus};
			# update short-link to current cluster
			foreach(split(/\s+/,$tnerruc{$jclus})) { 
				$current{$_}=$kclus; 
				# print "* * * current $_ -> $kclus\n";
			}
			$tnerruc{$kclus}.=$tnerruc{$jclus};
			delete $tnerruc{$jclus};
			$size{$jclus}=0;
		}
		# append iprot
		# print "* append $iprot -> $kclus\n";
		$cluster{$iprot}=$kclus;
		$protclus{$kclus}.="$iprot ";
		$current{$iclus}=$kclus;
		$tnerruc{$kclus}.="$iclus ";
		$size{$kclus}+=1;
		if($size{$kclus}>$maxsize) { $maxsize=$size{$kclus}; }
	} else  { # iprot is new unique
		# print "* $kclus * new $iprot -> $iclus\n";
	        $cluster{$iprot}=$iclus; 
	     	$protclus{$iclus}="$iprot ";
		$current{$iclus}=$iclus;
		$tnerruc{$iclus}="$iclus ";
		$size{$iclus}=1;
	}
}

##############################################################################
##############################################################################

sub diffnrdb {
	local(%ix,%xi,%rep,%hash,%newix,$oldix,%newrep,%oldrep,%name,%lenseq);
	# read in new nrdb -> %name %seq %lenseq %hash
	print "reading new: $courant_nrdb\n" unless $quiet;
	open(IN,"$courant_nrdb");
	$nprot=0; $t=0; $m=0;
	while(<IN>) {
	        chop;
	        if(/^>/) { 
			# save old entry
			if($nprot>0) { &save_entry; }
	                ($name)=/^>(.*)/; 
	                $seq=''; 
	                $nprot++; 
			$ix{$nprot}=$nprot;
	                if($nprot=~/0000$/) { print "new...$nprot\n" unless $quiet; }
	        }
	        else { $seq.=$_; }
	}
	close(IN);
	# save last entry
	&save_entry;
	print "$nprot sequences saved; $t keys; $m multiple hash key occupancies\n";

	# read in old nrdb -> %newix{oldseqno}=new-seqno %oldix{new-seqno}=oldseqno
	if(-e $rdbfile) {
	  $ndel=0;
	  print "reading old: $rdbfile\n" unless $quiet;
	  open(IN,"<$rdbfile");
	  open(OUT,">deleted.rdb"); print OUT $rdbheader;
	  while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip rdb header
	  $n=0;
	  while(<IN>) {
		chop($line=$_);
		$n++;
		($i,$rep,$l,$name,$seq)=split(/\t+/,$line);
		$rep{$i}=$rep;
		if($n=~/0000$/) { print "old...$n\n" unless $quiet; }
		($checksum)=&do_checksum($l,$seq);
		(@x)=split(/\:/,$hash{$checksum});
		# print "old-new $i $l $checksum vs. @x \n";
		$keep=0;
		foreach(@x) {
			if($seq eq $xseq{$_}) {
				$newix{$i}=$_;
				$oldix{$_}=$i;
				$keep=1;
				# print "old $i ==  $_ ($l) <=> ($lenseq{$_})\n";
				last;
			}
		}
		if(!$keep) { 
			$ndel++;
			print OUT "$i\t$rep\t$l\t$name\t$seq\n";
		}
	  }
	  close(IN); close(OUT);
	  print "$ndel deleted entries in deleted.rdb\n";
	}

	# transfer reps
	foreach $ix (keys(%lenseq)) {# ix keys=courant-seqno
		$oldix=$oldix{$ix};  
		if($oldix) { $oldrep=$rep{$oldix}; } else { $oldrep=0; }
		if($oldrep) { $newrep{$ix}=$newix{$oldrep}; } else { $newrep{$ix}=0; }
		# print "transfer $ix ($lenseq{$ix})-> $oldix -> $oldrep -> $newix{$oldrep} \n";
	}
	
	# renumber sorted by length -> %ix{courant-seqno}=sorted-seqno
	foreach $i (keys(%lenseq)) {
		if($newrep{$i}==$i) { $lenseq{$i}+=0.2; }
		elsif($newrep{$i}>0) { $lenseq{$i}+=0.1; }
	}
	$i=0;
	foreach (sort { $lenseq{$b} <=> $lenseq{$a} } keys(%lenseq) ) {
		$i++;
		$ix{$_}=$i;
		$xi{$i}=$_;
		# print "sorted $_ -> $i $lenseq{$_} \n";
	}

	# write renumbered sequences/reps to new.rdb, old.rdb
	print "writing renumbered sequences/reps to new.rdb, old.rdb\n";
	open(OUT,">old.rdb");
	open(OUT1,">new.rdb");
	$nnew=0; $nold=0;
	print OUT $rdbheader; print OUT1 $rdbheader;
	foreach $i (sort { $a <=> $b } keys(%ix) ) {
		$ix=$xi{$i};
		$_=$lenseq{$ix}; ($l)=/^(\d+)/;
		if($newrep{$ix}) {
			$nold++;
			print OUT "$i\t$ix{$newrep{$ix}}\t$l\t$name{$ix}\t$xseq{$ix}\n";
		} else {
			$nnew++;
			print OUT1 "$i\t$ix{$newrep{$ix}}\t$l\t$name{$ix}\t$xseq{$ix}\n";
		}
	}
	close(OUT); close(OUT1);
	print "$nnew new sequences in new.rdb\n$nold old sequences in old.rdb\n";
	undef(%ix); undef(%xi);
	undef(%rep); undef(%hash);
	undef(%newix); undef(%oldix);
	undef(%newrep); undef(%oldrep);
	undef(%xseq); undef(%name); undef(%lenseq);
}


sub save_entry {
	$xseq{$nprot}=$seq; 
	$name{$nprot}=$name;
	$l=length($seq);
	$lenseq{$nprot}=$l;
	($checksum)=&do_checksum($l,$seq);
	if(!defined($hash{$checksum})) { $hash{$checksum}="$nprot"; $t++; } 
	else { $hash{$checksum}.=":$nprot"; $m++; }
	# print "saved $nprot\t$checksum\t$hash{$checksum}\t$l\t$name\n";
}

sub do_checksum {
	local ($l,$seq)=@_;
	local($checksum,$k);
	$checksum=0;
	$k=$l/20; if($k<2) { $k=2; }
	while($seq) {
		$checksum+=unpack("%32C*",$seq);
		$seq=substr($seq,$k);
	}
	$checksum%=429496796; # 2^32;
	$checksum.=":$l";
	return($checksum);
}

###############################################################################
###############################################################################

sub pool_loop { # redefine representatives for dissolved groupies in set2
	# print "this is pool_loop dealing with $#pool sequences and $#rep2 rep2s: @rep2\n";
	# pentahash rep2s within island2
	local(%hash);
	undef(%hash);
	foreach $rep2 (@rep2) { # "clean" reps passed through rep_loop
		next if(!$rev{$rep2}); # has been merged before
		&load_m($seq{$rep2},$ktup2,$rep2,*hash);
	}
	# sort pool by length
	@x=@pool; @pool=sort { $a <=> $b } @x;
	while($pool[$[]) {
		$iprot=shift(@pool);
		&mates_m($iprot);
		# unique iprot within island2
		if(!$skip) {  # "clean" rep2 -> add to hash table!
			push(@newreps,$iprot);
			$rev{$iprot}="$iprot"; 
			# print "new rep $iprot from pool_loop $#newreps newreps\n" unless $quiet;
			&load_m($seq{$iprot},$ktup2,$iprot,*hash);
		}
	}
	undef(%hash);
	return;
}

##############################################################################
##############################################################################

sub rdb2clusters {
	open(IN,'<new.rdb');
	while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip rdb header
	while(<IN>) {
		chop;
		($iprot,$dummy1,$l,$name)=split(/\t/);
		$lenseq{$iprot}=$l;
		$name{$iprot}=$name;
	}
	close(IN);
	open(IN,'<old.rdb');
	while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip rdb header
	while(<IN>) {
		chop;
		($iprot,$dummy1,$l,$name)=split(/\t/);
		$lenseq{$iprot}=$l;
		$name{$iprot}=$name;
	}
	close(IN);
	
	$nclus=0;
	foreach $rep (sort { $rep{$a} <=> $rep{$b} } keys(%rev)) {
		(@prot)=sort { $a <=> $b } split(/\:/,$rev{$rep});
		$nclus++;
		# summarize information from names
		$info=&expert(@prot);
		print OUT "Cluster $nclus with $#prot members $rep $lenseq{$rep} $name{$rep} $info\n";
		print OUT1 "\>$name{$rep} $info\n$seq{$rep}\n";
		foreach(@prot) {
			print OUT "$nclus\t$_\t$lenseq{$_}\t$name{$_}\n";
		}
	}
	
}


sub expert {
	local(@prot)=@_;
	$size=$#prot;
	shift(@prot); # remove rep
	# aim: select informative annotation from cluster groupies
	$info='';
	# report first swissprot or pdb
	foreach(@prot) { if($name{$_}=~/\:swiss\|/) { $info.=$name{$_}; last; } }
	foreach(@prot) { if($name{$_}=~/\:pdb\|/) { $info.=$name{$_}; last; } }
	# exclude junk
	if(!$info) {
		foreach(@prot) {
			$x=$name{$_}; ($x)=~tr/[A-Z]/[a-z]/;
			next if($x=~/hypothetical/);
			next if($x=~/homolog/);
			next if($x=~/ orf /);
			next if($x=~/putative/);
			next if($x=~/open reading frame/);
			next if($x=~/ polyprotein /);
			next if($x=~/ cosmid /);
			next if($x=~/product\: \"[a-z]+[0-9]+\"/);
			$info.=$name{$_}; last;
		}
	}
	if($info) { return("\/\/info(N\=$size): $info"); } else { return(''); }
}

##############################################################################
##############################################################################


