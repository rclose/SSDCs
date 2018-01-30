# quorumSubsample 4.3 by John Alroy
# see end of file for version history
use Getopt::Long; 
# use strict;
# use warnings;

my $usage = "usage:quorumSubsample4-2.pl [options] <arguments...>
options:
--PATH <string>
--FILES <string>
--TRIALS <number>
--METHOD <string>
--QUORUM <string>
--RANK <string>
--SCALE <string>
--MAXMA <number>
--EXACT <string>
--ABUND <string>
--BYCOLLECTION <string>
--SINGLETONS <string>
--EXCLUDE <string>
--BIGGEST <string>
--SEQUENTIAL <string>
--COLLSPERREF <number>
--REFQUOTA <number>
--DISPERSE <string>
--INTERPOLATE <string>
--USEFAILED <string>
--TRIM <string>
--MINCOLLBYBIN <number>
--MATCHMAX <string>
--DEORPHAN <string>
";

# location of input and output data files
my $PATH = "./";

# if multiple files are listed, subsampled diversity for each one will
#  be reported both jointly and separately
my @FILES = ();

# number of subsampling trials
my $TRIALS = 1000;

# options CR, UW, or O2W (default SQS)
my $METHOD = "SQS";

# quorum level
my $QUORUM = 0.4;

# CR, UW or O2W quota
my $QUOTA = 100;

# taxonomic rank (required, options genus or species)
my $RANK = "species";

# time scale (options stages, subepochs, FR2, Peters, a file with numbers,
#  or a number; default is 11 Myr bins)
my $SCALE = "timebins.txt";

# a number in Ma is required if SCALE is numerical
my $MAXMA = "650";

# draw all occurrences and count taxa seen at the quorum level(s)
#  (option no, default yes)
my $EXACT = "yes";

# subsample abundance data (only works for CR and SQS; default no)
my $ABUND = "no";

# draw entire collections (option no, default yes)
my $BYCOLLECTION = "yes";
# definition of singletons (option one reference and not recommended,
#  default one occurrence)
my $SINGLETONS = "occurrence";
# exclude most common (dominant) genus (option yes, default no)
my $EXCLUDE = "no";
# exclude taxa in most diverse collection from singleton count
#  (option yes, default no)
my $BIGGEST = "no";

# draw all collections in a reference before going to the next one (optional
#  and recommended)
my $SEQUENTIAL = "yes";
# fixed number of collections drawn per reference (default no limit,
#  but recommended)
my $COLLSPERREF = 3;
# minimum number of references per bin (default none, and not recommended)
my $REFQUOTA = 0;
# disperse sampling among references with a throwback algorithm (optional
#  and not recommended)
my $DISPERSE = "no";
# use interpolated boundary estimates in PaleoDB download file (default no,
#  but recommended)
my $INTERPOLATE = "no";

# print data for bins sometimes under quota (default no)
my $USEFAILED = "no";
# print data only for bins in the range of those that meet the subsampling
#  target (default no)
my $TRIM = "no";

# number of collections a bin must include to be analyzed
my $MINCOLLBYBIN = 1;
# intervals must include as many collections as the maximum drawn
#  (default no, and not recommended)
my $MATCHMAX = "no";
# assign collections spanning multiple bins to the ones including more than
#  half of their age estimate limits (default no; options yes or a fraction)
my $DEORPHAN = "no"; 
if ( $DEORPHAN > 0 && $DEORPHAN < 0.5 || $DEORPHAN > 1 )	{
	print "\nExiting because DEORPHAN must be 'yes' or a number between 0.5 and 1.0\n\n";
	exit;
}
if ( $DEORPHAN =~ /y/i )	{
	$DEORPHAN = 0.5;
}


GetOptions (
	"PATH=s" => \$PATH, 
	"FILES=s" => \@FILES, 
	"TRIALS=i" => \$TRIALS, 
	"METHOD=s" => \$METHOD,
	"QUORUM=s" => \$QUORUM, 
	"QUOTA=s" => \$QUOTA,
	"RANK=s" => \$RANK, 
	"SCALE=s" => \$SCALE,
	"MAXMA=i" => \$MAXMA,
	"EXACT=s" => \$EXACT,
	"ABUND=s" => \$ABUND,
	"BYCOLLECTION=s" => \$BYCOLLECTION,
	"SINGLETONS=s" => \$SINGLETONS,
	"EXCLUDE=s" => \$EXCLUDE,
	"BIGGEST=s" => \$BIGGEST,
	"SEQUENTIAL=s" => \$SEQUENTIAL,
	"COLLSPERREF=i" => \$COLLSPERREF,
	"REFQUOTA=i" => \$REFQUOTA,
	"DISPERSE=s" => \$DISPERSE,
	"INTERPOLATE=s" => \$INTERPOLATE,
	"USEFAILED=s" => \$USEFAILED,
	"TRIM=s" => \$TRIM,
	"MINCOLLBYBIN=i" => \$MINCOLLBYBIN,
	"MATCHMAX=s" => \$MATCHMAX,
	"DEORPHAN=s" => \$DEORPHAN
);

print "OPTIONS:\n";
print "path: $PATH\n";
print "files: @FILES\n";
print "trials: $TRIALS\n";
print "method: $METHOD\n";
if ( $METHOD eq "SQS" ) {
	print "quorum: $QUORUM\n";	
} elsif ($METHOD eq "CR" || $METHOD eq "UW" || $METHOD eq "O2W") {
	print "quota: $QUOTA\n";
	print "abund: $ABUND\n";
}
print "rank: $RANK\n";
print "scale: $SCALE\n";
print "maxma: $MAXMA\n";
print "exact: $EXACT\n";
print "bycollection: $BYCOLLECTION\n";
print "singletons: $SINGLETONS\n";
print "exclude: $EXCLUDE\n";
print "biggest: $BIGGEST\n";
print "sequential: $SEQUENTIAL\n";
print "collsperref: $COLLSPERREF\n";
print "refquota: $REFQUOTA\n";
print "disperse: $DISPERSE\n";
print "interpolate: $INTERPOLATE\n";
print "usefailed: $USEFAILED\n";
print "trim: $TRIM\n";
print "mincollbybin: $MINCOLLBYBIN\n";
print "matchmax: $MATCHMAX\n";
print "deorphan: $DEORPHAN\n";




if ( $RANK !~ /^(genus|species)$/ )	{
	print "\nExiting because RANK must be 'species' or 'genus'\n\n";
	exit;
}


if ( $SCALE !~ /[^0-9]/ && $MAXMA == 0 )	{
	print "\nExiting because MAXMA has not been set\n\n";
	exit;
}

if ( $EXACT !~ /n/i && $SINGLETONS =~ /r/i )	{
	print "\nWarning: SINGLETONS is forced to equal 'occurrences' if you choose the EXACT option\n";
}

if ( $EXACT !~ /n/i && $DISPERSE =~ /y/i )	{
	print "\nWarning: DISPERSE doesn't work if you choose the EXACT option\n";
}


my ($nbin,%lookup,@binlist,@base,@top,$numeric);
if ( ! $SCALE )	{
	%BINS = ("Cambrian" => 4, "Ordovician" => 5, "Silurian" => 2, "Devonian" => 5, "Carboniferous" => 5, "Permian" => 4, "Triassic" => 4, "Jurassic" => 6, "Cretaceous" => 8, "Cenozoic" => 6);

	for my $p ( "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Cenozoic" )	{
		for my $b ( 1..$BINS{$p} )	{
			$nbin++;
			$lookup{$p." ".$b} = $nbin;
			my $i = $p;
			if ( $i =~ /^C/ )	{
				$i =~ s/Cambrian/Cm/;
				$i =~ s/Carboniferous/C/;
				$i =~ s/Cretaceous/K/;
				$i =~ s/Cenozoic/Cz/;
			} else	{
				$i =~ s/(.).*/$1/;
			}
			$binlist[$nbin] = $i.$b;
		}
	}
} elsif ( $SCALE !~ /stages|epochs|fr2|peters/i && $SCALE =~ /[A-Za-z].*[^0-9]/ )	{
	open IN,"<./$SCALE";
	while (<IN>)	{
		$nbin++;
	}
	close IN;
	# assumes there is no header, columns are tab-delimited, the first
	#  column gives the base of the bin in Ma, and the second gives the
	#  bin name
	my $id = $nbin;
	open IN,"<./$SCALE";
	while (<IN>)	{
		s/\n//;
		my @words = split /\t/,$_;
		$base[$id] = sprintf "%.3f",$words[0];
		$top[$id-1] = $base[$id];
		if ( $words[1] )	{
			$binlist[$id] = $words[1];
		} else	{
			$binlist[$id] = $nbin;
		}
		$id--;
	}
	close IN;
	# rename the scale to allow tests for numeric values later
	$SCALE = "custom";
} elsif ( $SCALE =~ /stages|epochs|fr2|peters/i )	{
	if ( $SCALE =~ /stages/ )	{
		open IN,"<./Gradstein_stages";
	} elsif ( $SCALE =~ /epochs/ )	{
		open IN,"<./Gradstein_subepochs";
	} elsif ( $SCALE =~ /fr2/i )	{
		open IN,"<./FR2_bins";
	} elsif ( $SCALE =~ /peters/ )	{
		open IN,"<./Peters_stages";
	}
	while (<IN>)	{
		s/\n//;
		s/"//g;
		$nbin++;
		$lookup{$_} = $nbin;
		$binlist[$nbin] = $_;
	}
	close IN;
	# intervals are in reverse order
	my @temp = @binlist;
	for my $i ( 0..$#binlist )	{
		$binlist[$i] = $temp[$nbin - $i + 1];
	}
	for my $b ( keys %lookup )	{
		$lookup{$b} = $nbin - $lookup{$b} + 1;
	}
} elsif ( $SCALE > 0 )	{
	$numeric = 1;
	$nbin = int( $MAXMA / $SCALE );
	for my $i ( 1..$nbin )	{
		$base[$nbin - $i + 1] = $i * $SCALE;
		$top[$nbin - $i + 1] = ( $i - 1 ) * $SCALE;
		$binlist[$nbin - $i + 1] = $i;
	}
}

if ( ! @binlist )	{
	print "\nExiting the scale '$SCALE' seems to be misformatted\n\n";
	exit;
}

my (@shortnames,%maxma,%minma,$fref,%incoll,$totalcolls,$lastcoll,$totaloccs,$binnedcolls,%bincolls,%binoccs,%binabund,%bybin,%bybinref,%collfromref,%id,%idname,%allseen,%abund,%raw);
for my $file ( @FILES )	{
	my $shortname = $file;
	$shortname =~ s/\.(txt|csv|tab)//;
	push @shortnames , $shortname;
	if ( ! open IN,"<$PATH/$file" )	{
		print "\nExiting because $file couldn't be found\n\n";
		exit;
	}
	$_ = <IN>;
	s/\n//;
	my @f = split /,/,$_;
	# coll no is column 0
	my ($fgen,$fsp,$fabund,$fbin,$fbin2,$fmax,$fmin,$fextant);
	for my $i ( 0..$#f )	{
		if ( $f[$i] =~ /occurrence(|s).genus_name/ )	{ #genus_name
			$fgen = $i;
		} elsif ( $f[$i] =~ /occurrence(|s).species_reso/ )	{ #species_reso
			$fspreso = $i;
		} elsif ( $f[$i] =~ /occurrence(|s).species_name/ )	{ #species_name
			$fsp = $i;
		} elsif ( $f[$i] =~ /occurrence(|s).abund_value/ )	{ #abund_value
			$fabund = $i;
		} elsif ( $f[$i] eq "stage" && ( $SCALE =~ /stages/i || $SCALE > 0 ) && $SCALE !~ /fr2|peters/i )	{ #
			$fbin = $i;
		} elsif ( $f[$i] eq "subepoch" && $SCALE =~ /epochs/i )	{
			$fbin = $i;
		} elsif ( $f[$i] eq "epoch" && $SCALE =~ /epochs/i )	{
			$fbin2 = $i;
		} elsif ( $f[$i] eq "FR2_bin" && $SCALE =~ /fr2/i )	{
			$fbin = $i;
		} elsif ( $f[$i] eq "Peters.interval" && $SCALE =~ /peters/i )	{
			$fbin = $i;
		} elsif ( $f[$i] eq "X10_my_bin" && $SCALE !~ /stages|epochs|fr2|peters/i && $SCALE =~ /[^0-9]/ )	{
			$fbin = $i;
		} elsif ( $f[$i] eq "interpolated_base" && $INTERPOLATE =~ /y/i )	{
			$fmax = $i;
		} elsif ( $f[$i] eq "interpolated_top" && $INTERPOLATE =~ /y/i )	{
			$fmin = $i;
		} elsif ( $f[$i] =~ /(interval_base|ma_max|max_ma)$/ && $INTERPOLATE !~ /y/i  )	{
			$fmax = $i;
		} elsif ( $f[$i] =~ /(interval_top|ma_min|min_ma)$/ && $INTERPOLATE !~ /y/i )	{
			$fmin = $i;
		} elsif ( $f[$i] =~ /collection(|s).reference_no/ )	{
			$fref = $i;
		} elsif ( $f[$i] eq $RANK."_extant" )	{
			$fextant = $i;
		}
	}

	if ( $RANK =~ /sp/i && ! $fsp )	{
		print "\nExiting because $file doesn't have a species name field\n\n";
		exit;
	} elsif ( ( $SCALE > 0 || $SCALE eq "custom" ) && ! $fmax )	{
		print "\nExiting because $file doesn't have a maximum estimate field\n\n";
		exit;
	} elsif ( ( $SCALE > 0 || $SCALE eq "custom" ) && ! $fmin )	{
		print "\nExiting because $file doesn't have a minimum age estimate field\n\n";
		exit;
	} elsif ( $SCALE == 0 && $SCALE ne "custom" && ! $fbin )	{
		print "\nExiting because $file doesn't have the right time scale field\n\n";
		exit;
	} elsif ( $fref == 0 && ( $COLLSPERREF > 0 || $REFQUOTA > 0 || $SEQUENTIAL =~ /y/i ) )	{
		print "\nExiting because $file doesn't have a collection reference number field\n\n";
		exit;
	}

	while (<IN>)	{
		s/\n//;
		# a comma followed by a space is assumed to be embedded junk
		$_ =~ s/, / /g;
		my @f = split /,/,$_;
		$f[$fbin] =~ s/"//g;
		my $within = $lookup{$f[$fbin]};
		if ( ! $within && $lookup{$f[$fbin2]} )	{
			$within = $lookup{$f[$fbin2]};
		}
		if ( $SCALE > 0 || $SCALE eq "custom" )	{
			my $range = $f[$fmax] - $f[$fmin];
			$within = "";
			for my $b ( 1..$nbin )	{
				if ( $f[$fmax] eq "" )	{
					last;
				}
				my ($upper,$lower) = ($f[$fmin],$f[$fmax]);
				if ( $top[$b] > $f[$fmin] )	{
					$upper = $top[$b];
				}
				if ( $base[$b] < $f[$fmax] )	{
					$lower = $base[$b];
				}
				if ( $top[$b] <= $f[$fmin] && $base[$b] >= $f[$fmax] )	{
					$within = $b;
					last;
				} elsif ( $DEORPHAN >= 0.5 && $lower - $upper >= $range * $DEORPHAN )	{
					$within = $b;
					last;
				}
			}
			if ( $within ne "" )	{
				if ( $f[$fbin] > 0 && $fbin > 0 )	{
					$stageIn{$within}{$f[$fbin]}++;
				}
			}
		} elsif ( $within && $f[$fmax] > 0 )	{
			if ( $f[$fmax] > $maxma{$within} )	{
				$maxma{$within} = $f[$fmax];
			}
			if ( $f[$fmin] > 0 && ( $f[$fmin] < $minma{$within} || $minma{$within} == 0 ) )	{
				$minma{$within} = $f[$fmin];
			}
		}
		$totaloccs++;
		my $thiscoll = $f[0];
		if ( $thiscoll != $lastcoll )	{
			$totalcolls++;
		}
		# forces by-occurrence sampling
		if ( $BYCOLLECTION =~ /n/i )	{
			$f[0] = $totaloccs;
		}
		if ( $within > 0 && ( $ABUND !~ /y/i || $f[$abund] > 0 ) )	{
			my $name = $f[$fgen];
			if ( $RANK =~ /^s/i )	{
				# knocks out sp. spp. indet. most informals
				if ( $f[$fsp] =~ /[^a-z]/ || $f[$fspreso] =~ /inf/ )	{
					next;
				}
				$name .= " ".$f[$fsp];
			}
			if ( $thiscoll != $lastcoll )	{
				$bincolls{$within}++;
				$binnedcolls++;
			}
			$incoll{$f[0]}++;
			if ( $incoll{$f[0]} == 1 )	{
				push @{$bybin{$within}} , $f[0];
				push @{$bybinref{$within}{$f[$fref]}} , $f[0];
				$collfromref{$f[0]} = $f[$fref];
			}
			if ( ! $id{$name} )	{
				$ngen++;
				$id{$name} = $ngen;
				$idname{$ngen} = $name;
				$infile{$ngen} = $shortname;
			}
			push @{$list{$f[0]}} , $id{$name};
			$binoccs{$within}++;
			if ( $ABUND =~ /y/i )	{
				$abund{$f[0]}{$id{$name}} = $f[$fabund];
				$binabund{$within} += $f[$fabund];
			}
			if ( ! $allseen{$name}{$within} )	{
				$raw{$within}++;
				$allseen{$name}{$within}++;
				if ( ! $first{$name} || $within < $first{$name} )	{
					$first{$name} = $within;
				}
				if ( ! $last{$name} || $within > $last{$name} )	{
					$last{$name} = $within;
				}
				$extant{$name} = $f[$fextant];
			}
		}
		$lastcoll = $thiscoll;
	}
	close IN;
}

print "\n";
if ( ! $fref )	{
	print "Warning: because the file doesn't include reference numbers, Good's u will be\n based on taxa found in one collection instead of one reference\n\n";
}

print "full data set: $totalcolls collections and $totaloccs occurrences\n";
my $temp = 0;
$temp += $binoccs{$_} foreach keys %binoccs;
printf "binned data set: $binnedcolls collections and %d occurrences\n",$temp;

my (%refs,%refseen);
for my $bin ( reverse 1..$nbin )	{
	my @byrefs = keys %{$bybinref{$bin}};
	$refs{$bin} = $#byrefs + 1;
	$refseen{$_}++ foreach keys %{$bybinref{$bin}};
}

if ( $fref > 0 )	{
	my @temp = keys %refseen;
	printf "binned data derive from %d references\n",$#temp+1;
}

my (%rawrt,%rawbc,%rawsingle,%gap,%notgap,%rawone,%rawtwo,%rawthree,%rawpart);
my $totalextant;
for my $g ( keys %first )	{
	if ( $extant{$g} =~ /y/i )	{
		$last{$g} = $nbin + 1;  # bug fix
		$allseen{$g}{$nbin+1}++;  # new line
		$totalextant++;
	}
	if ( $first{$g} == $last{$g} )	{
		$rawsingle{$first{$g}}++;
	}
	for my $bin ( $first{$g}..$last{$g} )	{
		$rawrt{$bin}++;
		if ( $bin != $first{$g} )	{
			$rawbc{$bin}++;
		}
		if ( $bin != $first{$g} && $bin != $last{$g} )	{
			if ( ! $allseen{$g}{$bin} )	{
				$gap{$bin}++;
			} else	{
				$notgap{$bin}++;
			}
		}
		if ( ! $allseen{$g}{$bin-1} && $allseen{$g}{$bin} && ! $allseen{$g}{$bin+1} )	{
			$rawone{$bin}++;
		} elsif ( $allseen{$g}{$bin-1} && $allseen{$g}{$bin} )	{
			$rawtwo{$bin}++;
		}
		if ( $allseen{$g}{$bin-1} && $allseen{$g}{$bin} && $allseen{$g}{$bin+1} )	{
			$rawthree{$bin}++;
		} elsif ( $allseen{$g}{$bin-1} && ! $allseen{$g}{$bin} && $allseen{$g}{$bin+1} )	{
			$rawpart{$bin}++;
		}
		# forward gap fillers
		if ( $allseen{$g}{$bin-1} && $allseen{$g}{$bin} && ! $allseen{$g}{$bin+1} && $allseen{$g}{$bin+2} )	{
			$rawXXOX{$bin}++;
		} elsif ( $allseen{$g}{$bin-1} && ! $allseen{$g}{$bin} && ! $allseen{$g}{$bin+1} && $allseen{$g}{$bin+2} )	{
			$rawXOOX{$bin}++;
		}
		# backward gap fillers
		if ( $allseen{$g}{$bin-2} && ! $allseen{$g}{$bin-1} && $allseen{$g}{$bin} && $allseen{$g}{$bin+1} )	{
			$rawXXOXb{$bin}++;
		} elsif ( $allseen{$g}{$bin-2} && ! $allseen{$g}{$bin-1} && ! $allseen{$g}{$bin} && $allseen{$g}{$bin+1} )	{
			$rawXOOXb{$bin}++;
		}
	}
}
$ranks = ( $RANK =~ /^s/i ) ? "species" : "genera";
if ( $totalextant > 0 )	{
	print "$totalextant of the $ngen binned $ranks are extant\n\n";
} else	{
	print "none of the $ngen binned $ranks are extant\n\n";
}

my (%everdrawn,%everseen,%failed,%lastfailed,%collsused,%collsdrawn,%occsdrawn,%singledrawn,%subone,%subtwo,%maxgdrawn,%freq1,%freq2,%freq3,%freqmax,%ccurve,%atq,%bpdbybin);
$|=1;
print "trials: ";
for my $t ( 1..$TRIALS )	{
	print "\rtrials: $t ";
	my (%seen,%lastused,%usedintrial,%trialcolls,%genera);
	for my $bin ( reverse 1..$nbin )	{
		if ( $bincolls{$bin} >= $MINCOLLBYBIN && ( ! $failed{$bin} || $USEFAILED =~ /y/i ) )	{
			my (@colls,%inref,%drawn,%drawnfromref);
			my %collsinref = ();
			for my $r ( keys %{$bybinref{$bin}} )	{
				my $n = $#{$bybinref{$bin}{$r}} + 1;
				$collsinref{$_} = $n foreach @{$bybinref{$bin}{$r}};
			}
			if ( ( $REFQUOTA > 0 || $SEQUENTIAL =~ /y/i ) && $refs{$bin} >= $REFQUOTA )	{
				my @byrefs = keys %{$bybinref{$bin}};
				my @subrefs = ();
				my $myquota;
				if ( $REFQUOTA > 0 )	{
					$myquota = $REFQUOTA;
				} else	{
					$myquota = $refs{$bin};
				}
				while ( $#subrefs + 1 < $myquota && $#byrefs > -1 )	{
					my $i = int( rand( $#byrefs + 1 ) );
					push @subrefs , $byrefs[$i];
					splice @byrefs , $i , 1;
				}
				for my $s ( @subrefs )	{
					my @temp = @{$bybinref{$bin}{$s}};
					if ( $COLLSPERREF > 0 && scalar(@temp) >= $COLLSPERREF )	{
						for my $i ( 1..$COLLSPERREF )	{
							my $x = int( rand( $#temp + 1 ) );
							push @colls , $temp[$x];
							splice @temp , $x , 1;
						}
					} else	{
						push @colls , @temp;
					}
					for my $c ( @{$bybinref{$bin}{$s}} )	{
						$inref{$_}{$s}++ foreach @{$list{$c}};
					}
				}
			} elsif ( $COLLSPERREF > 0 )	{
				for my $r ( keys %{$bybinref{$bin}} )	{
					my @temp = @{$bybinref{$bin}{$r}};
					if ( $COLLSPERREF > 0 && scalar(@temp) >= $COLLSPERREF )	{
						for my $i ( 1..$COLLSPERREF )	{
							my $x = int( rand( $#temp + 1 ) );
							push @colls , $temp[$x];
							splice @temp , $x , 1;
						}
					} else	{
						push @colls , @temp;
					}
				}
			} else	{
				@colls = @{$bybin{$bin}};
				for my $r ( keys %{$bybinref{$bin}} )	{
					for my $c ( @{$bybinref{$bin}{$r}} )	{
						$inref{$_}{$r}++ foreach @{$list{$c}};
					}
				}
			}

			# these variables need to be recomputed during
			#  every trial because of the REFQUOTA algorithm
			@{$trialcolls{$bin}} = @colls;
			my $nocc = 0;
			my $ncoll = $#colls + 1;
			my %freq = ();

			for my $c ( @colls )	{
				if ( $ABUND =~ /y/i )	{
					$nocc += $abund{$c}{$_} foreach @{$list{$c}};
					$freq{$_} += $abund{$c}{$_} foreach @{$list{$c}};
				} else	{
					$nocc += $#{$list{$c}} + 1;
					$freq{$_}++ foreach @{$list{$c}};
				}
			}
			my ($maxocc,$maxc,$maxg) = (0,0,"");
			for my $g ( keys %freq )	{
				if ( $freq{$g} == 1 )	{
					$freq1{$bin}[$t]++;
				} elsif ( $freq{$g} == 2 )	{
					$freq2{$bin}[$t]++;
				} elsif ( $freq{$g} == 3 )	{
					$freq3{$bin}[$t]++;
				}
			}

			# classical rarefaction of occurrences
			if ( $METHOD eq "CR" )	{
				my @occs = ();
				for my $c ( @colls )	{
					if ( $ABUND =~ /y/i )	{
						for my $g ( @{$list{$c}} )	{
							for my $i ( 1..int($abund{$c}{$g}) )	{
								push @occs , $g;
							}
						}
					} else	{
						push @occs , $_ foreach @{$list{$c}};
					}
				}
				if ( $#occs + 1 >= $QUOTA )	{
					for $q ( 1..$QUOTA )	{
						my $i = int( rand( $#occs + 1 ) );
						my $g = $occs[$i];
						splice @occs , $i , 1;
						$drawn{$g}++;
						$everdrawn{$idname{$g}}{$bin}++;
						$occsdrawn{$bin}[$t]++;
						my @temp = keys %drawn;
						$ccurve{$bin}{$q} += $#temp + 1;
						$atq{$bin}{$q}++;
					}
				} else	{
					$failed{$bin}++;
					$lastfailed{$bin} = $t;
					%drawn = ();
				}
				push @{$genera{$bin}} , $_ foreach keys %drawn;
			# unweighted list subsampling
			} elsif ( $METHOD eq "UW" )	{
				for $q ( 1..$QUOTA )	{
					if ( $#colls == -1 )	{
						$failed{$bin}++;
						$lastfailed{$bin} = $t;
						%drawn = ();
						last;
					}
					my $c = int( rand( $#colls + 1 ) );
					$drawn{$_}++ foreach @{$list{$colls[$c]}};
					$everdrawn{$idname{$_}}{$bin}++ foreach @{$list{$colls[$c]}};
					$collsdrawn{$bin}[$t]++;
					$occsdrawn{$bin}[$t] += $#{$list{$colls[$c]}};
					splice @colls , $c , 1;
					my @temp = keys %drawn;
					$ccurve{$bin}{$q} += $#temp + 1;
					$atq{$bin}{$q}++;
				}
				push @{$genera{$bin}} , $_ foreach keys %drawn;
			# occurrences squared-weighted list subsampling
			} elsif ( $METHOD eq "O2W" )	{
				my $sumo2 = 0;
				while ( $sumo2 < $QUOTA )	{
					if ( $#colls == -1 )	{
						$failed{$bin}++;
						$lastfailed{$bin} = $t;
						%drawn = ();
						last;
					}
					my $c = int( rand( $#colls + 1 ) );
					my $o2 = ( $#{$list{$colls[$c]}} + 1 )**2;
					if ( $o2 + $sumo2 > $QUOTA && $o2 + $sumo2 - $QUOTA > $QUOTA - $sumo2 )	{
						last;
					}
					$drawn{$_}++ foreach @{$list{$colls[$c]}};
					$everdrawn{$idname{$_}}{$bin}++ foreach @{$list{$colls[$c]}};
					$collsdrawn{$bin}[$t]++;
					$occsdrawn{$bin}[$t] += $#{$list{$colls[$c]}};
					$sumo2 += $o2;
					splice @colls , $c , 1;
				}
				push @{$genera{$bin}} , $_ foreach keys %drawn;
			# default is SQS
			} else	{
				# SQS requires at least two collections
				if ( $bincolls{$bin} < 2 )	{
					next;
				}
				my ($notDominant,$singletons,$doubletons,$bpd) = (0,0,0,0);
				my ($maxocc,$maxc,$maxg) = (0,0,"");
				for my $c ( @colls )	{
					if ( $#{$list{$c}} + 1 > $maxocc )	{
						$maxocc = $#{$list{$c}} + 1;
						$maxc = $c;
					}
				}
			# treat singletons appearing in the largest collection
			#  as non-singletons
				my %inmax = ();
				if ( $BIGGEST =~ /y/i )	{
					for my $g ( @{$list{$maxc}} )	{
						$inmax{$g}++;
					}
				}
			# count singletons and doubletons and identify dominant
				delete $freq{''};
				my %thisfreq;
				for my $g ( keys %freq )	{
					$notDominant++;
					my @refs = keys %{$inref{$g}};
					$thisfreq{$g} = $#refs + 1;
					if ( $SINGLETONS !~ /ref/i )	{
						$thisfreq{$g} = $freq{$g};
					}
					if ( $thisfreq{$g} > 1 || $inmax{$g} > 0 )	{
						if ( $thisfreq{$g} == 2 )	{
							$doubletons++;
						}
						if ( $freq{$g} / $nocc > $bpd )	{
							$bpd = $freq{$g} / $nocc;
							$maxg = $g;
						}
					} else	{
						$singletons++;
					}
				}
			# include most common (dominant) taxon
				if ( $EXCLUDE !~ /y/i )	{
					$maxg = "";
				} else	{
					$notDominant--;
					if ( $thisfreq{$maxg} == 2 )	{
						$doubletons--;
					}
				}
				$freqmax{$bin} = $freq{$maxg};
				$bpdbybin{$bin} += $bpd;

			# exact SQS algorithm that requires sampling all
			#  occurrences and counting taxa at quorum level(s)
				if ( $EXACT !~ /n/i )	{
					my (@temp,@pool,%freqnow,%firstseen,%hitat,@drawnfrom,$lastfrom,$collsnow,$occsnow,$maxnow,$allnow,$singlenow,@fromnow,@subnow,$u,$lastu,$occsdrawn);
				# if SEQUENTIAL is true the collections are
				#  already in random order (and clustered by
				#  reference)
				# if not, randomize them here
					if ( $SEQUENTIAL !~ /y/i )	{
						my @temp;
						while ( @colls )	{
							my $x = int( rand( $#colls + 1 ) );
							push @temp, ( splice @colls , $x , 1 );
						}
						@colls = @temp;
					}
				# randomize occurrences within collections
				#  one way or another
					for my $c ( @colls )	{
						my @temp = @{$list{$c}};
						while ( @temp )	{
							my $x = int( rand( $#temp + 1 ) );
							push @pool , splice( @temp , $x , 1 );
							push @drawnfrom , $c;
						}
					}
					while ( @pool )	{
						my $g = shift @pool;
						my $from = shift @drawnfrom;
						if ( $from != $lastfrom )	{
							$lastfrom = $from;
							$collsnow++;
						}
						$freqnow{$g}++;
						$occsnow++;
						if ( $g eq $maxg )	{
							$maxnow++;
						}
						push @{$hitat{$g}} , $occsnow;
						if ( $freqnow{$g} == 1 )	{
							$allnow++;
							if ( ( $inmax{$g} == 0 || $freq{$g} > 1 ) && $g ne $maxg )	{
								$singlenow++;
							}
							$firstseen{$g} = $occsnow;
						} elsif ( $freqnow{$g} == 2 && ( $inmax{$g} == 0 || $freq{$g} > 1 ) && $g ne $maxg )	{
							$singlenow--;
						}
						$u = 0;
						if ( $occsnow > $maxnow )	{
							$u = 1 - $singlenow / ( $occsnow - $maxnow );
						}
						if ( $TRIALS <= 100 )	{
							push @{$uyielding[$bin][$allnow]} , $u;
							$didyield[$bin][$allnow]{$t}++;
						}
						if ( $allnow > $maxyielding )	{
							$maxyielding = $allnow;
						}
				# record current diversity
				# the last draw couldn't have been a singleton,
				#  so diversity couldn't have changed
				# record each sample size at which the quorum
				#  is exceeded
						if ( $u >= $QUORUM && $lastu < $QUORUM )	{
							push @fromnow , $collsnow;
							push @subnow , $occsnow;
						}
						$lastu = $u;
					}
				# identify taxa that were drawn before or at
				#  the median occurrence count yielding a quorum
					if ( @subnow )	{
						my ($i,$j) = ( int ( $#subnow / 2 ),int( ( $#subnow + 1 ) / 2 ) );
						$occsdrawn{$bin}[$t] = ( $subnow[$i] + $subnow[$j] ) / 2;
						$collsdrawn{$bin}[$t] = ( $fromnow[$i] + $fromnow[$j] ) / 2;
						my $o = int( $occsdrawn{$bin}[$t] );
						for my $g ( keys %firstseen )	{
							if ( $firstseen{$g} <= $o )	{
				# compute the number of times each taxon was
				#  drawn in order to identify the dominant in
				#  the subsample 
								my @temp = @{$hitat{$g}};
								while ( $temp[$#temp] > $o )	{
									pop @temp;
								}
								$drawn{$g} = $#temp + 1;
								$everdrawn{$idname{$g}}{$bin} += $#temp + 1;
							}
						}
					} else	{
						$failed{$bin}++;
						$lastfailed{$bin} = $t;
						$maxbybin{$bin} -= $maxsum;
						$bpdbybin{$bin} -= $bpd;
						if ( $USEFAILED =~ /y/i )	{
							%drawn = %freqnow;
						} else	{
							$occsdrawn{$bin}[$t] = -1;
						}
					}
				}

			# inexact SQS algorithm employing estimated frequencies
			# NOTE: Good's u sensu stricto isn't used here at all
			#  because the correction factor is subtracted from the
			#  frequencies instead of multiplied by them, and the
			#  equation's numerator and denominator are both
			#  different; the former involves the doubleton counts
			#  and the latter is the number of genera and not of
			#  occurrences
				else	{
					my $correction;
					if ( $notDominant > 0 )	{
						$correction = ( $singletons + $doubletons / 2 ) / $notDominant;
					} else	{
						$failed{$bin}++;
						$lastfailed{$bin} = $t;
					}
					my %prop = ();
					my $maxsum = 0;
					if ( $nocc > $freq{$maxg} )	{
						$prop{$_} = ( $freq{$_} - $correction ) / ( $nocc - $freq{$maxg} ) foreach keys %freq;
						delete $prop{$maxg};
						$maxsum += $prop{$_} foreach keys %prop;
					}
					$maxbybin{$bin} += $maxsum;

					my ($sumprop,$stopped) = (0,0);
					while ( $sumprop < $QUORUM )	{
						if ( $#colls == -1 || $nocc - $freq{$maxg} < 1 )	{
							$failed{$bin}++;
							$lastfailed{$bin} = $t;
							$maxbybin{$bin} -= $maxsum;
							$bpdbybin{$bin} -= $bpd;
							last;
						}
						my $c = -1;
				# draw collections with a probability inversely
				#  proportional to the number of collections
				#  per reference
						if ( $DISPERSE =~ /y/i && $SEQUENTIAL !~ /y/i )	{
							while ( $c == -1 )	{
								$c = int( rand( $#colls + 1 ) );
								if ( rand() > 1 / $collsinref{$colls[$c]} )	{
									$c = -1;
								}
							}
						}
				# draw all collections in a given reference
				#  before going to the next one
				# the collections were grouped by reference
				#  earlier in the reference randomization
				#  section, so simply draw them in order
						elsif ( $SEQUENTIAL =~ /y/i )	{
							$c = 0;
						}
				# conventional algorithm ignoring references
						else	{
							$c = int( rand( $#colls + 1 ) );
						}
						my $fromcoll = 0;
						for my $g ( @{$list{$colls[$c]}} )	{
							if ( ! $drawn{$g} && $g ne $maxg )	{
								$sumprop += $prop{$g};
								if ( $sumprop >= $QUORUM )	{
									$stopped++;
								}
							}
							$drawn{$g}++;
							$everdrawn{$idname{$g}}{$bin}++;
							$drawnfromref{$g}{$collfromref{$colls[$c]}}++;
					# random throwback algorithm
							if ( rand() > $prop{$g} && $sumprop >= $QUORUM )	{
								$sumprop -= $prop{$g};
								$everdrawn{$idname{$g}}{$bin} -= $drawn{$g};
								delete $drawn{$g};
								delete $drawnfromref{$g};
							} else	{
								if ( $ABUND =~ /y/i )	{
									$occsdrawn{$bin}[$t] += $abund{$colls[$c]}{$g};
								} else	{
									$occsdrawn{$bin}[$t]++;
								}
								$fromcoll++;
							}
							if ( $stopped > 0 )	{
								last;
							}
						}
						if ( $fromcoll > 0 )	{
							$collsdrawn{$bin}[$t]++;
						}
						if ( $stopped > 0 )	{
							last;
						}
						splice @colls , $c , 1;
					}
				}
				for my $g ( keys %drawn )	{
					push @{$genera{$bin}} , $g;
					my @refs = keys %{$drawnfromref{$g}};
					my $thisfreq = $#refs + 1;
					if ( $SINGLETONS !~ /ref/i )	{
						$thisfreq = $drawn{$g};
					}
					if ( $thisfreq == 1 && $g ne $maxg && ( $inmax{$g} == 0 || $freq{$g} > 1 ) )	{
						$singledrawn{$bin}[$t] += $drawn{$g};
					}
					if ( $drawn{$g} == 1 )	{
						$subone{$bin}[$t]++;
					} elsif ( $drawn{$g} == 2 )	{
						$subtwo{$bin}[$t]++;
					}
				}
				$maxgdrawn{$bin}[$t] = $drawn{$maxg};
			}
		}
	}

	# compute tallies for this trial
	my (%first,%last);
	for my $g ( keys %extant )	{
		if ( $extant{$g} =~ /y/i )	{
			$last{$id{$g}} = $nbin + 1;
			$seen{$nbin+1}{$id{$g}}++;
		}
	}
	for my $bin ( reverse 1..$nbin )	{
		if ( $bincolls{$bin} >= $MINCOLLBYBIN && ( $lastfailed{$bin} != $t || ! $failed{$bin} || $USEFAILED =~ /y/i ) )	{
			for my $c ( @{$trialcolls{$bin}} )	{
				$collsused{$bin}++;
				for my $g ( @{$list{$c}} )	{
					if ( $lastused{$g} != $bin )	{
						$genused{$bin}++;
						$usedintrial{$bin}++;
						$lastused{$g} = $bin;
					}
				}
			}
			for my $g ( @{$genera{$bin}} )	{
				if ( $g > 0 && ! $seen{$bin}{$g} )	{
					$sib{$bin}[$t]++;
					$byfile{$bin}{$infile{$g}}++;
					$first{$g} = $bin;
					if ( $last{$g} eq "" )	{
						$last{$g} = $bin;
					}
				}
				if ( $g > 0 && ! $seen{$bin}{$g} && $lastfailed{$bin+1} != $t && $lastfailed{$bin+2} != $t )	{
					if ( $seen{$bin+1}{$g} )	{
						$two{$bin+1}[$t]++;
					}
					if ( $seen{$bin+1}{$g} && $seen{$bin+2}{$g} )	{
						$three{$bin+1}[$t]++;
						$sumthree++;
						$binsumthree{$bin+1}++;
					} elsif ( ! $seen{$bin+1}{$g} && $seen{$bin+2}{$g} )	{
						$part{$bin+1}[$t]++;
						$sumpart++;
						$binsumpart{$bin+1}++;
					}
					# forward gap fillers
					if ( $seen{$bin+1}{$g} && ! $seen{$bin+2}{$g} && $seen{$bin+3}{$g} )	{
						$XXOX{$bin+1}[$t]++;
					} elsif ( ! $seen{$bin+1}{$g} && ! $seen{$bin+2}{$g} && $seen{$bin+3}{$g} )	{
						$XOOX{$bin+1}[$t]++;
					}
					# backward gap fillers
					if ( ! $seen{$bin+1}{$g} && $seen{$bin+2}{$g} && $seen{$bin+3}{$g} )	{
						$XXOXb{$bin+2}[$t]++;
					} elsif ( ! $seen{$bin+1}{$g} && ! $seen{$bin+2}{$g} && $seen{$bin+3}{$g} )	{
						$XOOXb{$bin+2}[$t]++;
					}
				}
				$seen{$bin}{$g}++;
				$everseen{$idname{$g}}{$bin}++;
			}
			for my $g ( keys %{$seen{$bin+1}} )	{
				if ( ! $seen{$bin}{$g} && ! $seen{$bin+2}{$g} )	{
					$one{$bin+1}[$t]++;
				}
			}
		} else	{
			$collsdrawn{$bin}[$t] = -1;
			$occsdrawn{$bin}[$t] = -1;
			$freq1{$bin}[$t] = -1;
			$freq2{$bin}[$t] = -1;
			$freq3{$bin}[$t] = -1;
			$sib{$bin}[$t] = -1;
			$one{$bin+1}[$t] = -1;
			$two{$bin+1}[$t] = -1;
			$three{$bin+1}[$t] = -1;
			$sumthree -= $binsumthree{$bin+1};
			$part{$bin+1}[$t] = -1;
			$sumpart -= $binsumpart{$bin+1};
			($XXOX{$bin+1}[$t],$XOOX{$bin+1}[$t],$XXOXb{$bin+1}[$t],$XOOXb{$bin+1}[$t]) = (-1,-1,-1,-1);
			($rt{$bin}[$t],$bc{$bin}[$t],$single{$bin}[$t]) = (-1,-1,-1);
			$singledrawn{$bin}[$t] = -1;
			$subone{$bin}[$t] = -1;
			$subtwo{$bin}[$t] = -1;
			$maxgdrawn{$bin}[$t] = -1;
		}
		if ( $lastfailed{$bin} == $t )	{
			my %freq;
			for my $c ( @{$bybin{$bin}} )	{
				$freq{$_}++ foreach @{$list{$c}};
			}
		}
	}
	for my $g ( keys %first )	{
		for my $bin ( $first{$g}..$last{$g} )	{
			$rt{$bin}[$t]++;
		}
		if ( $first{$g} == $last{$g} )	{
			$single{$first{$g}}[$t]++;
		} else	{
			for my $bin ( $first{$g}+1..$last{$g} )	{
				$bc{$bin}[$t]++;
			}
		}
	}
}
print "\n";

my ($firstbin,$lastbin) = ($nbin,1);
if ( $TRIM =~ /y/i )	{
	($firstbin,$lastbin) = ("","");
	for $b ( reverse 1..$nbin )	{
		if ( $bincolls{$b} >= $MINCOLLBYBIN && ( $failed{$b} < $TRIALS || $USEFAILED =~ /y/i ) )	{
			if ( ! $firstbin )	{
				$firstbin = $lastbin;
			}
			$lastbin = $b;
		}
	}
	if ( ! $firstbin )	{
		$firstbin = $lastbin;
	}
}

if ( $METHOD =~ /CR|UW/ )	{
	open OUT,">./sqs/output/quorumSubsampleCurves4.txt";
	print OUT "items";
	for $b ( reverse $lastbin..$firstbin )	{
		print OUT "\t$binlist[$b]";
	}
	print OUT "\n";
	for my $q ( 1..$QUOTA )	{
		my ($step,$step2) = (0,1);
		if ( $q > 1 )	{
			$step = int( log( $q - 1 ) / log( 10 ) * 10 );
			$step2 = int( log( $q ) / log( 10 ) * 10 );
		}
		if ( $step < $step2 || $q == $QUOTA )	{
			print OUT "$q";
			for $b ( reverse $lastbin..$firstbin )	{
				if ( $atq{$b}{$q} > 0 )	{
					printf OUT "\t%.1f",$ccurve{$b}{$q}/$atq{$b}{$q};
				} else	{
					print OUT "\tNA";
				}
			}
			print OUT "\n";
		}
	}
	close OUT;
} elsif ( $METHOD eq "SQS" && $EXACT !~ /n/i && $TRIALS <= 100 )	{
	# print median u values for diversity levels that were reached at
	#  least 99% of the time
	open OUT,">./sqs/output/quorumSubsampleCurves4.txt";
	print OUT "diversity";
	for $b ( reverse $lastbin..$firstbin )	{
		print OUT "\t$binlist[$b]";
	}
	print OUT "\n";
	for $d ( 1..$maxyielding )	{
		print OUT $d;
		for $b ( reverse $lastbin..$firstbin )	{
			@us = @{$uyielding[$b][$d]};
			@us = sort { $a <=> $b } @us;
			@did = keys %{$didyield[$b][$d]};
			if ( @us && $#did + 1 >= 0.99 * $TRIALS )	{
				my ($i,$j) = ( int ( $#us / 2 ),int( ( $#us + 1 ) / 2 ) );
				printf OUT "\t%.4f",( $us[$i] + $us[$j] ) / 2;
			} else	{
				print OUT "\tNA";
			}
		}
		print OUT "\n";
	}
	close OUT;
}

my %trials;
for my $bin ( reverse 1..$nbin )	{
	$trials{$bin} = $TRIALS;
	if ( $USEFAILED =~ /y/i )	{
		$trials{$bin} = $TRIALS - $failed{$bin};
	}
}

if ( $MATCHMAX =~ /y/i )	{
	my $maxdrawn;
	for my $bin ( reverse 1..$nbin )	{
		if ( $collsdrawn{$bin} / $trials{$bin} > $maxdrawn )	{
			$maxdrawn = $collsdrawn{$bin} / $trials{$bin};
		}
	}
	for my $bin ( reverse 1..$nbin )	{
		if ( $bincolls{$bin} < $maxdrawn || ( $collsused{$bin} > 0 && $collsused{$bin} / $trials{$bin} < $maxdrawn ) )	{
			$maxbybin{$bin} = "";
			$bpdbybin{$bin} = "";
			$sib{$bin} = "";
			$one{$bin+1} = "";
			$two{$bin+1} = "";
			$three{$bin+1} = "";
			$sumthree -= $binsumthree{$bin+1};
			$part{$bin+1} = "";
			$sumpart -= $binsumpart{$bin+1};
			($XXOX{$bin+1},$XOOX{$bin+1},$XXOXb{$bin+1},$XOOXb{$bin+1}) = ("","","","");
			($rt{$bin+1},$bc{$bin+1},$single{$bin+1}) = ("","","");
			$failed{$bin} = $TRIALS;
		}
	}
}

open OUT,">./sqs/output/quorumSubsample4.txt";
print OUT "Bin name\t";
if ( $SCALE > 0 )	{
	print OUT "Stages\t";
}
if ( $SCALE ne "custom" && $maxma{$b} > 0 )	{
	print OUT "Base\tMidpoint\tDuration\t";
} elsif ( $SCALE eq "custom" )	{
	print OUT "Base\tMidpoint\t";
}
if ( $fref > 0 )	{
	print OUT "References\t";
}
print OUT "Collections\t";
if ( $REFQUOTA > 0 || $COLLSPERREF > 0 )	{
	print OUT "Collections used\t";
}
if ( $METHOD !~ /CR/ )	{
	print OUT "Collections drawn\t";
}
print OUT "Occurrences\t";
if ( $ABUND =~ /y/i )	{
	print OUT "Specimens\t";
	print OUT "Specimens drawn\t";
} else	{
	print OUT "Occurrences drawn\t";
}
if ( $RANK =~ /^s/i )	{
	print OUT "Raw species\t";
} else	{
	print OUT "Raw genera\t";
}
if ( ( $REFQUOTA > 0 || $COLLSPERREF > 0 ) && $RANK =~ /^s/i )	{
	print OUT "Species used\t";
} elsif ( $REFQUOTA > 0 || $COLLSPERREF > 0 )	{
	print OUT "Genera used\t";
}
print OUT "Raw one timers\tRaw two timers\tRaw three timers\tRaw part timers";
print OUT "\tRaw forward gap fillers\tRaw backward gap fillers";
print OUT "\tRaw range throughs\t\Raw boundary crossers\tRaw single intervals\tGaps";
print OUT "\tGap statistic\tChao-2";
print OUT "\tSubsampled diversity\tThree timer diversity";
if ( $#FILES > 0 )	{
	for my $file ( @shortnames )	{
		print OUT "\t$file";
	}
}
print OUT "\tOne timers\tTwo timers\tThree timers\tPart timers\tThree timer sampling";
print OUT "\tForward gap fillers\tBackward gap fillers";
print OUT "\tRange throughs\tBoundary crossers\tSingle intervals";
print OUT "\tDrawn once";
print OUT "\tSubsampled Chao-2";
print OUT "\tDominance";
print OUT "\t\"Good's u\"";
if ( $METHOD =~ /^(SQS|)$/ )	{
	print OUT "\tMaximum quorum";
	print OUT "\tQuorum drawn";
}
if ( $USEFAILED =~ /y/i )	{
	print OUT "\tFailed trials";
}
print OUT "\n";
my $grandsprob = 1;
if ( $sumthree > 0 )	{
	$grandsprob = $sumthree / ( $sumthree + $sumpart );
	printf "overall average three timer sampling stat: %.3f\n",$grandsprob;
}
my @under;
for $b ( reverse $lastbin..$firstbin )	{
	if ( $SCALE == 0 )	{
		print OUT "\"$binlist[$b]\"";
		if ( $SCALE ne "custom" && $maxma{$b} > 0 )	{
			printf OUT "\t%.2f\t%.2f\t%.2f",$maxma{$b},( $maxma{$b} + $minma{$b} ) / 2,$maxma{$b} - $minma{$b};
		} elsif ( $SCALE eq "custom" )	{
			printf OUT "\t%.2f\t%.2f",$base[$b],( $base[$b] + $top[$b] ) / 2;
		} else	{
			print OUT "\tNA\tNA\tNA";
		}
	} else	{
		printf OUT "bin $b";
		my $stages = ( $stageIn{$b} ) ? '"'.join(', ',sort(keys(%{$stageIn{$b}}))).'"' : "NA";
		printf OUT "\t$stages";
		printf OUT "\t%.2f\t%.2f",$base[$b],( $base[$b] + $top[$b] ) / 2;
	}
	if ( ! $fref )	{
		printf OUT "\t$refs{$b}";
	} elsif ( $refs{$b} > 0 )	{
		printf OUT "\t$refs{$b}";
		print OUT "\t$bincolls{$b}";
	} else	{
		print OUT "\tNA\tNA";
	}
	if ( $notgap{$b} + $gap{$b} > 0 )	{
		$gapstat{$b} = sprintf( "%.3f",$notgap{$b} / ( $notgap{$b} + $gap{$b} ) );
	} else	{
		($gap{$b},$notgap{$b},$gapstat{$b}) = ("NA","NA","NA");
	}
	my $moccsdrawn = mean(\@{$occsdrawn{$b}});
	if ( $bincolls{$b} >= $MINCOLLBYBIN && ( $bincolls{$b} > 1 || $METHOD =~ /CR|UW|O2W/i ) && $moccsdrawn > 0 && $failed{$b} < $TRIALS && ( $failed{$b} == 0 || $USEFAILED =~ /y/i ) )	{
		my $mcollsdrawn;
		if ( $REFQUOTA > 0 || $COLLSPERREF > 0 )	{
			$collsused{$b} /= $trials{$b};
			printf OUT "\t%.2f",$collsused{$b};
		}
		if ( $METHOD !~ /CR/ )	{
			$mcollsdrawn = mean(\@{$collsdrawn{$b}});
			printf OUT "\t%.2f",$mcollsdrawn;
		}
		print OUT "\t$binoccs{$b}";
		if ( $ABUND =~ /y/i )	{
			print OUT "\t$binabund{$b}";
		}
		printf OUT "\t%.2f",$moccsdrawn;
		$genused{$b} /= $trials{$b};
		my $mfreq1 = mean(\@{$freq1{$b}});
		my $mfreq2 = mean(\@{$freq2{$b}});
		my $mfreq3 = mean(\@{$freq3{$b}});
		my $msib = mean(\@{$sib{$b}});
		my $mone = mean(\@{$one{$b}});
		my $mtwo = mean(\@{$two{$b}});
		my $mthree = mean(\@{$three{$b}});
		my $mpart = mean(\@{$part{$b}});
		# raw sampled diversity and genera used
		print OUT "\t$raw{$b}";
		if ( $REFQUOTA > 0 || $COLLSPERREF > 0 )	{
			printf OUT "\t%.2f",$genused{$b};
		}
		printf OUT "\t%d\t%d\t%d\t%d",$rawone{$b},$rawtwo{$b},$rawthree{$b},$rawpart{$b};
		# raw gap fillers
		printf OUT "\t%d",$rawXXOX{$b} + $rawXOOX{$b};
		printf OUT "\t%d",$rawXXOXb{$b} + $rawXOOXb{$b};
		# raw range-through, boundary crossers, single-interval taxa,
		#  and gap taxa
		printf OUT "\t%d\t%d\t%d\t%d",$rawrt{$b},$rawbc{$b},$rawsingle{$b},$gap{$b};
		# raw gap statistic
		print OUT "\t$gapstat{$b}";
		# raw Chao-2 estimator
		if ( $REFQUOTA > 0 || $COLLSPERREF > 0 )	{
			if ( $mfreq2 > 0 )	{
				printf OUT "\t%.2f",$genused{$b} + ( $mfreq1**2 / ( 2 * $mfreq2 ) );
			} else	{
				print OUT "\tNA";
			}
		} else	{
			if ( $mfreq2 > 0 )	{
				printf OUT "\t%.2f",$raw{$b} + ( $mfreq1**2 / ( 2 * $mfreq2 ) );
			} else	{
				print OUT "\tNA";
			}
		}
		# subsampled SIB
		if ( $msib > 0 )	{
			printf OUT "\t%.2f",$msib;
		} else	{
			print OUT "\tNA";
		}
		if ( $mpart > 0 && $msib > 0 )	{
			my $sprob = 1;
			if ( $mthree > 0 && $b != $nbin )	{
				$sprob = $mthree / ( $mthree + $mpart );
			}
		# SIB corrected with sampling stat
			printf OUT "\t%.2f",$msib * $grandsprob / $sprob;
			if ( $#FILES > 0 )	{
				for my $file ( @shortnames )	{
					$byfile{$b}{$file} /= $trials{$b};
					printf OUT "\t%.2f",$byfile{$b}{$file} * $grandsprob / $sprob;
				}
			}
			printf OUT "\t%.2f\t%.2f\t%.2f\t%.2f",$mone,$mtwo,$mthree,$mpart;
			printf OUT "\t%.3f",$sprob;
			printf OUT "\t%.2f",mean(\@{$XXOX{$b}}) + mean(\@{$XOOX{$b}});
			printf OUT "\t%.2f",mean(\@{$XXOXb{$b}}) + mean(\@{$XOOXb{$b}});
			printf OUT "\t%.2f\t%.2f\t%.2f",mean(\@{$rt{$b}}),mean(\@{$bc{$b}}),mean(\@{$single{$b}});
		# subsampled Chao-2 estimator
			printf OUT "\t%.2f",mean(\@{$subone{$b}});
			my ($one,$two) = (mean(\@{$subone{$b}}),mean(\@{$subtwo{$b}}));
			if ( $two > 0 )	{
				printf OUT "\t%.2f",$msib + ( $one**2 / ( 2 * $two ) );
			} else	{
				print OUT "\tNA";
			}

		} else	{
			print OUT "\tNA";
			if ( $#FILES > 0 )	{
				for my $file ( @shortnames )	{
					if ( $msib > 0 )	{
						$byfile{$b}{$file} /= $trials{$b};
						printf OUT "\t%.2f",$byfile{$b}{$file};
					} else	{
						print OUT "\tNA";
					}
				}
			}
			print OUT "\tNA\tNA\tNA\tNA\tNA";
			print OUT "\tNA\tNA\tNA\tNA";
			print OUT "\tNA\tNA\tNA";
		}
		# dominance
		$bpdbybin{$b} /= $trials{$b};
		printf OUT "\t%.4f",$bpdbybin{$b};
		# Good's u (basic equation)
		printf OUT "\t%.3f",1 - $mfreq1 /  ( $binoccs{$b} - $freqmax{$b} );
		# maximum quorum and quorum drawn
		if ( $METHOD =~ /^(SQS|)$/ )	{
			$maxbybin{$b} /= $trials{$b};
			printf OUT "\t%.3f",$maxbybin{$b};
			my $msingle = mean(\@{$singledrawn{$b}});
			my $mmax = mean(\@{$maxgdrawn{$b}});
			if ( $msingle < $moccsdrawn - $mmax )	{
				printf OUT "\t%.3f",1 - $msingle / ( $moccsdrawn - $mmax );
			} else	{
				print OUT "\t0.000";
			}
		}
		if ( $USEFAILED =~ /y/i )	{
			printf OUT "\t%d",$failed{$b};
		}
		print OUT "\n";
	} else	{
		if ( ! $numeric )	{
			push @under , $binlist[$b];
		} else	{
			push @under , $binlist[$b-1];
		}
		# collections used
		if ( $REFQUOTA > 0 || $COLLSPERREF > 0 )	{
			print OUT "\tNA";
		}
		if ( $METHOD !~ /CR/ )	{
			print OUT "\tNA";
		}
		# occurrences
		if ( $binoccs{$b} > 0 )	{
			print OUT "\t$binoccs{$b}";
		} else	{
			print OUT "\tNA";
		}
		if ( $ABUND =~ /y/i && $binabund{$b} > 0 )	{
			print OUT "\t$binabund{$b}";
		} elsif ( $ABUND =~ /y/i )	{
			print OUT "\tNA";
		}
		# occurrences drawn
		print OUT "\tNA";
		# raw sampled diversity and genera used
		if ( $raw{$b} > 0 )	{
			print OUT "\t$raw{$b}";
		} else	{
			print OUT "\tNA";
		}
		if ( $REFQUOTA > 0 || $COLLSPERREF > 0 )	{
			print OUT "\tNA";
		}
		if ( $rawtwo{$b} > 0 )	{
			printf OUT "\t%d\t%d\t%d\t%d",$rawone{$b},$rawtwo{$b},$rawthree{$b},$rawpart{$b};
			printf OUT "\t%d",$rawXXOX{$b} + $rawXOOX{$b};
			printf OUT "\t%d",$rawXXOXb{$b} + $rawXOOXb{$b};
		} else	{
			print OUT "\tNA\tNA\tNA\tNA\tNA\tNA";
		}
		# raw range-through, boundary crossers, single-interval taxa,
		#  and gap taxa
		printf OUT "\t%d\t%d\t%d\t%d",$rawrt{$b},$rawbc{$b},$rawsingle{$b},$gap{$b};
		# raw gap statistic, raw Chao-2, subsampled SIB, and
		#  3T diversity estimate
		print OUT "\tNA\tNA\tNA\tNA";
		if ( $#FILES > 0 )	{
			for my $file ( @shortnames )	{
				print OUT "\tNA";
			}
		}
		# 1T, 2T, 3T, PT, sampling, and gap fillers
		print OUT "\tNA\tNA\tNA\tNA\tNA\tNA\t";
		# RT, BC, single intervals
		print OUT "\tNA\tNA\tNA";
		# drawn once, subsampled Chao-2, dominance
		print OUT "\tNA\tNA\tNA";
		# Good's u, maximum quorum, quorum drawn
		print OUT "\tNA";
		if ( $METHOD =~ /^(SQS|)$/ )	{
			print OUT "\tNA\tNA";
		}
		if ( $USEFAILED =~ /y/i )	{
			print OUT "\tNA";
		}
		print OUT "\n";
	}
}
close OUT;

open OUT,">./sqs/output/quorumSubsampleRanges4.txt";
print OUT "genus";
for $b ( $lastbin..$firstbin )	{
	print OUT "\t$binlist[$b]";
}
print OUT "\n";
my @genera = keys %everseen;
@genera = sort { $a cmp $b } @genera;
for my $g ( @genera )	{
	print OUT $g;
	for $b ( $lastbin..$firstbin )	{
		if ( $everseen{$g}{$b} == 0 )	{
			$n = "0";
		} elsif ( $TRIALS <= 10 )	{
			$n = sprintf "%.2f",$everseen{$g}{$b} / $TRIALS;
		} else	{
			$n = sprintf "%.3f",$everseen{$g}{$b} / $TRIALS;
		}
		print OUT "\t$n";
	}
	print OUT "\n";
}
close OUT;

open OUT,">./sqs/output/quorumSubsampleOccurrences4.txt";
print OUT "genus";
for $b ( $lastbin..$firstbin )	{
	print OUT "\t$binlist[$b]";
}
print OUT "\n";
for my $g ( @genera )	{
	print OUT $g;
	for $b ( $lastbin..$firstbin )	{
		if ( $everdrawn{$g}{$b} == 0 )	{
			$n = "0";
		} elsif ( $TRIALS <= 10 )	{
			$n = sprintf "%.2f",$everdrawn{$g}{$b} / $TRIALS;
		} else	{
			$n = sprintf "%.3f",$everdrawn{$g}{$b} / $TRIALS;
		}
		print OUT "\t$n";
	}
	print OUT "\n";
}
close OUT;

if ( $#under == 0 )	{
	print "$under[0] is under quota\n";
} elsif ( $#under == 1 )	{
	print "$under[0] and $under[1] are under quota\n";
} elsif ( $#under > 0 )	{
	$under[$#under] = "and ".$under[$#under];
	printf "%d intervals are under quota (".join(', ',@under).")\n",$#under + 1;
}

print "\n";

sub mean	{
	my @x = @{$_[0]};
	my $mean;
	for my $i ( 1..$TRIALS )	{
		$mean += $x[$i];
	}
	if ( $mean == 0 )	{
		return 0;
	}
	$mean = $mean / $TRIALS;
	# needed to prevent printf rounding errors
	$mean += 0.0000001;
	return $mean;
}

# version history:
# first draft 19-20.4.09
# version 1.0 release date 20.6.09
# version 2.0 release date 13.9.09
# new features: Pielou's J correction replaced with exclude-most-common-taxon
#  rule; throwback algorithm added; Gradstein stages scale can be used; TRIM
#  option added; outputs more stats, including occurrences drawn and counts of
#  taxa with one, two, or three occurrences; output file column names revised;
#  prints taxa-by-bins output file 
# version 3.0 release date 12.3.10
# new features: coverage now estimated with a version of Good's u that uses
#  counts of taxa found in a single reference instead of taxa with only one
#  occurrence; abundance data and Gradstein subepochs scale can be used;
#  outputs collection and occurrence totals; prints additional totals to screen
# version 3.1 release date 9.7.10
# new features: bug fixes in computation of quorumSubsampleCurves3.txt file
#  and output of two timer counts
# version 3.2 release date 26.7.10 
# new feature: INTERPOLATE option
# version 3.3 release date 11.12.10
# new feature: DISPERSE option
# version 4.0: 19-21.8.11
# new features:
#  improved inexact subsampling algorithm based on sqs version 3.1
#  new EXACT and SEQUENTIAL algorithms
#  BYCOLLECTION, SINGLETONS, EXCLUDE, and BIGGEST made optional
#  means of count variables are geometric instead of arithmetic
#  subsampling curves are printed when EXACT is used
# version 4.1: 3.3.14
# new features:
#  custom time scales can be used
#  multiple updates and bug fixes in time scale and occurrence parsing routines
#  gap filler statistics and lots of additional basic stats are now output
#  switched back to taking arithmetic means because counts can be considered to
#   have been generated by a Poisson process
# version 4.2: 24.6.14
#  a couple of minor bug fixes
# version 4.3: 29.7.14
#  a couple more minor bug fixes

