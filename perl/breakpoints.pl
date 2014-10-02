#!/bin/perl

### Written by Andrey Rozenberg (jaera at yandex.com), Ruhr-Universität Bochum
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.
### You should have received a copy of the GNU General Public License
### along with this program. If not, see <http://www.gnu.org/licenses/>.

## Use as follows: breakpoints.pl -l [linear genomes] -c [circular genomes]

use strict;
use warnings;

my @fnames;
my %chroms;
my %presence;
my %genes;
my %strands;
my %ends;

my $linear = 1;

sub strand {
	my $strand = shift;
	return 1 if !$strand;
	my $strand1 = substr($strand, 0, 1);
	return  1 if ($strand1 eq '+' || $strand1 eq 'f');
	return -1 if ($strand1 eq '-' || $strand1 eq 'r');
	die("Unrecognized strand: '$strand'\n");
}

foreach my $fname (@ARGV) {
	if ($fname eq '-c') {
		$linear = 0;
		next;
	}
	if ($fname eq '-l') {
		$linear = 1;
		next;
	}
	open(FP, "<", $fname) or die("Couldn't open '$fname' for reading\n");
	my $g1;
	my $g2;
	my $chrom1 = "";
	my $chrom2;
	my $strand;
	while (<FP>) {
		($chrom2, $g2, $strand) = split;
		++$presence{$g2};
		die("'$g2' already defined in the genome '$fname'\n") if defined($chroms{$fname}{$g2});
		$chroms{$fname}{$g2} = $chrom2;
		$strands{$fname}{$g2} = strand($strand);
		$genes{$fname}{$chrom2}{$g2} = $.;
		if ($linear && $chrom1 ne $chrom2) {
			$ends{$fname}{$chrom1}[1] = $g1 if defined($g1);
			$ends{$fname}{$chrom2}[0] = $g2;
		}
	} continue {
		$chrom1 = $chrom2;
		$g1 = $g2;
	}
	$ends{$fname}{$chrom1}[1] = $g1 if $linear;
	close(FP);
	push(@fnames, $fname);
}

die("No files specified\n") if $#fnames < 0;

while ( my ($gene, $count) = each %presence ) {
	next if $count == $#fnames + 1;
	foreach my $fname (@fnames) {
		next if !defined($chroms{$fname}{$gene});
		my $chrom = $chroms{$fname}{$gene};
		delete $chroms{$fname}{$gene};
		delete $strands{$fname}{$gene};
		delete $genes{$fname}{$chrom}{$gene};
	}
}

undef %presence;

sub end2end {
	my $ename = shift;
	my $l = shift;
	my $r = shift;
	my $sl = shift;
	my $sr = shift;
	foreach my $fname (@fnames) {
		next if $fname eq $ename;
		next if !defined($ends{$fname});
		my $chroml = $chroms{$fname}{$l};
		my $chromr = $chroms{$fname}{$r};
		next if $chroml eq $chromr;
		return 0 if
			( $sl ==  $strands{$fname}{$l} && $ends{$fname}{$chroml}[1] eq $l   ||
			  $sl == -$strands{$fname}{$l} && $ends{$fname}{$chroml}[0] eq $l ) &&
			( $sr ==  $strands{$fname}{$r} && $ends{$fname}{$chromr}[0] eq $r   ||
			  $sr == -$strands{$fname}{$r} && $ends{$fname}{$chromr}[1] eq $r );
	}
	return 1;
}

my %pairs;
my %include_pairs;

sub vergence {
	my $fname = shift;
	my $g1 = shift;
	my $g2 = shift;
	if ($strands{$fname}{$g1} * $strands{$fname}{$g2} > 0) {
		return ($g1, $g2, 1, 1) if $strands{$fname}{$g1} > 0;
		return ($g2, $g1, 1, 1);
	}
	else {
		my $ng1 = ($g1,$g2)[$g1 ge $g2];
		my $ng2 = ($g1,$g2)[$g1 lt $g2];
		return ($ng1, $ng2, 1, -1) if $strands{$fname}{$g1} > 0;
		return ($ng1, $ng2, -1, 1);
	}
}

sub check_pair {
	my $fname = shift;
	my $g1 = shift;
	my $g2 = shift;
	my ($l, $r, $sl, $sr) = vergence($fname, $g1, $g2);
	my $pair = sprintf "%s%s|%s%s", $sl < 0 ? "-" : "", $l, $sr < 0 ? "-" : "", $r;
	$include_pairs{$pair} = end2end($fname, $l, $r, $sl, $sr) if !defined($include_pairs{$pair});
	$pairs{$pair}{$fname} = 1 if $include_pairs{$pair};
}

print "Pair";

foreach my $fname (@fnames) {
	while ( my ($chrom, $gs) = each %{$genes{$fname}}) {
		my @geneset = sort { ${$gs}{$a} <=> ${$gs}{$b} } keys %{$gs};
		my $g1 = $geneset[0];
		foreach my $g2 (@geneset[1..$#geneset]) {
			check_pair($fname, $g1, $g2);
		} continue {
			$g1 = $g2;
		}
		check_pair($fname, $g1, $geneset[0]) if !defined($ends{$fname});
	}
	printf "\t%s", $fname;
}
print "\n";

my %BPs;

foreach my $pair ( keys %pairs ) {
	print $pair;
	for my $f (0..$#fnames) {
		my $fname = $fnames[$f];
		my $state = defined($pairs{$pair}{$fname});
		printf "\t%d", $state;
		for my $f2 ($f + 1..$#fnames) {
			my $fname2 = $fnames[$f2];
			$BPs{"$fname|$fname2"} += ($state != defined($pairs{$pair}{$fname2}));
		}
	}
	print "\n";
}

while ( my ($pair, $BP) = each %BPs) {
	printf STDERR "%s\t%d\n", $pair, $BP;
}
