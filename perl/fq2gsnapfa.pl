#!/usr/bin/perl

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
### along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;

die("Specify at least two input fastq-files") if ($#ARGV < 1);

my $reads1;
my $reads2;
my $row1;
my $row2;
my $tmp;
foreach my $file (@ARGV) {
	die("File $file doesn't exist") unless (-e $file);
}

open($reads1, "<", $ARGV[0]) or die("Couldn't open file ".$ARGV[0]);
open($reads2, "<", $ARGV[1]) or die("Couldn't open file ".$ARGV[1]);

my $phase = 1;
my $name;
while (!eof($reads1)) {
	$row1 = readline($reads1);
	$row2 = readline($reads2);
	if ($row1 eq "\n" or $row2 eq "\n") { next; }
	if ($phase == 3) {
		die("Out of phase: expecting '+', found: $row1 and $row2") if ($row1 ne "+\n" or $row2 ne "+\n");
	}
	elsif ($phase == 1) {
		die("Out of phase: expecting \@xxx, found: $row1 and $row2") if (substr($row1, 0, 1) ne '@' or substr($row2, 0, 1) ne '@');
		$name = substr($row1, 1, index($row1, ' ') - 1);
		die("Incongruence: $row1 vs. $row2") if ($name ne substr($row1, 1, index($row2, ' ') - 1));
		print ">$name\n";
	}
	elsif ($phase == 2) {
		print "$row1$row2";
	}
	elsif ($phase == 4) {
		print "+\n$row1$row2";
	}
	die ("The two paired files have different compositions") if ( (eof($reads1) and !eof($reads2)) or (eof($reads2) and !eof($reads1)) );
	$phase = $phase == 4 ? 1 : $phase+1;
}
close($reads1);
close($reads2);
if ($#ARGV > 1) {
	splice(@ARGV, 0, 2);
	foreach my $filename (@ARGV) {
		my $file;
		open($file, "<", $filename) or die("Couldn't open file $filename");
		my $phase = 4;
		while (!eof($file)) {
			$row1 = readline($file);
			if ($row1 eq "\n") { next; }
			if ($phase == 4) {
				$phase = 1;
				die("Out of phase: expecting \@xxx, found: $row1") if (substr($row1, 0, 1) ne '@');
				$name = substr($row1, 1, index($row1, ' ') - 1);
				print ">$name\n";
			}
			else {
				$phase++;
				if ($phase == 3) {
					die("Out of phase: expecting \@xxx, found: $row1") if ($row1 ne "+\n");
					print "+\n";
				}
				else {
					print $row1;
				}
			}
		}
		close($file);
	}
}
print STDERR "Done\n";
