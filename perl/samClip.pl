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

my @cols;
my $join = 0;
while (my $line = <>) {
	if (substr($line, 0, 1) ne '@') {
		@cols = split(/\t/, $line, 12);
		if ($cols[5] =~ /(\d+)S$/) {
			$cols[9] = substr($cols[9], 0, -$1);
			$cols[10] = substr($cols[10], 0, -$1);
			$cols[5] =~ s/S$/H/;
			$join = 1;
		}
		if ($cols[5] =~ /^(\d+)S/) {
			$cols[9] = substr($cols[9], $1);
			$cols[10] = substr($cols[10], $1);
			$cols[5] =~ s/^(\d+)S/$1H/;
			$join = 1;
		}
		if ($join) {
			print join("\t", @cols);
			$join = 0;
			next;
		}
	}
	print $line;
}
print STDERR "Done\n";
