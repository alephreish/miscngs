#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

use Getopt::Std;

$Getopt::Std::STANDARD_HELP_VERSION = 1;
my $version = "biotags.pl v. 0.1\n";
my $help    = "Use as:
  biotags.pl -i {input} [-p {primary}] -t {tags} [-H] [-d {delimiter}] > stdout

 -i - input file to parse
 -p - comma-separated list of primary tags to consider
 -t - comma-separated list of tags to output (includes all explicitely named tags,
      as well as 'start', 'strand', 'spliced_seq' etc. and reserved tags 'translate'
      and 'join')
 -H - whether to include header in the output
 -d - field delimiter (tab is the default)

from https://github.com/har-wradim/miscngs/
";

our($opt_v, $opt_h, $opt_i, $opt_t, $opt_p, $opt_H, $opt_d);

die $help    if not getopts('hvHi:t:p:d:') or $opt_h;
die $version if $opt_v;

sub VERSION_MESSAGE {
	print { shift } $version;
}

sub HELP_MESSAGE {
	print { shift } $help;
}

# check input args
die "No input file specified\n" if not    $opt_i;
die "Input file not found\n"    if not -s $opt_i;
die "No tags specified\n" if not $opt_t;

$opt_d = "\t" if not defined $opt_d or $opt_d eq '\t';
my %primaries;
if (defined $opt_p) {
	$primaries{$_} = 1 for split /,/, $opt_p;
}
my @tags = split /,/, $opt_t;

# print the header if requested
printf "%s\n", join $opt_d, @tags if defined $opt_H;

# open the input and iterate over individual seqs/features
my $in = Bio::SeqIO->new(-file => $opt_i) or die "$!\n";
for my $seq ($in->next_seq) {
	for my $feature ($seq->get_SeqFeatures) {
		next if $opt_p and not $primaries{$feature->primary_tag};
		my @vals = ();
		push @vals, get_tag_values($feature, $_) for @tags;
		printf "%s\n", join $opt_d, @vals;
	}
}

# get value for a given feature/tag pair
sub get_tag_values {
	my $feature = shift;
	my $tag = shift;

	# first check explicit tags
	if ($feature->has_tag($tag)) {
		return join ';', $feature->get_tag_values($tag);
	}

	# then object members
	if ($feature->can($tag)) {
		my $val = $feature->$tag;
		return $val->seq if ref($val) eq 'Bio::PrimarySeq';
		return $val;
	}

	# and reserved tags
	if ($tag eq 'translate') {
		return translate_feature($feature);
	}
	if ($tag eq 'join') {
		return join_feature($feature);
	}

	return '';
}

# extract join'ed segment coordinates
sub join_feature {
	my $feature = shift;
	my $loc = $feature->location;
	return sprintf "%s..%s", $loc->start, $loc->end if not $loc->isa('Bio::Location::SplitLocationI');
	my @f = ();
	push @f, sprintf "%s..%s", $_->start, $_->end for ($loc->sub_Location);
	return join ",", @f;
}

# translate the feature anew
sub translate_feature {
	my $feature = shift;
	my $seq = $feature->spliced_seq;
	my $table = get_tag_values($feature, 'transl_table');
	return $seq->translate(-codontable_id => $table)->seq if $table;
	return '';
}

