#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

use Getopt::Std;

$Getopt::Std::STANDARD_HELP_VERSION = 1;
my $version = "biotags.pl v. 1.0\n";
my $help    = "Use as:
  biotags.pl -i {input} [-p {primary}] -t {tags} [-H] [-d {delimiter}] > stdout

 -i - input file to parse
 -p - comma-separated list of primary tags to consider
 -T - comma-separated list of sequence tags (includes 'primary_id', 'desc' etc.,
      as well as 'translate', 'revcom' etc.)
 -t - comma-separated list of feature tags (includes all explicitely named tags,
      as well as 'start', 'strand', 'spliced_seq', 'get_all_tags' etc. and
      a reserved tag 'join')
 -H - whether to include header in the output
 -d - field delimiter (tab is the default)

from https://github.com/har-wradim/miscngs/
";

our($opt_v, $opt_h, $opt_i, $opt_T, $opt_t, $opt_p, $opt_H, $opt_d);

die $help    if not getopts('hvHi:T:t:p:d:') or $opt_h;
die $version if $opt_v;

sub VERSION_MESSAGE {
	print { shift } $version;
}

sub HELP_MESSAGE {
	print { shift } $help;
}

# check input args
die "No input file specified\n" if not    $opt_i;
print STDERR "No tags specified\n" if not $opt_t and not $opt_T;

$opt_d = "\t" if not defined $opt_d or $opt_d eq '\t';
my %primaries;
if (defined $opt_p) {
	$primaries{$_} = 1 for split /,/, $opt_p;
}

my @stags = ();
my @ftags = ();

@stags = split /,/, lc $opt_T if $opt_T;
@ftags = split /,/, lc $opt_t if $opt_t;

# print the header if requested
printf "%s\n", join $opt_d, (@stags, @ftags) if defined $opt_H;

# open the input and iterate over individual seqs/features
my $in = Bio::SeqIO->new(-file => $opt_i) or die "$!\n";
while (my $seq = $in->next_seq) {
	if ($opt_t) {
		for my $feature ($seq->get_SeqFeatures) {
			next if $opt_p and not $primaries{$feature->primary_tag};
			my @vals = ();
			push @vals, get_seq_values($seq, $_)     for @stags;
			push @vals, get_tag_values($feature, $_) for @ftags;
			printf "%s\n", join $opt_d, @vals if @vals;
		}
	}
	else {
		my @vals = ();
		push @vals, get_seq_values($seq, $_) for @stags;
		printf "%s\n", join $opt_d, @vals if @vals;
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
		my @vals = $feature->$tag;
		return $vals[0]->seq if ref($vals[0]) eq 'Bio::PrimarySeq';
		return join ";", @vals;
	}

	# and reserved tags
	if ($tag eq 'join') {
		return join_feature($feature);
	}
	return '';
}

# get value for a given seq/tag pair
sub get_seq_values {
	my $seq = shift;
	my $tag = shift;
	my %replace = (
		locus      => 'display_id',
		definition => 'desc',
		accession  => 'accession_number',
		organism   => 'species',
		source     => 'species'
	);
	$tag = $replace{$tag} if defined $replace{$tag};
	if ($seq->can($tag)) {
		my @vals = $seq->$tag;
		my $ref = ref($vals[0]);
		return $vals[0]->seq      if $ref eq 'Bio::Seq::RichSeq';
		return $vals[0]->binomial if $ref eq 'Bio::Species';
		return join ";", @vals;
	}
	return '';
}

# extract join'ed segment coordinates
sub join_feature {
	my $feature = shift;
	my $loc = $feature->location;
	return location_range($loc) if not $loc->isa('Bio::Location::SplitLocationI');
	my @f = ();
	push @f, location_range($_) for ($loc->sub_Location);
	return join ",", @f;
}

sub location_range {
	my $loc = shift;
	return sprintf "%s%s..%s%s", pos_type($loc->start_pos_type), $loc->start, pos_type($loc->end_pos_type), $loc->end;
}

sub pos_type {
	my $pos = shift;
	my %vals = ( BEFORE => '<', AFTER => '>', EXACT => '' ); # WITHIN => '', BETWEEN => ''
	return $vals{$pos} if defined $vals{$pos};
	return '';
}
