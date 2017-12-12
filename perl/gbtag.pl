#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

use Getopt::Std;

$Getopt::Std::STANDARD_HELP_VERSION = 1;
my $version = "gbtag.pl v. 0.1\n";
my $help    = "Use as:
  gbtag.pl -g {gb} -p {primary} -t {tags} -H > stdout

 -g - gb file to parse
 -p - (comma-separated list of) primary tag(s) to consider
 -t - comma-separated list of tags to output (includes all explicitely named tags,
      as well as 'start', 'end', 'strand', 'spliced_seq' etc. and a reserved tag 'translate')
 -H - whether to include header in the output
 -d - field delimiter (tab is the default)
";

our($opt_v, $opt_h, $opt_g, $opt_t, $opt_p, $opt_H, $opt_d);

die $help    if !getopts('hvHg:t:p:d:') || $opt_h;
die $version if $opt_v;

sub VERSION_MESSAGE {
	my ($fp) = @_;
	print { $fp } $version;
}

sub HELP_MESSAGE {
	my ($fp) = @_;
	print { $fp } $help;
}

die "Invalid gb specified\n" if not defined $opt_g or not -s $opt_g;
die "No tags specified\n" if not $opt_t;
$opt_d = "\t" if not defined $opt_d or $opt_d eq '\t';
my %primaries;
if (defined $opt_p) {
	$primaries{$_} = 1 for split /,/, $opt_p;
}

my @tags = split /,/, $opt_t;
printf "%s\n", join $opt_d, @tags if defined $opt_H;

my $gb = Bio::SeqIO->new(-file => $opt_g) or die "$!\n";
for my $seq ($gb->next_seq) {
	for my $feature ($seq->get_SeqFeatures) {
		next if $opt_p && not $primaries{$feature->primary_tag};
		my @vals = ();
		push @vals, get_tag_values($feature, $_) for @tags;
		printf "%s\n", join $opt_d, @vals;
	}
}

sub get_tag_values {
	my $feature = shift;
	my $tag = shift;
	if ($feature->can($tag)) {
		my $val = $feature->$tag;
		ref($val) eq 'Bio::PrimarySeq' ? return $val->seq : return $val;
	}
	if ($feature->has_tag($tag)) {
		return join ';', $feature->get_tag_values($tag);
	}
	if ($tag eq 'translate') {
		return translate_feature($feature);
	}
	return '';
}

sub translate_feature {
	my $feature = shift;
	my $seq = $feature->spliced_seq;
	my $table = get_tag_values($feature, 'transl_table');
	$table ? return $seq->translate(-codontable_id => $table)->seq : return '';
}

