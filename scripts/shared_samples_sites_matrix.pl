#!/usr/bin/env perl
# Perl script for making edit level matrix

use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Cwd;


# Defining variables to hold the arguments
my ($path_edit_files, $output_dir, $tissue, $mincov, $minsamps);

# Parsing command-line options
GetOptions(
    'path_to_edit_files=s'    => \$path_edit_files,
    'output_dir=s'            => \$output_dir,
	'tissue=s'				  => \$tissue,
	'min_coverage=s'          => \$mincov,
	'min_samples=s'           => \$minsamps,
) or die "Error in command line arguments: $!\n";

# Adding "/" to thr end of path if it's not there already
$path_edit_files .= "/" unless $path_edit_files =~ /\/$/;

print "Path to Edit Files: $path_edit_files\n";
print "Output Directory (from args): $output_dir\n";
print "Tissue: $tissue\n";
print "Min Coverage: $mincov\n";
print "Min Samples: $minsamps\n";


my %sitehash;
my %totalhash;
my %lvlhash;


foreach my $file (glob (join "", $path_edit_files, "*.rnaeditlevel.gz")) { # It takes any file ends in "*.rnaeditlevel.gz" in this path
	my $filename = (split/\//,$file)[-1]; # Added by Me
	print "Analyzing: $filename\n"; # Added by Me
	my @names = split(/\./,$filename);
	my $sample = $names[0];
	my $strain = join '-', (split(/-/, $sample))[0], (split(/-/, $sample))[1];

	# open (my $INPUT, "<", $file);
	my $INPUT = IO::Uncompress::Gunzip->new($file) or die "Failed to open $file: $GunzipError\n"; # Added by Me
	# <$INPUT>; # Added by Me - Skip the first line (header) or use line 47 [next if ($chr eq '#chrom');]
	while(<$INPUT>) {
		chomp;
		my @fields = split;
		# my ($chr, $pos, $cov, $edit, $lvl) = ($fields[0], $fields[1], $fields[3], $fields[4],$fields[5]); # It has mistake
		my ($chr, $pos, $cov, $edit, $lvl) = ($fields[0], $fields[1], $fields[2], $fields[3],$fields[4]);

		# Skip header line
		next if ($chr eq '#chrom');

		# Only keep chr1 to chr22
		next unless ($chr =~ /^chr([1-9]|1[0-9]|2[0-2])$/);

		my $site = join ':', $chr, ($pos-1), $pos;
		my $ratio = join '/', $edit, $cov;
		if ($cov >= $mincov) {
			$sitehash{$strain}{$site} = $ratio;
			if ($totalhash{$site}) {
				$totalhash{$site}++;
				$lvlhash{$site} = join ',', $lvlhash{$site}, $ratio;
			} else {
				$totalhash{$site} = 1;
				$lvlhash{$site} = $ratio;
			}
		}
	}
	close $INPUT;
}

# Prepare output file
my $output_file = "${output_dir}${tissue}.edMat.${mincov}cov.${minsamps}samps.noXYM.txt"; # Added by Me
unlink $output_file if -e $output_file; # Added by Me
print "Output file: $output_file\n";
open(my $OUTPUT, ">", $output_file) or die "Could not open file '$output_file': $!\n"; # Added by Me

print $OUTPUT "chrom";
# foreach my $strain (keys %sitehash) {
foreach my $strain (sort keys %sitehash) { # Sorted: Perl hash order is random, so unsorted columns change order every run
	print $OUTPUT " $strain";
}
print $OUTPUT "\n";
# foreach my $site (keys %totalhash) {
foreach my $site (sort keys %totalhash) { # Sorted: stable row order across runs
	if ($totalhash{$site} >= $minsamps) {
		my @lvls = split(/\,/, $lvlhash{$site});
		# @lvls = sort {$b <=> $a} @lvls;
		# @lvls = sort { $b cmp $a } @lvls; # Added by Me

		
			print $OUTPUT "$site";
			# foreach my $strain (keys %sitehash) {
			foreach my $strain (sort keys %sitehash) { # Sorted: must match the header order above
				if ($sitehash{$strain}{$site}) {
					print $OUTPUT " $sitehash{$strain}{$site}";
				} else {
					print $OUTPUT " 0/0";
				}
			}
			print $OUTPUT "\n";
		
	}
}
# system("gzip $output_file") == 0 or die "gzip failed: $?\n";
close $OUTPUT;

### END