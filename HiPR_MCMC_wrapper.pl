#!/usr/bin/perl -w
# Wrapper to run HiPR MCMC algorithm
use strict;
use Getopt::Long qw( :config posix_default no_ignore_case );
use Pod::Usage;
use IPC::Open3;
use Symbol qw(gensym);
use IO::File;
use Forks::Super;
$Forks::Super::MAX_PROC = 1;
$Forks::Super::ON_BUSY = 'block';
use File::HomeDir qw(home);
use File::Basename;
use Cwd 'abs_path';
my $SCRIPTDIR = dirname(abs_path(__FILE__));

# flags and variables

my $HOMEDIR = home();
my $hiprDir="$SCRIPTDIR";
my $HIPRmcmc="$hiprDir/HiPR_MCMC";

my $ostype=$^O;
if ( $ostype eq "darwin" )
{
   $HIPRmcmc="$hiprDir/HiPR_MCMC.mac";
}

print "$HIPRmcmc\n"; 

my $readsFile = "";
my $structureFile = "";
my $ratesFile = "";
my $rmin = 40;
my $rmax = 80;
my $n = 100000;
my $numCPU = 1;
my $numStartingPoints = 0;
my $burnIn = 10000;
my $samp = 100;
my $t = "s";
my $model = 1;
my $mcConstant = 999.0;
my $likelihoodOnly = 0;
my $outDir = "HiPR_output";
my $man = 0;
my $help = 0;
my $id = "";

my $doStructure = 0;
my $doParameters = 0;

# Parse options and print usage if there is a syntax error,
# or if usage was explicitly requested.
GetOptions('id=s' => \$id, 'rmin=i' => \$rmin, 'rmax=i' => \$rmax, 'n=i' => \$n, 'numCPU=i' => \$numCPU, 'burn_in=i' => \$burnIn, 'samp=i' => \$samp, 't=s' => \&type_handler, 'm=i' => \$model, 'mc=f' => \$mcConstant, 'outdir=s' => \$outDir, 'likelihood_only' => \$likelihoodOnly, 'help|?' => \$help, 'man' => \$man) or pod2usage(-verbose => 1);

sub type_handler {
	$t = $_[1];
#	die("Invalid type $t") if ($t =~ /[^sd]/);
	die("Invalid type $t\n") unless (($t eq "s") || ($t eq "d") || ($t eq "sd") || ($t eq "ds"));
	if ($t =~ /s/) {
		$doStructure = 1;
	}
	if ($t =~ /d/) {
		$doParameters = 1;
	}
}

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

my $remaining = @ARGV;
if ($remaining < 3) {
	printf STDERR "\nERROR: Missing required arguments(s).\n\n";
	pod2usage(-verbose => 1);
}
$readsFile = shift;
$ratesFile = shift;
$structureFile = shift;
#if ($structureFile =~ /\/([^\/]+)\.rnafold/) {
#if ($structureFile =~ /\/([^\/]+)\.([^.]+)$/) {
#if ($structureFile =~ /([^\/]+)\.([^.]+)$/) {
#	$id = $1;
#}

if (@ARGV) {
	printf STDERR "\nWarning: Extra arguments/options detected. Exiting!\n\n";
	pod2usage(-verbose => 1);
}

$outDir =~ s/\/$//;
if (!(-d "$outDir")) {
	mkdir($outDir);
}
$Forks::Super::MAX_PROC = $numCPU;

# separate starting points
open(STRUCT, "$structureFile") || die "Unable to read from $structureFile: $!\n";
while (my $line = <STRUCT>) {
	$numStartingPoints++;
	chomp($line);
	my $line2 = <STRUCT>; chomp($line2);
        my $startStructFile=sprintf(">%s/%s\.starting\.%.3d\.rnafold", $outDir, $id, $numStartingPoints);
	open(TOUT, $startStructFile) || die "Unable to write $startStructFile: $!\n";
	print TOUT "$line\n$line2\n";
	close(TOUT);
}
close(STRUCT);

# run MCMC and get posterior if all checks passed
print STDERR sprintf("Starting MCMC with %d CPUs ...\n", $numCPU);
my $i = 1;
my @pids;
my $pids = @pids;
for (my $i=1; $i<=$numStartingPoints; $i++) {
	my $starting_fn = sprintf("%s/%s\.starting\.%.3d\.rnafold", $outDir, $id, $i);
	my $out_mcmc_fn = sprintf("%s/%s\.mcmc\.%.3d\.txt", $outDir, $id, $i);
	my $out_err_fn = sprintf("%s/%s\.mcmc\.%.3d\.err", $outDir, $id, $i);
	my $out_log_fn = sprintf("%s/%s\.mcmc\.%.3d\.log", $outDir, $id, $i);
	
	my $likelihoodOnlyFlag = ($likelihoodOnly) ? "-likelihood_only" : "";
	my $cmd = "$HIPRmcmc $starting_fn $readsFile $ratesFile -rmin $rmin -rmax $rmax -r 0 -n $n -t $t -m $model -mc $mcConstant -log $out_log_fn $likelihoodOnlyFlag -out $out_mcmc_fn";
	
	open(ERR, ">$out_err_fn") || die "Unable to write $out_err_fn: $!\n";
	my $out_err = "";
	my $pid = Forks::Super::fork { cmd => $cmd, stderr => \$out_err };
	print STDERR "Running MCMC for sequence $i ...\n";
        print STDERR "$cmd\n";
        print STDERR "outdir=$outDir\n";
        print STDERR "id=$id\n";
        print STDERR "i=$i\n";
	print ERR "$out_err";
	close(ERR);
}
waitall;

 
=head1 NAME
 
HiPR - Bayesian Markov chain Monte Carlo inference of RNA secondary structures.
 
=head1 SYNOPSIS
 
HiPR [options] reads_file rates_file structure_file
 
=head1 ARGUMENTS

=over 8
 
=item B<reads_file>

File containing either DMS-seq or dsRNA/ssNRA-seq reads.

=item B<rates_file>

File containing enzyme digestion rates. The first line should contain four comma-delimited rates corresponding to digestion probabilities 3' of paired [A,C,U,G] nucleotides.
The second line should contain rates 3' of unpaired [A,C,U,G] positions.

=item B<structure_file>

File containing the sequence and structure of the RNA of interest. The provided dot-paren structure is used as the initial estimate.

=back

=head1 OPTIONS

=over 8

=item B<-rmin> INT

Minimum read length. Fragments shorter than this will not be included in likelihood calcuations. Defaults to 40.

=item B<-rmax> INT

Maximum read length. Fragments longer than this will not be included in likelihood calcuations. Defaults to 80.

=item B<-n> INT

Number of MCMC iterations to run. Defaults to 100000.

=item B<-numCPU> INT

Maximum number of CPUs to use. Defaults to 1.

=item B<-burn_in> INT

Number of MCMC iterations to discard as burn-in when computing posterior. Defaults to 10000.

=item B<-samp> INT

Sampling frequency of MCMC iterations after burn-in. Defaults to 100.

=item B<-t> [sd]

Selects which parameter to optimize. B<'s'> indicates structure and B<'d'> indicates digestion rates. Defaults to B<'s'> for structure estimation.

=item B<-m> INT

Select likelihood model to use (1=dsRNA/ssRNA-seq, 2=DMS-seq)

=item B<-mc> DOUBLE

Proportionality constant used in Metropolis-Hastings acceptance criterion (Default = 999.0)

=item B<-likelihood_only>

Only compute the initial log-likelihood (do not perform MCMC iterations)

=item B<-outdir dir>

Directory to output results. Defaults to 'HiPR_output/'.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

Estimate secondary structure and base pairing posteriors for a given RNA sequence based on the distribution of read fragments along the locus.

This program requires a file containing the sequence and initial secondary structure of an RNA of interest, a file containing DMS-seq or dsRNA-seq/ssRNA-seq reads, and a file containing the initial estimates of per-nucleotide modification rates. A Bayesian MCMC algorithm is then used to estimate the base pairing posterior that best fits the observed sequencing reads. The results and intermediate files are written to a directory (HiPR_output/ by default).


The output file posterior.txt contains the base pairing posteriors at each nucleotide position, one entry per line.

The output file structure.txt contains the consensus secondary structure.

=cut

