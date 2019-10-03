#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use IPC::Cmd qw(can_run);
use Scalar::Util qw(looks_like_number);


my $usage = 
"\nUsage: perl $0 [arguments]
   
   REQUIRED ARGUMENTS
   
   Character matrix input:
   
         --matrix     - file containing character matrix in Nexus format

   Output:
   
         --output     - base name for output files (extensions will be
                        appended)

   PAUP* executable:
   
         --paup       - name (and path if necessary) of PAUP* executable 


   OPTIONAL ARGUMENTS

   Bootstrap:
   
         --bs         - number of bootstrap pseudoreplicates to generate
 
   Retain bootstrap replicates:
   
         --save_reps  - Add this flag to retain quartet tree files and MRP
                        matrices from each bootstrap replicate.
         
   Parallelization:
   
         --forks      - number of forks to run in parallel (default = 1).
                        This option requires the Parallel::ForkManager
                        Perl module.

   Tree search method:
   
         --search     - Specify \"bandb\" to change the final tree search
                        method from the default of a heuristic search with
                        TBR to branch and bound [default: tbr]
   
   EXAMPLE
         perl $0
            --matrix=sample_data/input.nex
            --output=myoutput
            --paup=/usr/local/paup/paup4a166_centos64
            --bs=1000
            --forks=40
            --search=bandb
\n\n";

print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";

our $MATRIX;
our $OUTPUT;
our $BS;
our $FORKS = 1;
our $PAUP;
our $SAVE_REPS;
our $SEARCH = "tbr";

GetOptions(
    'matrix=s'  => \$MATRIX,
    'output=s'  => \$OUTPUT,
    'bs=i'  => \$BS,
    'forks=i'  => \$FORKS,
    'paup=s'  => \$PAUP,
    'search=s'  => \$SEARCH,
    'save_reps'  => \$SAVE_REPS
);

$MATRIX or die ("$usage\nERROR: Must specify input character matrix with --matrix\n\n");
$OUTPUT or die ("$usage\nERROR: Must specify output file name with --output\n\n");
$PAUP or die ("$usage\nERROR: Must specify PAUP* executable name with --paup\n\n");
can_run ($PAUP) or die ("\n$usage\n\nERROR: Could not run PAUP* executable ($PAUP). Provide executable name (and path if necessary) with --paup.\n\n");

unless ($SEARCH eq "tbr" or $SEARCH eq "bandb"){
	die ("\n$usage\n\nERROR: The only acceptable options for the --search paramter are tbr and bandb.\n\n");
}

looks_like_number($FORKS) or die ("$usage\nERROR: Value provided with --forks must be numerical.\n\n");

if ($BS){
	looks_like_number($BS) or die ("$usage\nERROR: Value provided with --bs must be numerical.\n\n");
}


unless ($FORKS == 1){
	require "Parallel/ForkManager.pm";
}


my @nexus_lines = file_to_array ($MATRIX);

my @matrix_lines;
my $in_matrix=0;
my $out_of_matrix;
my $ntax;
my $nchar;
foreach (@nexus_lines){
	chomp $_;
	if ($out_of_matrix){
		if ($_ =~ /^\s*matrix\s*$/i){
			die ("\nERROR: Identified more than one MATRIX sections in $MATRIX\. This script assumes only a single character matrix per input file.\n\n")
		}
	}
	
	if ($in_matrix){
		$_ =~ /^[\s\;]*$/ and next;
		$_ =~ /^\s*[\#\[]/ and next; #skip lines that start with # or [, under the assumption that they are comments 
		if ($_ =~ /end\;/i){
			$in_matrix = 0;
			$out_of_matrix = 1;
		}else{
			push (@matrix_lines, $_);
		}
	}
	
	unless ($out_of_matrix){
		if ($_ =~ /dimensions/i){
			if ($_ =~ /ntax\=(\d+)/i){
				$ntax = $1;
			}
			if ($_ =~ /nchar\=(\d+)/i){
				$nchar = $1;
			}
		}elsif($_ =~ /^\s*matrix\s*$/i){	
			$in_matrix = 1;
		}
	}
}

@matrix_lines or die ("\nERROR: failed to extract character matrix from $MATRIX\. Check that formating is consistent with sample file.\n\n");
if ($ntax){
	my $line_count = scalar (@matrix_lines);
	$line_count == $ntax or print STDERR ("\nWARNING: Identified $line_count lines in character matrix but stated ntax is $ntax\n\n");
}


my %char_hash;

my $char_count;
my @taxa_list;

foreach (@matrix_lines){
	my $taxon;
	if ($_ =~ /^([\w\-]+)\s/){
		$taxon = $1;
		push (@taxa_list, $taxon);
	}else{
		die ("\nERROR: Could not parse taxon name from the following line:\n $_\n\n");
	}
	
	my $char_line = substr ($_, length ($taxon));
	
	$char_line =~ s/\s//g;
	$char_line =~ s/\[[^\[\]]+\]//g; #remove any strings enclosed in square brackets (assumed to be comments or additional annotation)

	if ($char_count){
		unless ($char_count = length ($char_line)){
			die ("\nERROR: Found character matrix lines with different numbers of characters.\n");
		}
	}else{
		$char_count = length ($char_line);
	}
	
	if (exists ($char_hash{$taxon})){
		die ("\nERROR: duplicate taxa with name $taxon\.\n");
	}

	$char_hash{$taxon} = $char_line;

}

if ($nchar){
	$char_count == $nchar or print STDERR ("\nWARNING: Identified $char_count characters in matrix but stated nchar is $nchar\n\n");
}


extract_and_analyze_quartets(\%char_hash, $FORKS, $OUTPUT, $char_count, $PAUP, 1);

print "\n" . (localtime) . "\nCompleted generation of quartet trees. Starting conversion to MRP matrix and performing tree search.\n\n";

run_mrp_search ("$OUTPUT\.quartets.tre", \@taxa_list, $OUTPUT, $PAUP, $SEARCH);
my $mrp_tree = parse_paup_output_single("$OUTPUT\.MRP4.nex");
my $FH_MRPTREE = open_output ("$OUTPUT\.MRP.tre");
print $FH_MRPTREE $mrp_tree;
close $FH_MRPTREE;

unlink("$OUTPUT\.MRP1.nex");
unlink("$OUTPUT\.MRP2.nex");
unlink("$OUTPUT\.MRP3.nex");
unlink("$OUTPUT\.MRP4.nex");

if ($BS){
	
	print "\n" . (localtime) . "\nStarting bootstrap analysis.\n\n";
	
	my %char_HoH;
	
	for (my $i = 1; $i <= $BS; ++$i){
		for (my $j = 0; $j < $char_count; ++$j){
			my $rand_pos = int(rand($char_count));
			foreach my $species (keys %char_hash){
				$char_HoH{$i}->{$species} .= substr ($char_hash{$species}, $rand_pos, 1);
			}
		}		
	}
	
	foreach my $rep (sort {$a <=> $b} keys %char_HoH){
		extract_and_analyze_quartets($char_HoH{$rep}, $FORKS, "$OUTPUT\_BS$rep", $char_count, $PAUP, 0);		
	}	


	print "\n" . (localtime) . "\nCompleted generation of quartet trees for al $BS pseudoreplicates. Starting conversions to MRP matrices and performing tree searches on pseudoreplicate datasets.\n\n";

	my $pm_bs;
	unless ($FORKS == 1){
		$pm_bs = Parallel::ForkManager->new($FORKS);
	}

	foreach my $rep (sort {$a <=> $b} keys %char_HoH){
		my $pid;
		unless ($FORKS == 1){
			$pid = $pm_bs->start and next; #this is the forking process that will start multiple calls in parallel
		}

		run_mrp_search ("$OUTPUT\_BS$rep\.quartets.tre", \@taxa_list, "$OUTPUT\_BS$rep", $PAUP, $SEARCH);
		unlink("$OUTPUT\_BS$rep\.MRP1.nex");
		unlink("$OUTPUT\_BS$rep\.MRP2.nex");
		unlink("$OUTPUT\_BS$rep\.MRP3.nex");
		unless ($SAVE_REPS){
			unlink ("$OUTPUT\_BS$rep\.quartets.tre");
			unlink ("$OUTPUT\_BS$rep\.MRP_search.nex");			
		}
			
		unless ($FORKS == 1){
			$pm_bs->finish;
		}
	}	
	unless ($FORKS == 1){
		$pm_bs->wait_all_children;
	}
	
	my $FH_BSTREES = open_output ("$OUTPUT\.MRP_bs_pseudoreplicates.tre");
	foreach my $rep (sort {$a <=> $b} keys %char_HoH){
		print $FH_BSTREES parse_paup_output_single("$OUTPUT\_BS$rep\.MRP4.nex");
		unlink("$OUTPUT\_BS$rep\.MRP4.nex");
	}
	close $FH_BSTREES;

	my $FH_CONSENSUS_COMMAND = open_output ("$OUTPUT\.consensus_run.nex");
	
	print $FH_CONSENSUS_COMMAND "BEGIN PAUP;
set crit=parsimony;
set autoclose=yes warnreset=no notifybeep=no warntsave=no maxtrees=1000 increase=auto;
Gettrees file=$OUTPUT\.MRP_bs_pseudoreplicates.tre mode=3 unrooted=yes duptrees=keep warntree=no;
ConTree all/grpfreq=yes showtree=no strict=no majrule=yes le50=yes append=no treefile=$OUTPUT\.MRP_bs_consensus.temp.tre;
END;\nquit;\n";
	close $FH_CONSENSUS_COMMAND;

	system ("$PAUP $OUTPUT\.consensus_run.nex > /dev/null");
	my $contree = parse_paup_output_single("$OUTPUT\.MRP_bs_consensus.temp.tre");

	my $FH_CONTREE = open_output ("$OUTPUT\.MRP_bs_consensus.tre");
	print $FH_CONTREE $contree;
	close $FH_CONTREE;

	unlink ("$OUTPUT\.consensus_run.nex");
	unlink ("$OUTPUT\.MRP_bs_consensus.temp.tre");
}


print "\n" . (localtime) . "\nAnalysis Complete. Final tree written to $OUTPUT\.MRP.tre. Output quartet trees written to $OUTPUT\.quartets.tre.\n\n";

sub run_mrp_search{
	
	my ($quartets_file, $taxa_array_ref, $output_basename, $paup_exe, $search) = @_;
	
	my @all_taxa = @{$taxa_array_ref};
	my $taxa_count = scalar (@all_taxa);
	
	my $FH_MRP = open_output("$output_basename\.MRP1.nex");
	
	print $FH_MRP "\#NEXUS 
Begin data;
Dimensions ntax=$taxa_count nchar=1;
Matrix\n";
	
	foreach (@all_taxa){
		print $FH_MRP "$_ 0\n";
	}

	print $FH_MRP ";
end;

BEGIN PAUP;
set autoclose=yes warnreset=no notifybeep=no warntsave=no maxtrees=100 increase=auto;
Gettrees file=$quartets_file unrooted=yes duptrees=keep warntree=no;
matrixrep file=$output_basename\.MRP2.nex;
END;
QUIT;
";

	close $FH_MRP;
	system ("$paup_exe $output_basename\.MRP1.nex > /dev/null");
	
	my $FH_MRP2 = open_file ("$output_basename\.MRP2.nex");
	my $FH_MRP3 = open_output ("$output_basename\.MRP3.nex");
	
	while (<$FH_MRP2>){
		$_ =~ /^Begin\ assumptions;/ and last;
		print $FH_MRP3 $_;
	}
	
	print $FH_MRP3 "begin paup;
exclude uninf;
export file=$output_basename\.MRP_search.nex format=nexus;
end;
quit;
";

	close $FH_MRP2;
	close $FH_MRP3;
	
	system ("$paup_exe $output_basename\.MRP3.nex > /dev/null");

	my $FH_MRPSEARCH = open_output_append ("$output_basename\.MRP_search.nex");

	print $FH_MRPSEARCH "\n\nBEGIN PAUP;
set crit=parsimony;
set autoclose=yes warnreset=no notifybeep=no warntsave=no maxtrees=5000 increase=auto;
";
	
	if ($search eq "tbr"){
		print $FH_MRPSEARCH "hsearch addseq=random swap=tbr multrees=yes nreps=100 nchuck=50 chuckscore=1;\n";
	}elsif ($search eq "bandb"){
		print $FH_MRPSEARCH "bandb multrees=yes;\n";		
	}else{
		die ("\nError: unrecognized setting for --search ($search)\n\n");
	}

	print $FH_MRPSEARCH "ConTree all/grpfreq=no showtree=no strict=yes majrule=no append=yes treefile=$output_basename\.MRP4.nex;
END;
quit;\n";

	close $FH_MRPSEARCH;

	system ("$paup_exe $output_basename\.MRP_search.nex > /dev/null");


}


sub extract_and_analyze_quartets{

	my ($matrix_hash_ref, $forks, $output_base, $character_count, $paup, $keep_last_log) = @_;

	my %matrix_hash = %{$matrix_hash_ref};

	my $FHO = open_output ("$output_base\.quartets.tre");
	my @taxa = sort { "\L$a" cmp "\L$b" } keys %matrix_hash;

	my $pm;
	unless ($FORKS == 1){
		$pm = Parallel::ForkManager->new($FORKS);
	}


	for (my $i = 0; $i < scalar(@taxa) - 3; ++$i){
		for (my $j = $i + 1; $j < scalar(@taxa) - 2; ++$j){
			for (my $k = $j + 1; $k < scalar(@taxa) - 1; ++$k){
				for (my $l = $k + 1; $l < scalar(@taxa); ++$l){
				
					my $pid;
					unless ($FORKS == 1){
						$pid = $pm->start and next; #this is the forking process that will start multiple calls in parallel
					}
					#print nexus file for given quartet
					print_nexus($taxa[$i], $taxa[$j], $taxa[$k], $taxa[$l], $output_base, "$i\_$j\_$k\_$l", $character_count, $matrix_hash_ref);
				
					#run paup
					if ($keep_last_log){
						system ("$paup $output_base\.$i\_$j\_$k\_$l\.temp.nex > $output_base\.last_quartet_log.txt");
					}else{
						system ("$paup $output_base\.$i\_$j\_$k\_$l\.temp.nex > /dev/null");
					}
					
					#parse paup output and print trees
					my $tree_string = parse_paup_output ("$output_base\.$i\_$j\_$k\_$l\.temp.tre");
					print $FHO $tree_string;
				
					#delete nexus quartet file and paup output
					unlink ("$output_base\.$i\_$j\_$k\_$l\.temp.nex");
					unlink ("$output_base\.$i\_$j\_$k\_$l\.temp.tre");

					unless ($FORKS == 1){
						$pm->finish;
					}
				}
			}	
		}
	}
	unless ($FORKS == 1){
		$pm->wait_all_children;
	}
}



sub print_nexus{
	
    my ($taxon1, $taxon2, $taxon3, $taxon4, $output, $rep_string, $num_chars, $hashRef) = @_;
	my $FH = open_output ("$output\.$rep_string\.temp.nex");
    my %hash_for_nexus = %{$hashRef};
    
    print $FH "BEGIN PAUP;
set autoclose=yes warnreset=no notifybeep=no warntsave=no maxtrees=100 increase=auto;
set crit=parsimony;
END;

#NEXUS
BEGIN DATA;
	DIMENSIONS  NTAX=4 NCHAR=$num_chars;
	FORMAT DATATYPE=STANDARD  SYMBOLS=\"0 1 2 3 4 5 6 7 8 9\" MISSING=? GAP=-  INTERLEAVE ;
MATRIX
";

print $FH "$taxon1 $hash_for_nexus{$taxon1}\n$taxon2 $hash_for_nexus{$taxon2}\n$taxon3 $hash_for_nexus{$taxon3}\n$taxon4 $hash_for_nexus{$taxon4}\n";

print $FH ";
end;

begin paup;
pset collapse=no;
alltrees;
Savetrees file=$output\.$rep_string\.temp.tre format=altnexus root=no;
end;
quit;";

	close $FH;

}


sub parse_paup_output{
	my $paup_tree_file = shift (@_);
	my @tree_lines = file_to_array($paup_tree_file);
	my @trees;
	foreach (@tree_lines){
		if ($_ =~ /tree.*\[\&U\]\ (.*\n)/){
			push (@trees, $1);
		}
	}
	
	my $return_string;
	if (scalar (@trees) == 1){
		$return_string = $trees[0] . $trees[0] . $trees[0] . $trees[0] . $trees[0] . $trees[0];
	}elsif (scalar (@trees) == 2){
		$return_string = $trees[0] . $trees[0] . $trees[0] . $trees[1] . $trees[1] . $trees[1];	
		#print STDERR "2 equally parsimonious trees:", substr ($trees[0], 0, -1), "\t",  substr ($trees[1], 0, -1), "\n";
	}
	elsif (scalar (@trees) == 3){
		$return_string = $trees[0] . $trees[0] . $trees[1] . $trees[1] . $trees[2] . $trees[2];	
		#print STDERR "3 equally parsimonious trees:", substr ($trees[0], 0, -1), "\t",  substr ($trees[1], 0, -1), "\t",  substr ($trees[2], 0, -1), "\n";
	}else{
		die ("\nERROR: PAUP search returned unexpected number of equally parsimonius trees (0 or >3)\n\n");
	}
	
	return $return_string;
}

sub parse_paup_output_single{
	my $paup_tree_file = shift (@_);
	my @tree_lines = file_to_array($paup_tree_file);
	my @trees;
	foreach (@tree_lines){
		if ($_ =~ /tree.*\[\&U\]\ (.*\n)/){
			push (@trees, $1);
		}
	}

	my $return_string;
	if (scalar (@trees) == 1){
		$return_string = $trees[0];
	}else{
		die ("\nERROR: PAUP search returned unexpected number of trees (0 or >1)\n\n");
	}

	return $return_string;

}

sub open_output {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output;

    unless(open($fh_output, ">$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output;
}

sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}

sub open_file {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh;

    unless(open($fh, $filename)) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh;
}

sub open_output_append {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output_append;

    unless(open($fh_output_append, ">>$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output_append;
}