# SDPquartets
A Perl script for implementing Split Decomposition with Parsimony and quartet assembly (SDPquartets) 

`SDPquartets.pl`

## Overview: 
This script automates the multi-step process of taking a binary character matrix (such as presence/absence of retroelement insertions), building parsimony trees for every possible quartet of species within the taxon sample, and using [matrix representation with parsimony (MRP)](https://www.sciencedirect.com/science/article/pii/105579039290035F) to infer an overall species tree based on the set of quartet trees. It also will perform bootstrap resampling of the original character matrix to estimate support for the inferred species relationships. All parsimony searches are implemented in [PAUP\*](https://paup.phylosolutions.com/).

## Requirements: 

This automation is implemented with a Perl script that has been designed for a Unix environment (Mac OSX or Linux, not Windows). It has been tested in Mac OSX 10.14.6 and Linux CentOS 6, but it should work in most Unix environments.

- **Perl modules:** The provided Perl script should be called by users (`SDPquartets.pl`). Perl is pre-installed in most Mac OSX and Linux distributions. If the `--forks` option is to be used for parallelization, the [Parallel::ForkManager Perl module](https://metacpan.org/pod/Parallel::ForkManager) must be installed.

- **PAUP\*:** The script calls the the [command-line version of PAUP\*](http://phylosolutions.com/paup-test/). It has been tested on paup4a166 but is expected to work on the latest version. It does NOT work with the older 4.0b10 version. If the PAUP\* binaries are installed in a directory in your PATH, you can simply refer to the executable name with the `--paup` argument. If not, provide the full path and executable name with that argument.

## Running SDPquartets.pl:
The script can be called from the command line to analyze an input character matrix (in Nexus format; see sample\_data subdirectory). The user must also specify an output name, the name/location of the Paup\* command-line executable, and desired options for bootstrapping, parallelization and search methods.

Usage: `perl SDPquartets.pl [arguments]`
   
   REQUIRED ARGUMENTS
   
   Character matrix input:
   
         --matrix     - file containing character matrix in Nexus format

   Output:
   
         --output     - base name for output files (extensions will be appended)

   PAUP\* executable:
   
         --paup       - name (and path if necessary) of PAUP* executable


   OPTIONAL ARGUMENTS

   Bootstrap:
   
         --bs         - number of bootstrap pseudoreplicates to generate
 
   Retain bootstrap replicates:
   
         --save_reps  - Add this flag to retain resampled bootstrap matrices
         				and quartet trees.
         
   Parallelization:
   
         --forks      - number of forks to run in parallel (default = 1).
                        This option requires the Parallel::ForkManager
                        Perl module. The optimal number of forks will depend
                        on the dataset and the hardware available. We have
                        found that forks can be set well in excess of the 
                        available number of cores and still produce further
                        speed-ups

   Tree search method:
   
         --search     - Specify "bandb" to change the final tree search method
                        from the default of a heuristic search with TBR to branch
                        and bound. If bandb is not turned on, a TBR search will
                        be run with 100 random addition sequences and up to 50
                        equally parsimonious trees held.  [default: tbr] 
   
   EXAMPLE
   
   `perl SDPquartets.pl --matrix=sample_data/input.nex --output=my_sample_output --paup=/usr/local/paup/paup4a166_centos64 --bs=1000 --forks=40 --search=bandb`
              
## Output:
The script will produce multiple output files with different extensions appended to the base name (e.g., "OUTPUT") provided with `--output`.

- **OUTPUT.MRP.tre:** The final parsimony species tree (or strict consensus of multiple equally parsimonious trees) in newick format.
- **OUTPUT.quartets.tre:** A file that contains every quartet tree (in newick format). To address the possibility of equally parsimonious ties for a given quartet of species, the script applies a weighting by writing 6 copies of any quartet tree that is the single most parsimonious resolution, 3 copies each for the best trees when there are two equally most-parsimonious trees, and 2 copies each when all three trees for a quartet of species are equally parsimonious.
- **OUTPUT.MRP_search.nex:** The Nexus file containing the data matrix created from quartet trees to perform the MRP search for the final species tree.
- **OUTPUT.last_quartet_log.txt:** This files contains the stdout log produced by PAUP\* during the last quartet-tree parsimony search. It can be checked for warning or errors if runs do not appear to be proceeding properly. Note that PAUP\* logs are not saved for other parts of the analysis such as the MRP search.

Additional output files are generated if the `--bs` option is used to run bootstrapping.

- **OUTPUT.MRP_bs_pseudoreplicates.tre:** A file with the final species tree (in newick format) for each bootstrap replicate
- **OUTPUT.MRP_bs_consensus.tre:** A file summarizing bootstrap support values (extended majority rule consensus of all bootstrap trees in newick format).

If bootstrapping is run and the `--save_reps` options is specified, two additional files are retained after the end of the run for every bootstrap replicate (where replicate number is indicated by XXX).

- **OUTPUT_BSXXX.quartets.tre:** Same as OUTPUT.quartets.tre (see above) except for the bootstrap resampling.
- **OUTPUT_BSXXX.MRP_search.nex:** Same as OUTPUT.MRP_search.nex (see above) except for the bootstrap resampling. Note that this is the MRP matrix summarizing quartet trees. The resampled matrices of the original binary characters (e.g., retroelement data) are only stored in memory during the run and not saved to file.
