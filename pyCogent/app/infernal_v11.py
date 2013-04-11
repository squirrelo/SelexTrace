#!/usr/bin/env python
"""
Provides an application controller for the commandline version of:
Infernal 1.1r2 only.
"""
from cogent.app.parameters import FlagParameter, ValuedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.rfam import MinimalRfamParser, ChangedSequence, \
    ChangedRnaSequence, ChangedDnaSequence
from cogent.parse.infernal_v11 import CmsearchParser
from cogent.core.moltype import DNA, RNA
from cogent.core.alignment import SequenceCollection, Alignment, DataError
from cogent.format.stockholm import stockholm_from_alignment
from cogent.struct.rna2d import ViennaStructure, wuss_to_vienna
from os import remove

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann, Joshua Shorenstein"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Joshua Shorenstein"
__email__ = "joshua.shorenstein@colorado.edu"
__status__ = "Development"

MOLTYPE_MAP = {'DNA':'--dna',\
                DNA:'--dna',\
               'RNA':'--rna',\
                RNA:'--rna',\
               }
               
SEQ_CONSTRUCTOR_MAP = {'DNA':ChangedDnaSequence,\
                        DNA:ChangedDnaSequence,\
                       'RNA':ChangedRnaSequence,\
                        RNA:ChangedRnaSequence,\
                       }

class Cmalign(CommandLineApplication):
    """cmalign application controller."""
    _options = {
 
    # -o <f> Save the alignment in Stockholm format to a file <f>. The default
    #   is to write it to standard output.
    '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' '),\
    
    # -l Turn on the global alignment algorithm. Default is local.
    '-g':FlagParameter(Prefix='-',Name='g'),\
    
    # Expert Options 
    
    # --optacc Align sequences using the Durbin/Holmes optimal accuracy 
    #   algorithm. This is default behavior, so this option is probably useless. 
    '--optacc':FlagParameter(Prefix='--',Name='optacc'),\
    
    # --cyk Do not use the Durbin/Holmes optimal accuracy alignment to align the 
    #   sequences, instead use the CYK algorithm which determines the optimally
    #   scoring alignment of the sequence to the model. 
    '--cyk':FlagParameter(Prefix='--',Name='cyk'),\
    
    # --sample Sample an alignment from the posterior distribution of
    #   alignments.
    '--sample':FlagParameter(Prefix='--',Name='sample'),\

    # --seed <n> Set the random number generator seed to <n>, where <n> is a 
    #   positive integer. This option can only be used in combination with 
    #   --sample. The default is to use time() to generate a different seed for
    #   each run, which means that two different runs of cmalign --sample on the
    #   same alignment will give slightly different results. You can use this
    #   option to generate reproducible results.
    '--seed':ValuedParameter(Prefix='-',Name='seed',Delimiter=' '),\
    
    # --notrunc Do not use the truncated alignment algorithm
    '--notrunc':FlagParameter(Prefix='--',Name='notrunc'),\
    
    # --sub Turn on the sub model construction and alignment procedure.
    '--sub':FlagParameter(Prefix='--',Name='sub'),\
    
    # --small Use the divide and conquer CYK alignment algorithm described in 
    #   SR Eddy, BMC Bioinformatics 3:18, 2002.
    '--small':FlagParameter(Prefix='--',Name='small'),\
    
    # --hbanded This option is turned on by default. Accelerate alignment by
    #   pruning away regions of the CM DP matrix that are deemed negligible by
    #   an HMM.
    '--hbanded':FlagParameter(Prefix='--',Name='hbanded'),\
    
    # --nonbanded Turns off HMM banding.
    '--nonbanded':FlagParameter(Prefix='--',Name='nonbanded'),\
    
    # --tau <x> Set the tail loss probability used during HMM band calculation
    #   to <x>.
    '--tau':ValuedParameter(Prefix='--',Name='tau',Delimiter=' '),\
    
    # --mxsize <x> Set the maximum allowable DP matrix size to <x> megabytes.
    '--mxsize':ValuedParameter(Prefix='--',Name='mxsize',Delimiter=' '),\
    
    # --dnaout Output the alignments as DNA sequence alignments.
    '--dnaout':FlagParameter(Prefix='--',Name='dnaout'),\
    
    # --sfile <f> output alignment score information to <f>
    '--sfile':ValuedParameter(Prefix='--',Name='sfile',Delimiter=' '),\

    # --tfile <f> output individual sequence parsetrees to <f>
    '--tfile':ValuedParameter(Prefix='--',Name='tfile',Delimiter=' '),\

    # --ifile <f> output information on per-sequence inserts to <f>
    '--ifile':ValuedParameter(Prefix='--',Name='ifile',Delimiter=' '),\

    # --efile <f> output information on per-sequence  EL inserts to <f>
    '--efile':ValuedParameter(Prefix='--',Name='efile',Delimiter=' '),\

    # --mapali <f> include alignment that CM came from in file <f>
    '--mapali':FlagParameter(Prefix='--',Name='mapali'),\

    # --regress <f> save regression test data to file <f>
    '--mapali':ValuedParameter(Prefix='--',Name='mapali',Delimiter=' '),\

    # --mapstr include structure in alignment file output by --mapali
    '--mapstr':FlagParameter(Prefix='--',Name='mapstr'),\

    # --informat <s> Make sure input sequence file matches input format <s>
    #Available formats: FASTA, Genbank
    '--informat':ValuedParameter(Prefix='--',Name='informat',Delimiter=' '),\

    # --outformat <f> Make sure input sequence file matches input format <s>
    #Available formats: Stockholm, Pfam, AFA (aligned FASTA), A2M, Clustal, PHYLIP
    '--outformat':ValuedParameter(Prefix='--',Name='outformat',Delimiter=' '),\

    # --noprob do not include posterior probabilities in alignment
    '--noprob':FlagParameter(Prefix='--',Name='noprob'),\

    # --cpu <n> Number of parallel threads to use while computing n >= 1 (default all cores)
    '--outformat':ValuedParameter(Prefix='--',Name='outformat',Delimiter=' '),\
    
    
    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmalign"
    _suppress_stderr=True

    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str
    
    def _tempfile_as_multiline_string(self, data):
        """Write a multiline string to a temp file and return the filename.

            data: a multiline string to be written to a file.

           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

        """
        filename = FilePath(self.getTmpFilename(self.TmpDir))
        data_file = open(filename,'w')
        data_file.write(data)
        data_file.close()
        return filename

    def _alignment_out_filename(self):
        
        if self.Parameters['-o'].isOn():
            refined_filename = self._absolute(str(\
                self.Parameters['-o'].Value))
        else:
            raise ValueError, 'No alignment output file specified.'
        return refined_filename

    def _get_result_paths(self,data):
        result = {}
        if self.Parameters['-o'].isOn():
            out_name = self._alignment_out_filename()
            result['Alignment'] = ResultPath(Path=out_name,IsWritten=True)
        
        return result

class Cmbuild(CommandLineApplication):
    """cmbuild application controller."""
    _options = {
    
    # -n <s> Name the covariance model <s>. (Does not work if alifile contains
    #   more than one alignment).
    '-n':ValuedParameter(Prefix='-',Name='n',Delimiter=' '),\
    
    # -F Allow cmfile to be overwritten. Normally, if cmfile already exists,
    #   cmbuild exits with an error unless the -A or -F option is set.
    '-F':FlagParameter(Prefix='-',Name='F'),\

    # -o <f> Direct summary data to file <f> instead of stdout
    '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' '),\
    
    # -O <s> Resave consensus column annotated MSA file to <f>
    '-O':ValuedParameter(Prefix='-',Name='O',Delimiter=' '),\
    
    # Expert Options
    # --fast assign columns with >= symfrac residues as consensus
    '--fast':FlagParameter(Prefix='--',Name='fast'),\

    # --hand Use reference coordinate annotation to specify consensus.
    '--hand':FlagParameter(Prefix='--',Name='hand'),\

    # --symfrac <x>  Fraction of non-gaps to require in a consensus column.
    # Default is 0.1
    '--symfrac':ValuedParameter(Prefix='--',Name='symfrac',Delimiter=' '),\

    # --noss Ignore secondary structure in input alignment
    '--noss':FlagParameter(Prefix='--',Name='noss'),\

    # --rsearch <f> Parameterize emission scores a la RSEARCH, using the
    #   RIBOSUM matrix in file <f>. For more information see the RSEARCH 
    #   publication (Klein and Eddy, BMC Bioinformatics 4:44, 2003). Actually,
    #   the emission scores will not exactly With --rsearch enabled, all
    #   alignments in alifile must contain exactly one sequence or the --call
    #   option must also be enabled.
    '--rsearch':ValuedParameter(Prefix='--',Name='rsearch',Delimiter=' '),\
    
    # --wgsc Use the Gerstein/Sonnhammer/Chothia (GSC) weighting algorithm.
    #   This is the default unless the number of sequences in the alignment
    #   exceeds a cutoff (see --pbswitch), in which case the default becomes
    #   the faster Henikoff position-based weighting scheme.
    '--wgsc':FlagParameter(Prefix='--',Name='wgsc'),\
    
    # --wblosum Use the BLOSUM filtering algorithm to weight the sequences,
    #   instead of the default GSC weighting.
    '--wblosum':FlagParameter(Prefix='--',Name='wblosum'),\
    
    # --wpb Use the Henikoff position-based weighting scheme. This weighting
    #   scheme is automatically used (overriding --wgsc and --wblosum) if the
    #   number of sequences in the alignment exceeds a cutoff (see --pbswitch).
    '--wpb':FlagParameter(Prefix='--',Name='wpb'),\
    
    # --wnone Turn sequence weighting off; e.g. explicitly set all sequence
    #   weights to 1.0.
    '--wnone':FlagParameter(Prefix='--',Name='wnone'),\
    
    # --wgiven Use sequence weights as given in annotation in the input
    #   alignment file. If no weights were given, assume they are all 1.0.
    #   The default is to determine new sequence weights by the Gerstein/
    #   Sonnhammer/Chothia algorithm, ignoring any annotated weights.
    '--wgiven':FlagParameter(Prefix='--',Name='wgiven'),\
    
    # --wid <x> Controls the behavior of the --wblosum weighting option by
    #   setting the percent identity for clustering the alignment to <x>.
    '--wid':ValuedParameter(Prefix='--',Name='wid',Delimiter=' '),\
    
    # --eent Use the entropy weighting strategy to determine the effective
    #   sequence number that gives a target mean match state relative entropy.
    '--eent':FlagParameter(Prefix='--',Name='eent'),\
    
    # --enone Turn off the entropy weighting strategy. The effective sequence
    #   number is just the number of sequences in the alignment.
    '--enone':FlagParameter(Prefix='--',Name='enone'),\
    
    # --ere <x> Set the target mean match state entropy as <x>. By default the
    #   target entropy 1.46 bits.
    '--ere':ValuedParameter(Prefix='--',Name='ere',Delimiter=' '),\
    
    # --eset <x> Set effective number of sequences for all models to <x>
    '--eset':ValuedParameter(Prefix='--',Name='eset',Delimiter=' '),\

    # --eminseq <x> Set min effective number of sequences for all models to <x>
    '--eminseq':ValuedParameter(Prefix='--',Name='eminseq',Delimiter=' '),\

    # --ehmmre <x> Set minimum relative entropy to <x>, requires --eent
    '--ehmmre':ValuedParameter(Prefix='--',Name='ehmmre',Delimiter=' '),\

    # --esigma <x> Set sigma parameter to <x>, requires --eent
    '--esigma':ValuedParameter(Prefix='--',Name='esigma',Delimiter=' '),\

    # --null <f> Read a null model from <f>. The null model defines the
    #   probability of each RNA nucleotide in background sequence, the default
    #   is to use 0.25 for each nucleotide.
    '--null':ValuedParameter(Prefix='--',Name='null',Delimiter=' '),\
    
    # --prior <f> Read a Dirichlet prior from <f>, replacing the default mixture
    #   Dirichlet.
    '--prior':ValuedParameter(Prefix='--',Name='prior',Delimiter=' '),\
    
    # --p7ere <x> for the filter p7 HM, set minimum relative entropy/pos to <x>
    '--p7ere':ValuedParameter(Prefix='--',Name='p7ere',Delimiter=' '),\

    # --p7ml define the filter p7 HMM as the ML p7 HMM
    '--p7ml':FlagParameter(Prefix='--',Name='p7ml'),\

    # --EmN <x> Number of sampled sequences to use for p7 local MSV calibration
    '--EmN':ValuedParameter(Prefix='--',Name='EmN',Delimiter=' '),\

    # --EvN <x> Number of sampled sequences to use for p7 local Viturbi
    # calibration
    '--EvN':ValuedParameter(Prefix='--',Name='EvN',Delimiter=' '),\

    # --ElfN <x> Number of sampled sequences to use for p7 local Forward
    # calibration
    '--ElfN':ValuedParameter(Prefix='--',Name='ElfN',Delimiter=' '),\

    # --EgfN <x> Number of sampled sequences to use for p7 glocal Forawrd
    # calibration
    '--EgfN':ValuedParameter(Prefix='--',Name='EgfN',Delimiter=' '),\

    # --refine <f> Attempt to refine the alignment before building the CM using
    #   expectation-maximization (EM). The final alignment (the alignment used
    #   to build the CM that gets written to cmfile) is written to <f>.
    '--refine':ValuedParameter(Prefix='--',Name='refine',Delimiter=' '),\

    # -l With --refine, turn on the local alignment algorithm, which allows the
    #   alignment to span two or more subsequences if necessary (e.g. if the
    #   structures of the query model and target sequence are only partially
    #   shared), allowing certain large insertions and deletions in the
    #   structure to be penalized differently than normal indels. The default is
    #   to globally align the query model to the target sequences.
    '-l':ValuedParameter(Prefix='-',Name='l',Delimiter=' '),\
    
    # --gibbs Modifies the behavior of --refine so Gibbs sampling is used
    #   instead of EM.
    '--gibbs':FlagParameter(Prefix='--',Name='gibbs'),\
    
    # --seed <n> Set the random seed to <n>, where <n> is a positive integer.
    #   This option can only be used in combination with --gibbs. The default is
    #   to use time() to generate a different seed for each run, which means
    #   that two different runs of cmbuild --refine <f> --gibbs on the same
    #   alignment will give slightly different results. You can use this option
    #   to generate reproducible results.
    '--seed':ValuedParameter(Prefix='--',Name='seed',Delimiter=' '),\
    
    
    # --cyk With --refine, align with the CYK algorithm.
    '--cyk':FlagParameter(Prefix='--',Name='cyk'),\
    
    # --notrunc With --refine, do not use truncated alignment algorithm
    '--notrunc':FlagParameter(Prefix='--',Name='notrunc'),\
    
    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmbuild"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str
    
    def _refine_out_filename(self):
        
        if self.Parameters['--refine'].isOn():
            refined_filename = self._absolute(str(\
                self.Parameters['--refine'].Value))
        else:
            raise ValueError, 'No refine output file specified.'
        return refined_filename
    
    def _cm_out_filename(self):
        
        if self.Parameters['-n'].isOn():
            refined_filename = self._absolute(str(\
                self.Parameters['-n'].Value))
        else:
            raise ValueError, 'No cm output file specified.'
        return refined_filename
    
    def _tempfile_as_multiline_string(self, data):
        """Write a multiline string to a temp file and return the filename.

            data: a multiline string to be written to a file.

           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

        """
        filename = FilePath(self.getTmpFilename(self.TmpDir))
        data_file = open(filename,'w')
        data_file.write(data)
        data_file.close()
        return filename
    
    def _get_result_paths(self,data):
        result = {}
        if self.Parameters['--refine'].isOn():
            out_name = self._refine_out_filename()
            result['Refined'] = ResultPath(Path=out_name,IsWritten=True)
        if self.Parameters['-n'].isOn():
            cm_name = self._cm_out_filename()
            result['CmFile'] = ResultPath(Path=cm_name,IsWritten=True)
        
        return result

    
class Cmcalibrate(CommandLineApplication):
    """cmcalibrate application controller."""
    _options = {
    # -L <x> Set the length of random sequence to search for the CM
    #  to <x> megabases (Mb) where 0.01<=x<=160. Default is 1.6
    '-L':ValuedParameter(Prefix='-',Name='L',\
        Delimiter=' '),\
    
    # -seed <n> Set the random number generator seed to <n>, where <n> is a
    #   positive integer. The default is to use time() to generate a different
    #   seed for each run, which means that two different runs of cmcalibrate on 
    #   the same CM will give slightly different E-value and HMM filter
    #   threshold parameters. You can use this option to generate reproducible
    #   results.
    '--seed':ValuedParameter(Prefix='--',Name='seed',Delimiter=' '),\
    
    # --nforecast <n> Predict the running time of the calibration for cmfile and
    #   <n> processors and exit, DO NOT perform the calibration.
    '--nforecast':ValuedParameter(Prefix='--',Name='nforecast',Delimiter=' '),\

    # --forecast Predict the running time of the calibration for cmfile and
    # exit, DO NOT perform the calibration.
    '--forecast':FlagParameter(Prefix='--',Name='forecast'),\

    # --memreq Predict the memory required for the calibration for cmfile and
    # exit, DO NOT perform the calibration.
    '--memreq':FlagParameter(Prefix='--',Name='memreq'),\

    # Expert Options
    
    # --gtailn <n> fit the top <n> hits/Mb in histogram for glocal modes.
    # Default is 250
    '--gtailn':ValuedParameter(Prefix='--',Name='gtailn',Delimiter=' '),\

    # --ltailn <n> fit the top <n> hits/Mb in histogram for local modes.
    # Default is 750
    '--ltailn':ValuedParameter(Prefix='--',Name='ltailn',Delimiter=' '),\

    # --tailp <n> set the fraction of histogram tail to fit to exponental
    # tail to <n>
    '--tailp':ValuedParameter(Prefix='--',Name='tailp',Delimiter=' '),\

    # --hfile <f> save fitted score histogram(s) to file <f>
    '--hfile':ValuedParameter(Prefix='--',Name='hfile',Delimiter=' '),\

    # --sfile <f> save survival plot to file <f>
    '--sfile':ValuedParameter(Prefix='--',Name='sfile',Delimiter=' '),\

    # --qqfile <f> save save Q-Q plit for score histograms to file <f>
    '--qqfile':ValuedParameter(Prefix='--',Name='qqfile',Delimiter=' '),\

    # --ffile <f> save lambdas for different tail fit probabilities
    # to file <f>
    '--ffile':ValuedParameter(Prefix='--',Name='ffile',Delimiter=' '),\

    # --xfile <f> save scores in fit tail to file <f>
    '--xfile':ValuedParameter(Prefix='--',Name='xfile',Delimiter=' '),\

    # --beta <x> Set tail loss probability for query dependent banding
    # (QDB) to <x>
    '--beta':ValuedParameter(Prefix='--',Name='beta',Delimiter=' '),\

    # --nonbanded Do not use QDB
    '--nonbanded':FlagParameter(Prefix='--',Name='nonbanded'),\

    # --nonull3 turn OFF the NULL3 post hoc additional null model
    '--nonull3':FlagParameter(Prefix='--',Name='nonull3'),\

    # --random use GC content of random null background model of CM
    '--random':FlagParameter(Prefix='--',Name='random'),\

    # --gc <f> Use GC content distribution from file <f>
    '--gc':ValuedParameter(Prefix='--',Name='gc',Delimiter=' '),\

    # --cpu <n> Use <n> parallel CPU workers for computation
    '--cpu':ValuedParameter(Prefix='--',Name='cpu',Delimiter=' '),\

    }

    _parameters = {}
    _parameters.update(_options)
    _command = "cmcalibrate"
    #_suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str

class Cmemit(CommandLineApplication):
    """cmemit application controller."""
    _options = {

    # -o <f> Send output to file <f>, not stdout
    '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' '),\

    # -N <n> Generate <n> sequences. Default is 10.
    '-N':ValuedParameter(Prefix='-',Name='N',Delimiter=' '),\

    # --rna Output as RNA sequences (Default)
    '--rna':FlagParameter(Prefix='--',Name='rna'),\

    # --dna Output as DNA sequences
    '--dna':FlagParameter(Prefix='--',Name='dna'),\
    
    # -u Write the generated sequences in unaligned format (FASTA). This is the
    # default, so this option is probably useless.
    '-u':FlagParameter(Prefix='-',Name='u'),\
    
    # -a Write the generated sequences in an aligned format (STOCKHOLM) with
    #   consensus structure annotation rather than FASTA.
    '-a':FlagParameter(Prefix='-',Name='a'),\

    # --outformat <s> used with -a to specify what output format to use
    #Default Stockholm format
    #Available formats: Stockholm,Pfam,AFA (aligned FASTA),A2M,Clustal,PHYLIP
    '--o':ValuedParameter(Prefix='--',Name='o',Delimiter=' '),\
    
    # -c Predict a single majority-rule consensus sequence instead of sampling
    #   sequences from the CM's probability distribution.
    '-c':FlagParameter(Prefix='-',Name='c'),\

    # -e <n> Embed emitted sequences in larger random sequences of length <n>
    '-e':ValuedParameter(Prefix='-',Name='e',Delimiter=' '),\

    # -iid Used with -e, generate sequencs as 25% ACGU
    '-e':FlagParameter(Prefix='-',Name='e'),\
    
    # -l Configure the CMs into local mode before emitting sequences. See the
    #   User's Guide for more information on locally configured CMs.
    '-l':FlagParameter(Prefix='-',Name='l'),\
    
    # --seed <n> Set the random seed to <n>, where <n> is a positive integer. The
    #   default is to use time() to generate a different seed for each run,
    #   which means that two different runs of cmemit on the same CM will give
    #   different results. You can use this option to generate reproducible
    #   results.
    '--seed':ValuedParameter(Prefix='--',Name='seed',Delimiter=' '),\
    
    # --u5p truncate unaligned sequences 5' end at random start position
    '--u5p':FlagParameter(Prefix='--',Name='u5p'),\

    # --u3p truncate unaligned sequences 3' end at random start position
    '--u3p':FlagParameter(Prefix='--',Name='u3p'),\

    # --a5p <n> Truncate 5' of aligned sequences at column <n>
    # set n to 0 for random
    '--a5p':ValuedParameter(Prefix='--',Name='a5p',Delimiter=' '),\

    # --a3p <n> Truncate 3' of aligned sequences at column <n>
    # set n to 0 for random
    '--a3p':ValuedParameter(Prefix='--',Name='a3p',Delimiter=' '),\

    # --idx <n> Start sequence numbering at <n>
    '--idx':ValuedParameter(Prefix='--',Name='idx',Delimiter=' '),\

    # --tfile <f> dump parsetrees to file <f>
    '--tfile':ValuedParameter(Prefix='--',Name='tfile',Delimiter=' '),\

    # --exp <x> exponentiate CM probabilities by <x>, then emit
    '--tfile':ValuedParameter(Prefix='--',Name='tfile',Delimiter=' '),\

    # --hmmonly Emit using HMM only, no CM
    '--hmmonly':FlagParameter(Prefix='--',Name='hmmonly'),\

    # --nohmmonly Emit using CM only, no HMM
    '--nohmmonly':FlagParameter(Prefix='--',Name='nohmmonly'),\

 
    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmemit"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str
    
class Cmscore(CommandLineApplication):
    """cmscore application controller."""
    _options = {
    
    # -n <n> Set the number of sequences to generate and align to <n>. This
    #   option is incompatible with the --infile option.
    '-n':ValuedParameter(Prefix='-',Name='n',Delimiter=' '),\
    
    # -l Turn on the local alignment algorithm, which allows the alignment to
    #   span two or more subsequences if necessary (e.g. if the structures of
    #   the query model and target sequence are only partially shared), allowing
    #   certain large insertions and deletions in the structure to be penalized
    #   differently than normal indels. The default is to globally align the
    #   query model to the target sequences.
    '-l':FlagParameter(Prefix='-',Name='l'),\
    
    # -s <n> Set the random seed to <n>, where <n> is a positive integer. The
    #   default is to use time() to generate a different seed for each run,
    #   which means that two different runs of cmscore on the same CM will give
    #   different results. You can use this option to generate reproducible
    #   results. The random number generator is used to generate sequences to
    #   score, so -s is incompatible with the --infile option which supplies
    #   the sequences to score in an input file.
    '-s':ValuedParameter(Prefix='-',Name='s',Delimiter=' '),\
    
    # -a Print individual timings and score comparisons for each sequence in
    #   seqfile. By default only summary statistics are printed.
    '-a':FlagParameter(Prefix='-',Name='a'),\
    
    # --sub Turn on the sub model construction and alignment procedure.
    '--sub':FlagParameter(Prefix='--',Name='sub'),\
    
    # --mxsize <x> Set the maximum allowable DP matrix size to <x> megabytes.
    '--mxsize':ValuedParameter(Prefix='--',Name='mxsize',Delimiter=' '),\
    
    # Expert Options
    
    # --emit Generate sequences to score by sampling from the CM.
    '--emit':FlagParameter(Prefix='--',Name='emit'),\
    
    # --random Generate sequences to score by sampling from the CMs null
    #   distribution. This option turns the --emit option off.
    '--random':FlagParameter(Prefix='--',Name='random'),\
    
    # --infile <f> Sequences to score are read from the file <f>. All the
    #   sequences from <f> are read and scored, the -n and -s options are
    #   incompatible with --infile.
    '--infile':ValuedParameter(Prefix='--',Name='infile',Delimiter=' '),\
    
    # --outfile <f> Save generated sequences that are scored to the file <f> in
    #   FASTA format. This option is incompatible with the --infile option.
    '--outfile':ValuedParameter(Prefix='--',Name='outfile',Delimiter=' '),\
    
    # --Lmin <n1> Must be used in combination with --random and --Lmax <n2>.
    #sets min length of ramdom sequences tp <n1>
    '--Lmin':ValuedParameter(Prefix='--',Name='Lmin',Delimiter=' '),\

    # --Lmax <n2> Must be used in combination with --random and --Lmin <n1>.
    #sets max length of ramdom sequences tp <n1>
    '--Lmax':ValuedParameter(Prefix='--',Name='Lmin',Delimiter=' '),\
    
    # --pad Must be used in combination with --emit and --search. Add <n> cm->W
    #   (max hit length) minus L (sequence <x> length) residues to the 5' and 3'
    #   end of each emitted sequence <x>.
    '--pad':FlagParameter(Prefix='--',Name='pad'),\
    
    # --hbanded Specify that the second stage alignment algorithm be HMM banded
    #   CYK. This option is on by default.
    '--hbanded':FlagParameter(Prefix='--',Name='hbanded'),\
    
    # --tau <x> For stage 2 alignment, set the tail loss probability used during
    #   HMM band calculation to <x>.
    '--tau':ValuedParameter(Prefix='--',Name='tau',Delimiter=' '),\
    
    # --aln2bands With --search, when calculating HMM bands, use an HMM
    #   alignment algorithm instead of an HMM search algorithm.
    '--aln2bands':FlagParameter(Prefix='--',Name='aln2bands'),\
    
    # --hsafe For stage 2 HMM banded alignment, realign any sequences with a
    #   negative alignment score using non-banded CYK to guarantee finding the
    #   optimal alignment.
    '--hsafe':FlagParameter(Prefix='--',Name='hsafe'),\
    
    # --nonbanded Specify that the second stage alignment algorithm be standard,
    #   non-banded, non-D&C CYK. When --nonbanded is enabled, the program fails
    #   with a non-zero exit code and prints an error message if the parsetree
    #   score for any sequence from stage 1 D&C alignment and stage 2 alignment
    #   differs by more than 0.01 bits. In theory, this should never happen as
    #   both algorithms are guaranteed to determine the optimal parsetree. For
    #   larger RNAs (more than 300 residues) if memory is limiting, --nonbanded
    #   should be used in combination with --scoreonly.
    '--nonbanded':FlagParameter(Prefix='--',Name='nonbanded'),\
    
    # --scoreonly With --nonbanded during the second stage standard non-banded
    #   CYK alignment, use the "score only" variant of the algorithm to save
    #   memory, and don't recover a parse tree.
    '--scoreonly':FlagParameter(Prefix='--',Name='scoreonly'),\
    
    # --viterbi Specify that the second stage alignment algorithm be Viterbi to
    #   a CM Plan 9 HMM.
    '--viterbi':FlagParameter(Prefix='--',Name='viterbi'),\
    
    # --search Run all algorithms in scanning mode, not alignment mode.
    '--search':FlagParameter(Prefix='--',Name='search'),\
    
    # --inside With --search Compare the non-banded scanning Inside algorithm to
    #   the HMM banded scanning Inside algorith, instead of using CYK versions.
    '--inside':FlagParameter(Prefix='--',Name='inside'),\
    
    # --forward With --search Compare the scanning Forward scoring algorithm
    #   against CYK.
    '--forward':FlagParameter(Prefix='--',Name='forward'),\
    
    # --taus <n> Specify the first alignment algorithm as non-banded D&C CYK,
    #   and multiple stages of HMM banded CYK alignment. The first HMM banded
    #   alignment will use tau=1E-<x>, which will be the highest value of tau
    #   used. Must be used in combination with --taue.
    '--taus':ValuedParameter(Prefix='--',Name='taus',Delimiter=' '),\
    
    # --taue <n> Specify the first alignment algorithm as non-banded D&C CYK,
    #   and multiple stages of HMM banded CYK alignment. The final HMM banded
    #   alignment will use tau=1E-<x>, which will be the lowest value of tau
    #   used. Must be used in combination with --taus.
    '--taue':ValuedParameter(Prefix='--',Name='taue',Delimiter=' '),\
    
    # --tfile <f> Print the parsetrees for each alignment of each sequence to
    #   file <f>.
    '--tfile':ValuedParameter(Prefix='--',Name='tfile',Delimiter=' '),\

    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmscore"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str

class Cmsearch(CommandLineApplication):
    """cmsearch application controller."""
    _options = {
    
    # -g Turn on the 'glocal' alignment algorithm, local with respect to the
    #   target database, and global with respect to the model. By default, the
    #   local alignment algorithm is used which is local with respect to both
    #   the target sequence and the model.
    '-g':FlagParameter(Prefix='-',Name='g'),\

    # -Z <x> Calculate E-values as if the target database size was <x> megabases
    #   (Mb). Ignore the actual size of the database. This option is only valid
    #   if the CM file has been calibrated. Warning: the predictions for timings
    #   and survival fractions will be calculated as if the database was of size
    #   <x> Mb, which means they will be inaccurate.
    '-Z':ValuedParameter(Prefix='-',Name='Z',Delimiter=' '),\

    # -o <f> Save the high-scoring alignments of hits to a file <f>. The default
    #   is to write them to standard output.
    '-o':ValuedParameter(Prefix='-',Name='o',Delimiter=' '),\

    # -A <f> Save multiple alignment of all significant hits to file <f>.
    '-A':ValuedParameter(Prefix='-',Name='A',Delimiter=' '),\

    # --tblout <f> Create a new output file <f> and print tabular results to
    #   it.
    '--tblout':ValuedParameter(Prefix='--',Name='tblout',Delimiter=' '),\
    
    # --acc Prefer accenssions over names in output
    '--acc':FlagParameter(Prefix='--',Name='acc'),\

    # --noali Don't output alignments (Used to make smaller output)
    '--noali':FlagParameter(Prefix='--',Name='noali'),\

    # --notextw Remove limit on ascii text width per line
    '--notextw':FlagParameter(Prefix='--',Name='notextw'),\

    # --textw <n> Limit ascii text width per line to <n>, where n > 120
    # Default is 120
    '--textw':ValuedParameter(Prefix='--',Name='textw',Delimiter=' '),\

    # -E <x> Set the E-value cutoff for the per-sequence/strand ranked hit list
    #   to <x>, where <x> is a positive real number.
    '-E':ValuedParameter(Prefix='-',Name='E',Delimiter=' '),\
    
    # -T <x> Set the bit score cutoff for the per-sequence ranked hit list to
    #   <x>, where <x> is a positive real number.
    '-T':ValuedParameter(Prefix='-',Name='T',Delimiter=' '),\

    # --incE <x> Consider E values <= <x> as significant. Default is 0.01
    '--incE':ValuedParameter(Prefix='--',Name='incE',Delimiter=' '),\
    
    # --incT <x> Consider scores >= <x> significant
    '--incT':ValuedParameter(Prefix='--',Name='incT',Delimiter=' '),\
    
    # --toponly Only search the top (Watson) strand of the sequences in seqfile.
    #   By default, both strands are searched.
    '--toponly':FlagParameter(Prefix='--',Name='toponly'),\
    
    # --bottomonly Only search the bottom (Crick) strand of the sequences in
    #   seqfile. By default, both strands are searched.
    '--bottomonly':FlagParameter(Prefix='--',Name='bottomonly'),\

    # --cpu <x> Use <x> parallel CPU workers for computation
    '--cpu':ValuedParameter(Prefix='--',Name='cpu',Delimiter=' '),\

    # --tformat <s> Assert that <seqdb> is in format <s>
    '--tformat':ValuedParameter(Prefix='--',Name='tformat',Delimiter=' '),\
    
    # Expert Options
    # --cut_nc Set the bit score cutoff as the NC cutoff value used by Rfam curators
    #   as the noise cutoff score.
    '--cut_nc':FlagParameter(Prefix='--',Name='cut_nc'),\
    
    # --ga Set the bit score cutoff as the GA cutoff value used by Rfam curators
    #   as the gathering threshold.
    '--cut_ga':FlagParameter(Prefix='--',Name='cut_ga'),\
    
    # --cut_tc Set the bit score cutoff as the TC cutoff value used by Rfam curators
    #   as the trusted cutoff.
    '--cut_tc':FlagParameter(Prefix='--',Name='cut_tc'),\

    # --max Turn all heuristic filtering off (Slow).
    '--max':FlagParameter(Prefix='--',Name='max'),\

    # --nohmm Skipp all HMM filtering, use only CM(Slow).
    '--nohmm':FlagParameter(Prefix='--',Name='nohmm'),\

    # --mid skip first two HMM filtering steps (SSV & Vit).
    '--mid':FlagParameter(Prefix='--',Name='mid'),\

    # --default run search space dependent pipeline. Switch for default behavior.
    '--default':FlagParameter(Prefix='--',Name='default'),\

    # --rfam set heuristic filters at R-fam level (Fast).
    '--rfam':FlagParameter(Prefix='--',Name='rfam'),\

    # --hmmonly use the HMM only, ignore the CM.
    '--hmmonly':FlagParameter(Prefix='--',Name='hmmonly'),\

    # --FZ <x> set filteres to dfaults used for search space of size <x> Mb
    '--FZ':ValuedParameter(Prefix='--',Name='FZ',Delimiter=' '),\

    # --Fmid <x> Set P-value threshold for HMM stages to <x>. Default 0.02
    # Used with --mid
    '--Fmid':ValuedParameter(Prefix='--',Name='FZ',Delimiter=' '),\

    # --notrunc Do not allow truncated hits at sequence termini
    '--notrunc':FlagParameter(Prefix='--',Name='notrunc'),\

    # --anytrunc Allow truncated hits anywhere in the sequence
    '--anytrunc':FlagParameter(Prefix='--',Name='anytrunc'),\

    # --nonull3 Turn off NULL3 post hoc aditional null model
    '--nonull3':FlagParameter(Prefix='--',Name='nonull3'),\

    # --mxsize <x> Set max allowed size of alignemnt DP matrices to <x> Mb
    #Default is 128
    '--mxsize':ValuedParameter(Prefix='--',Name='mxsize',Delimiter=' '),\
    
    # --smxsize <x> Set max allowed size of search DP matrices to <x> Mb
    #Default is 128
    '--smxsize':ValuedParameter(Prefix='--',Name='smxsize',Delimiter=' '),\

    # --cyk Use the CYK algorithm for the final round of searching.
    '--cyk':FlagParameter(Prefix='--',Name='cyk'),\

    # --acyk Align hits with the CYK algorithm. Not optimal for accuracy.
    '--acyk':FlagParameter(Prefix='--',Name='cyk'),\

    # --wcx <x> Set expected max hit length to <x> * cm->clen
    '--wcx':ValuedParameter(Prefix='--',Name='wcx',Delimiter=' '),\

    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmsearch"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str
    
    def _tabfile_out_filename(self):
        
        if self.Parameters['--tblout'].isOn():
            tabfile_filename = self._absolute(str(\
                self.Parameters['--tblout'].Value))
        else:
            raise ValueError, 'No tabfile output file specified.'
        return tabfile_filename
    
    def _tempfile_as_multiline_string(self, data):
        """Write a multiline string to a temp file and return the filename.

            data: a multiline string to be written to a file.

           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

        """
        filename = FilePath(self.getTmpFilename(self.TmpDir))
        data_file = open(filename,'w')
        data_file.write(data)
        data_file.close()
        return filename
    
    def _get_result_paths(self,data):
        result = {}
        if self.Parameters['--tblout'].isOn():
            out_name = self._tabfile_out_filename()
            result['SearchResults'] = ResultPath(Path=out_name,IsWritten=True)
        
        return result

class Cmstat(CommandLineApplication):
    """cmstat application controller."""
    _options = {
    
    # -E <x> Print bitscores that correspond to E-value threshold of <x>
    '-E':ValuedParameter(Prefix='-',Name='E',Delimiter=' '),\

    # -P <x> Print bitscores that correspond to E-value threshold of <x>
    '-P':ValuedParameter(Prefix='-',Name='P',Delimiter=' '),\

    # -T <x> Print E-values that correspond to bitscore threshold of <x>
    '-T':ValuedParameter(Prefix='-',Name='T',Delimiter=' '),\

    # -Z <x> Set database size to <x> Mb for E-value calculations
    '-Z':ValuedParameter(Prefix='-',Name='Z',Delimiter=' '),\

    # --key <s> Print statistics for CM with name or accession <s>
    '--key':ValuedParameter(Prefix='--',Name='key',Delimiter=' '),\

    # --cut_ga Print E-values that correspond to GA bit score thresholds
    '--cut_ga':FlagParameter(Prefix='--',Name='cut_ga'),\

    # --cut_nc Print E-values that correspond to NC bit score thresholds
    '--cut_nc':FlagParameter(Prefix='--',Name='cut_nc'),\

    # --cut_tc Print E-values that correspond to TC bit score thresholds
    '--cut_tc':FlagParameter(Prefix='--',Name='cut_tc'),\

    # --hmmonly print fitler HMM bit scores and E-values, not CM derived
    '--hmmonly':FlagParameter(Prefix='--',Name='hmmonly'),\

    # --nohmmonly print CM bit scores an E-values, not HMM derived
    '--nohmmonly':FlagParameter(Prefix='--',Name='nohmmonly'),\
    
    }
    _parameters = {}
    _parameters.update(_options)
    _command = "cmstat"
    _suppress_stderr=True
    
    def getHelp(self):
        """Method that points to the Infernal documentation."""
        
        help_str = \
        """
        See Infernal documentation at:
        http://infernal.janelia.org/
        """
        return help_str

def cmbuild_from_alignment(aln, structure_string, refine=False, calibrate=False, \
    return_alignment=False,params=None):
    """Uses cmbuild to build a CM file given an alignment and structure string.
    
        - aln: an Alignment object or something that can be used to construct
            one.  All sequences must be the same length.
        - structure_string: vienna structure string representing the consensus
            stucture for the sequences in aln.  Must be the same length as the
            alignment.
        - refine: refine the alignment and realign before building the cm.
            (Default=False)
        - calibrate: whether or not to run cmcalibrate automatically on the cm file
            (Default=False)
        - return_alignment: Return (in Stockholm format) alignment file used to
            construct the CM file.  This will either be the original alignment
            and structure string passed in, or the refined alignment if --refine 
            was used. (Default=False)
            - Note.  This will be a string that can either be written to a file
                or parsed.
    """
    aln = Alignment(aln)
    if len(structure_string) != aln.SeqLen:
        raise ValueError, """Structure string is not same length as alignment.  
        Structure string is %s long. Alignment is %s long."""%(len(structure_string),\
        aln.SeqLen)
    else:
        struct_dict = {'SS_cons':structure_string}

    #Make new Cmbuild app instance.
    app = Cmbuild(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    
    #turn on refine flag if True.
    if refine:
        app.Parameters['--refine'].on(get_tmp_filename(app.WorkingDir))
        
    #Get alignment in Stockholm format
    aln_file_string = stockholm_from_alignment(aln,GC_annotation=struct_dict)
    
    #get path to alignment filename
    aln_path = app._input_as_multiline_string(aln_file_string)
    cm_path = aln_path.split('.txt')[0]+'.cm'
    app.Parameters['-n'].on(cm_path)
    
    filepaths = [cm_path,aln_path]
    
    res = app(filepaths)

    #calibrate model if needed
    if calibrate:
        calapp = Cmcalibrate(WorkingDir='/tmp')
        calres = calapp(cm_path)
    #don't clean up if calibrated because need temp file binary data for calibration
    cm_file = ""
    if calibrate:
        cm_file = res['CmFile']
    else:
        cm_file = res['CmFile'].read()
    
    if return_alignment:
        #If alignment was refined, return refined alignment and structure,
        # otherwise return original alignment and structure.
        #return path to file if 
        if refine:
            aln_file_string = res['Refined'].read()
            if not calibrate:
                res.cleanUp()
        return cm_file, aln_file_string
    #Just return cm_file
    elif calibrate:
        return cm_file
    else:
        res.cleanUp()
        return cm_file


def cmbuild_from_file(stockholm_file_path, refine=False,calibrate=False, \
    return_alignment=False, params=None):
    """Uses cmbuild to build a CM file given a stockholm file.
    
        - stockholm_file_path: a path to a stockholm file.  This file should
            contain a multiple sequence alignment formated in Stockholm format. 
            This must contain a sequence structure line:
                #=GC SS_cons <structure string>
        - refine: refine the alignment and realign before building the cm.
            (Default=False)
        - calibrate: whether or not to run cmcalibrate automatically on the cm file
            (Default=False)
        - return_alignment: Return alignment and structure string used to
            construct the CM file.  This will either be the original alignment
            and structure string passed in, or the refined alignment if
            --refine was used. (Default=False)
    """
    #get alignment and structure string from stockholm file.
    info, aln, structure_string = \
        list(MinimalRfamParser(open(stockholm_file_path,'U'),\
            seq_constructor=ChangedSequence))[0]
    
    #call cmbuild_from_alignment.
    res = cmbuild_from_alignment(aln, structure_string, refine=refine, \
        calibrate=calibrate, return_alignment=return_alignment,params=params)
    return res

def calibrate_file(cm_file_path, params=None):
    '''calls cmcalibrate on CM file passed

        -cm_file_path: Path to the CM built by cmbuild
    '''
    if params is None:
        params = {}  

    app = Cmcalibrate(WorkingDir='/tmp', params=params)
    res = app(cm_file_path)
    return res['ExitStatus']


def cmalign_from_alignment(aln, structure_string, seqs, moltype,\
    include_aln=True,refine=False, return_stdout=False,params=None,\
    cmbuild_params=None):
    """Uses cmbuild to build a CM file, then cmalign to build an alignment.
    
        - aln: an Alignment object or something that can be used to construct
            one.  All sequences must be the same length.
        - structure_string: vienna structure string representing the consensus
            stucture for the sequences in aln.  Must be the same length as the
            alignment.
        - seqs: SequenceCollection object or something that can be used to
            construct one, containing unaligned sequences that are to be aligned 
            to the aligned sequences in aln.
        - moltype: Cogent moltype object.  Must be RNA or DNA.
        - include_aln: Boolean to include sequences in aln in final alignment.
            (Default=True)
        - refine: refine the alignment and realign before building the cm.
            (Default=False)
        - return_stdout: Boolean to return standard output from infernal.  This
            includes alignment and structure bit scores and average
            probabilities for each sequence. (Default=False)
    """
    #NOTE: Must degap seqs or Infernal well seg fault!
    seqs = SequenceCollection(seqs,MolType=moltype).degap()
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    
    cm_file, aln_file_string = cmbuild_from_alignment(aln, structure_string,\
        refine=refine,return_alignment=True,params=cmbuild_params)
    
    if params is None:
        params = {}    
    params.update({MOLTYPE_MAP[moltype]:True})
    
    app = Cmalign(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    app.Parameters['--informat'].on('FASTA')
    
    #files to remove that aren't cleaned up by ResultPath object
    to_remove = []    
    #turn on --withali flag if True.
    if include_aln:
        app.Parameters['--withali'].on(\
            app._tempfile_as_multiline_string(aln_file_string))
        #remove this file at end
        to_remove.append(app.Parameters['--withali'].Value)
    
    seqs_path = app._input_as_multiline_string(int_map.toFasta())
    cm_path = app._tempfile_as_multiline_string(cm_file)
    
    #add cm_path to to_remove
    to_remove.append(cm_path)
    paths = [cm_path,seqs_path]

    app.Parameters['-o'].on(get_tmp_filename(app.WorkingDir))
    
    res = app(paths)
    
    info, aligned, struct_string = \
        list(MinimalRfamParser(res['Alignment'].readlines(),\
            seq_constructor=SEQ_CONSTRUCTOR_MAP[moltype]))[0]
    
    #Make new dict mapping original IDs
    new_alignment={}
    for k,v in aligned.NamedSeqs.items():
        new_alignment[int_keys.get(k,k)]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    
    std_out = res['StdOut'].read()
    #clean up files
    res.cleanUp()
    for f in to_remove: remove(f)
    
    if return_stdout:
        return new_alignment, struct_string, std_out
    else:
        return new_alignment, struct_string
    

def cmalign_from_file(cm_file_path, seqs, moltype, alignment_file_path=None,\
    include_aln=False,return_stdout=False,params=None):
    """Uses cmalign to align seqs to alignment in cm_file_path.
        
        - cm_file_path: path to the file created by cmbuild, containing aligned
            sequences. This will be used to align sequences in seqs.
        - seqs: unaligned sequendes that are to be aligned to the sequences in
            cm_file.
        - moltype: cogent.core.moltype object.  Must be DNA or RNA
        - alignment_file_path: path to stockholm alignment file used to create
            cm_file.
            __IMPORTANT__: This MUST be the same file used by cmbuild
            originally.  Only need to pass in this file if include_aln=True.
            This helper function will NOT check if the alignment file is correct
            so you must use it correctly.
        - include_aln: Boolean to include sequences in aln_file in final
            alignment. (Default=False)
        - return_stdout: Boolean to return standard output from infernal.  This
            includes alignment and structure bit scores and average
            probabilities for each sequence. (Default=False)
    """
    #NOTE: Must degap seqs or Infernal well seg fault!
    seqs = SequenceCollection(seqs,MolType=moltype).degap()
    
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    
    if params is None:
        params = {}
    params.update({MOLTYPE_MAP[moltype]:True})
    
    app = Cmalign(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    app.Parameters['--informat'].on('FASTA')
        
    #turn on --withali flag if True.
    if include_aln:
        if alignment_file_path is None:
            raise DataError, """Must have path to alignment file used to build CM if include_aln=True."""
        else:
            app.Parameters['--withali'].on(alignment_file_path)
                
    seqs_path = app._input_as_multiline_string(int_map.toFasta())
    paths = [cm_file_path,seqs_path]
    
    app.Parameters['-o'].on(get_tmp_filename(app.WorkingDir))
    res = app(paths)
    
    info, aligned, struct_string = \
        list(MinimalRfamParser(res['Alignment'].readlines(),\
            seq_constructor=SEQ_CONSTRUCTOR_MAP[moltype]))[0]
    
    
    #Make new dict mapping original IDs
    new_alignment={}
    for k,v in aligned.items():
        new_alignment[int_keys.get(k,k)]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    std_out = res['StdOut'].read()
    res.cleanUp()
    if return_stdout:
        return new_alignment, struct_string, std_out
    else:
        return new_alignment, struct_string
    
def cmsearch_from_alignment(aln, structure_string, seqs, moltype, cutoff=0.0,\
    refine=False,params=None):
    """Uses cmbuild to build a CM file, then cmsearch to find homologs.
    
        - aln: an Alignment object or something that can be used to construct
            one.  All sequences must be the same length.
        - structure_string: vienna structure string representing the consensus
            stucture for the sequences in aln.  Must be the same length as the
            alignment.
        - seqs: SequenceCollection object or something that can be used to
            construct one, containing unaligned sequences that are to be
            searched.
        - moltype: cogent.core.moltype object.  Must be DNA or RNA
        - cutoff: bitscore cutoff.  No sequences < cutoff will be kept in
            search results. (Default=0.0).  Infernal documentation suggests
            a cutoff of log2(number nucleotides searching) will give most
            likely true homologs.
        - refine: refine the alignment and realign before building the cm.
            (Default=False)
    """
    #NOTE: Must degap seqs or Infernal well seg fault!
    seqs = SequenceCollection(seqs,MolType=moltype).degap()
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    
    cm_file, aln_file_string = cmbuild_from_alignment(aln, structure_string,\
        refine=refine,calibrate=True, return_alignment=True)
    
    app = Cmsearch(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    app.Parameters['--tformat'].on('FASTA')
    app.Parameters['-T'].on(cutoff)
    
    to_remove = []
    
    seqs_path = app._input_as_multiline_string(int_map.toFasta())
    cm_path = cm_file.name

    #calibrate CM before search
    calapp = Cmcalibrate(WorkingDir='/tmp')
    calres = calapp(cm_path)

    paths = [cm_path,seqs_path]
    to_remove.append(cm_path)
    
    app.Parameters['--tblout'].on(get_tmp_filename(app.WorkingDir))
    res = app(paths)
    search_results = list(CmsearchParser(res['SearchResults'].readlines()))
    #put headers back 
    if search_results:
        for i,line in enumerate(search_results):
            label = line[0]
            search_results[i][0]=int_keys.get(label,label)
    
    res.cleanUp()
    for f in to_remove:remove(f)
    return search_results

def cmsearch_from_file(cm_file_path, seqs, moltype, cutoff=0.0, params=None):
    """Uses cmbuild to build a CM file, then cmsearch to find homologs.
    
        - cm_file_path: path to the file created by cmbuild, containing aligned
            sequences. This will be used to search sequences in seqs.
        - seqs: SequenceCollection object or something that can be used to
            construct one, containing unaligned sequences that are to be
            searched.
        - moltype: cogent.core.moltype object.  Must be DNA or RNA
        - cutoff: bitscore cutoff.  No sequences < cutoff will be kept in
            search results. (Default=0.0).  Infernal documentation suggests
            a cutoff of log2(number nucleotides searching) will give most
            likely true homologs.
    """
    #NOTE: Must degap seqs or Infernal well seg fault!
    seqs = SequenceCollection(seqs,MolType=moltype).degap()

    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seqs.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=RNA)
    
    app = Cmsearch(InputHandler='_input_as_paths',WorkingDir='/tmp',\
        params=params)
    app.Parameters['--tformat'].on('FASTA')
    app.Parameters['-T'].on(cutoff)
    
    seqs_path = app._input_as_multiline_string(int_map.toFasta())

    paths = [cm_file_path,seqs_path]
    
    app.Parameters['--tblout'].on(get_tmp_filename(app.WorkingDir))
    res = app(paths)
    
    search_results = list(CmsearchParser(res['SearchResults'].readlines()))
    
    if search_results:
        for i,line in enumerate(search_results):
            label = line[0]
            search_results[i][0]=int_keys.get(label,label)

    res.cleanUp()

    return search_results

