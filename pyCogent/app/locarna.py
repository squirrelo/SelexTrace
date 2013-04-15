#!/usr/bin/env python
# cogent.app.locarna
"""
Provides an application controller for the commandline version of:
LocARNA 1.7.5
"""
import shutil
from cogent.app.parameters import FlagParameter, ValuedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, get_tmp_filename
from cogent.parse.clustal import ClustalParser
from cogent.core.moltype import DNA, RNA
from cogent.core.alignment import SequenceCollection, Alignment, DataError
from os import remove

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Jeremy Widmann, Joshua Shorenstein"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Joshua Shorenstein"
__email__ = "joshua.shorenstein@colorado.edu"
__status__ = "Development"

               
class LocARNA(CommandLineApplication):
    """LocARNA application controller."""
    _options = {
    #Options:

    # === Scoring parameters ===
    # -m,--match=<score>(50)           Match score
    '-m':ValuedParameter(Prefix='-',Name='m',Delimiter='='),\
    
    # -M,--mismatch=<score>(0)         Mismatch score
    '-M':ValuedParameter(Prefix='-',Name='M',Delimiter='='),\
    
    # --ribosum-file=<f>(RIBOSUM85_60) Ribosum file
    '--ribosum-file':ValuedParameter(Prefix='--',Name='ribosum-file',\
        Delimiter='='),\
    
    # --use-ribosum=<bool>(true)       Use ribosum scores
    '--use-ribosum':ValuedParameter(Prefix='--',Name='use-ribosum',\
        Delimiter='='),\
    
    # -i,--indel=<score>(-350)         Indel score
    '-i':ValuedParameter(Prefix='-',Name='i',Delimiter='='),\
    
    # --indel-opening=<score>(-500)    Indel opening score
    '--indel-opening':ValuedParameter(Prefix='--',Name='indel-opening',\
        Delimiter='='),\
    
    # -s,--struct-weight=<score>(180)  Maximal weight of 1/2 arc match
    '-s':ValuedParameter(Prefix='-',Name='s',Delimiter='='),\
    
    # -e,--exp-prob=<prob>             Expected probability
    '-e':ValuedParameter(Prefix='-',Name='e',Delimiter='='),\
    
    # -t,--tau=<factor>(0)             Tau factor in percent
    '-t':ValuedParameter(Prefix='-',Name='t',Delimiter='='),\
    
    # -E,--exclusion=<score>(0)        Exclusion weight
    '-E':ValuedParameter(Prefix='-',Name='E',Delimiter='='),\
    
    # --stacking         Use stacking terms (needs stack-probs by RNAfold -p2)
    '--stacking':FlagParameter(Prefix='--',Name='stacking'),\

    # === Type of locality ===
    # --struct-local=<bool>(false)     Structure local
    # NOTE!  Must set to 0 for False and 1 for True
    '--struct-local':ValuedParameter(Prefix='--',Name='struct-local',\
        Delimiter='='),\
    
    # --sequ-local=<bool>(false)       Sequence local
    '--sequ-local':ValuedParameter(Prefix='--',Name='sequ-local',Delimiter='='),\

    # === Controlling output ===
    # -w,--width=<columns>(120)        Output width
    '-w':ValuedParameter(Prefix='-',Name='w',Delimiter='='),\
    
    # --clustal=<file>                 Clustal output
    '--clustal':ValuedParameter(Prefix='--',Name='clustal',Delimiter='='),\
    
    # --pp=<file>                      PP output
    '--pp':ValuedParameter(Prefix='--',Name='pp',Delimiter='='),\
    
    # -L,--local-output                Output only local sub-alignment
    '-L':FlagParameter(Prefix='-',Name='L'),\
    
    # -P,--pos-output                  Output only local sub-alignment positions
    '-P':FlagParameter(Prefix='-',Name='P'),\
    
    # --write-structure                Write guidance structure in output
    '--write-structure':FlagParameter(Prefix='-',Name='write-structure'),\

    # === Heuristics for speed accuracy trade off ===
    # -p,--min-prob=<prob>(0.0005)     Minimal probability
    '-p':ValuedParameter(Prefix='-',Name='p',Delimiter='='),\
    
    # -D,--max-diff-am=<diff>(-1)      Maximal difference for sizes of matched 
    # arcs
    '-D':ValuedParameter(Prefix='-',Name='D',Delimiter='='),\
    
    # -a,--min-am-prob=<amprob>(0.0005) Minimal Arc-match probability
    '-a':ValuedParameter(Prefix='-',Name='a',Delimiter='='),\
    
    # -b,--min-bm-prob=<bmprob>(0.0005) Minimal Base-match probability
    '-b':ValuedParameter(Prefix='-',Name='b',Delimiter='='),\

    # === Special sauce options ===
    # -k,--kbest=<k>(-1)               Find k-best alignments
    '-k':ValuedParameter(Prefix='-',Name='k',Delimiter='='),\

    # === Options for controlling MEA score (experimental, under construction)
    # --mea-alignment                  Do MEA alignment
    '--mea-alignment':FlagParameter(Prefix='-',Name='mea-alignment'),\
    
    # --probcons-file=<file>           Probcons parameter file
    '--probcons-file':ValuedParameter(Prefix='--',Name='probcons-file',\
        Delimiter='='),\
    
    # --match-prob-method=<int>(0)     Method for computation of match probs
    '--match-prob-method':ValuedParameter(Prefix='--',Name='match-prob-method',\
        Delimiter='='),\
        
    # --temperature=<int>(300)         Temperature for PF-computation
    '--temperature':ValuedParameter(Prefix='--',Name='temperature',\
        Delimiter='='),\
        
    # --pf-struct-weight=<weight>(200) Structure weight in PF-computation
    '--pf-struct-weight':ValuedParameter(Prefix='--',Name='pf-struct-weight',\
        Delimiter='='),\
    
    # --mea-alpha=<weight>(100)        Weight alpha for MEA
    '--mea-alpha':ValuedParameter(Prefix='--',Name='mea-alpha',Delimiter='='),\
    
    # --mea-beta=<weight>(100)         Weight beta for MEA
    '--mea-beta':ValuedParameter(Prefix='--',Name='mea-beta',Delimiter='='),\
    
    # --mea-gamma=<weight>(100)        Weight gamma for MEA
    '--mea-gamma':ValuedParameter(Prefix='--',Name='mea-gamma',Delimiter='='),\
    
    # --write-match-probs=<file>       Write match probs to file (don't align!)
    '--write-match-probs':ValuedParameter(Prefix='--',Name='write-match-probs',\
        Delimiter='='),\
        
    # --read-match-probs=<file>        Read match probabilities from file
    '--read-match-probs':ValuedParameter(Prefix='--',Name='read-match-probs',\
        Delimiter='='),\
        
    # --write-arcmatch-scores=<file>   Write arcmatch scores (don't align!)
    '--write-arcmatch-scores':ValuedParameter(Prefix='--',\
        Name='write-arcmatch-scores',Delimiter='='),\
        
    # --read-arcmatch-scores=<file>    Read arcmatch scores
    '--read-arcmatch-scores':ValuedParameter(Prefix='--',\
        Name='read-arcmatch-scores',Delimiter='='),\

    # === Constraints ===
    # --noLP                           No lonely pairs
    '--noLP':FlagParameter(Prefix='--',Name='noLP'),\
    
    # --anchorA=<string>()             Anchor constraints sequence A
    '--anchorA':ValuedParameter(Prefix='--',Name='anchorA',Delimiter='='),\
    
    # --anchorB=<string>()             Anchor constraints sequence B
    '--anchorB':ValuedParameter(Prefix='--',Name='anchorB',Delimiter='='),\
    
    # --ignore-constraints             Ignore constraints in pp-file
    '--ignore-constraints':FlagParameter(Prefix='--',\
        Name='ignore-constraints'),\
    
    # Not Supported:
    # -h,--help                        This help
    # -V,--version                     Version info
    # -v,--verbose                     Verbose
    }
    
    _synonyms = {
                '--match':'-m',
                '--mismatch':'-M',
                '--indel':'-i',
                '--struct-weight':'-s',
                '--exp-prob':'-e',
                '--tau':'-t',
                '--exclusion':'-E',
                '--width':'-w',
                '--local-output':'-L',
                '--pos-output':'-P',
                '--min-prob':'-p',
                '--max-diff-am':'-D',
                '--min-am-prob':'-a',
                '--min-bm-prob':'-b',
                '--kbest':'-k',
                }
    
    _parameters = {}
    _parameters.update(_options)
    _command = "locarna"
    _suppress_stderr=True

    def getHelp(self):
        """Method that points to the LocARNA documentation."""
        
        help_str = \
        """
        See LocARNA documentation at:
        http://www.bioinf.uni-freiburg.de/Software/LocARNA/
        """
        return help_str

    def _get_result_paths(self,data):
        raise NotImplementedError

class MLocarna(LocARNA):
    """MLocarna application controller.
    
        NOTE: Documentation is sparse! Man page says (regarding options):
            "as well as most of the options of locarna (see man page of
                locarna)"
            About MLocarna having the same options as locarna.  This is not true
            for all options.  Need to explicitly test all to see which, but
            don't feel like it now.  Use at own risk until this is determined.
    """
    
    _options = {
        #Options:

        #--probabilistic
        '--probabilistic':FlagParameter(Prefix='--',Name='probabilistic'),\

        #--consistency-transformation
        '--consistency-transformation':FlagParameter(Prefix='--',Name='consistency-transformation'),\

        #--cpus
        '--cpus':ValuedParameter(Prefix='--',Name='cpus',Delimiter='='),\

        # --write-structure                Write guidance structure in output
        '--write-structure':FlagParameter(Prefix='-',Name='write-structure'),\
        
        #-P, --tree-min-prob=<f>       Minimal prob for constructing guide tree
        '-P':ValuedParameter(Prefix='-',Name='P',Delimiter='='),\
        
        # --noLP                        No lonely pair
        '--noLP':FlagParameter(Prefix='--',Name='noLP'),\
        
        # --LP                          Allow lonely pairs
        '--LP':FlagParameter(Prefix='--',Name='LP'),\
        
        #--treefile                    Guide tree file
        '--treefile':ValuedParameter(Prefix='--',Name='treefile',\
            Delimiter='='),\

        #--tgtdir                      Target directory
        '--tgtdir':ValuedParameter(Prefix='--',Name='tgtdir',Delimiter='='),\
        
        # --skip-pp                     Skip computation of pair probs if
        # already existing
        '--skip-pp':FlagParameter(Prefix='--',Name='skip-pp'),\
        
        # --extlib                      Use library extension -- experimental
        # feature
        '--extlib':FlagParameter(Prefix='--',Name='extlib'),\
        
        #--iterate                     Refine iteratively after progressive
        # alignment (experimental)
        '--iterate':FlagParameter(Prefix='--',Name='iterate'),\
        
        #--iterations                  Refine iteratively for given number
        # of iterations (experimental)
        '--iterations':ValuedParameter(Prefix='--',Name='iterations',\
            Delimiter='='),\
        
        # --ignore-constraints          Ignore constraints even if given
        '--ignore-constraints':FlagParameter(Prefix='--',\
            Name='ignore-constraints'),\
        
        #LocARNA options that work with MLocARNA
        
        # === Type of locality ===
        # --struct-local=<bool>(false)     Structure local
        # NOTE!  Must set to 0 for False and 1 for True
        '--struct-local':ValuedParameter(Prefix='--',Name='struct-local',\
            Delimiter='='),\
        
        # --sequ-local=<bool>(false)       Sequence local
        '--sequ-local':ValuedParameter(Prefix='--',Name='sequ-local',\
            Delimiter='='),\

    }
    
    _synonyms = {
                '--match':'-m',
                '--mismatch':'-M',
                '--indel':'-i',
                '--struct-weight':'-s',
                '--exp-prob':'-e',
                '--tau':'-t',
                '--exclusion':'-E',
                '--width':'-w',
                '--local-output':'-L',
                '--tree-min-prob':'-P',
                '--min-prob':'-p',
                '--max-diff-am':'-D',
                '--min-am-prob':'-a',
                '--min-bm-prob':'-b',
                '--kbest':'-k',
                }
    
    _parameters = {}
    _parameters.update(_options)
    _command = "mlocarna"
    _suppress_stderr=True
    
    def _get_result_paths(self,data):
        result = {}
        if self.Parameters['--tgtdir'].isOn():
            out_path = self._absolute(str(\
                self.Parameters['--tgtdir'].Value))
        else:
            out_path = self._input_filename+'.out'
            
        result['ClustalAlignment'] = \
            ResultPath(Path=out_path+'/results/result.aln',IsWritten=True)
        #result['Structure'] = \
        #    ResultPath(Path=out_path+'/alifold.out',IsWritten=True)
        if  self.Parameters['--probabilistic'].isOn():
            result['ProbabilisticAlignment'] = \
            ResultPath(Path=out_path+'/results/result_prog.aln',IsWritten=True)
        return result

def create_locarnap_alignment(seqs,moltype,struct=False,params=None):
    """Returns mlocarna results given an unaligned SequenceCollection.
    
        - seqs: A SequenceCollection object or something that behaves like one.
        - moltype: cogent.core.moltype object.
        -struct: Boolean whether or not to also output vienna structure string
    """
    #Construct SequenceCollection object.
    seqs = SequenceCollection(seqs,MolType=moltype)

    #need to make int map.
    int_map, int_keys = seqs.getIntMap()
    #construct SequenceCollection object from int map to use functionality
    int_map = SequenceCollection(int_map, MolType=moltype)
    
    #Create application.
    app = MLocarna(InputHandler='_input_as_multiline_string',params=params)
    #Get temporary directory to write all mlocarna files.
    mlocarna_dir = get_tmp_filename(suffix='')
    app.Parameters['--tgtdir'].on(mlocarna_dir)

    #set parameters to run locarna-p
    app.Parameters['--write-structure'].on()
    app.Parameters['--probabilistic'].on()
    app.Parameters['--consistency-transformation'].on()
    res = app(int_map.toFasta())
    #get the structure from the results if necessary
    if struct:
        structfile = open(res['ProbabilisticAlignment'].name, 'U')
        structure = ""
        newstrline = True
        for line in structfile:
            line = line.strip()
            #read in structure lines of alignment (--write-structure)
            if len(line) > 0 and (line[0] == "." or line[0] == "("):
                #only append if new structure aspect, since struct is 
                #written both above and below blocks in alignment
                if newstrline:
                    structure += line
                    newstrline = not newstrline
                else:
                    newstrline = not newstrline
        


    aligned = dict(ClustalParser(res['ClustalAlignment']))
    
    #Make new dict mapping original IDs
    new_alignment={}
    for k,v in aligned.items():
        new_alignment[int_keys.get(k,k)]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)

    #Clean up after MlocARNA
    res.cleanUp()
    shutil.rmtree(mlocarna_dir)

    #output alignment and structure if asked for, else outout just alignment
    if struct:
        return new_alignment, structure
    else:
        return new_alignment