import HyPhy
import os
import re

gap_prefix = re.compile('^[-]+')
gap_suffix = re.compile('[-]+$')

def get_boundaries (str):
    # return a tuple giving indices of subsequence without gap prefix and suffix
    res = [0,len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])
    if right:
        res[1] = len(str) - len(right[0])

    return res


nuc_alphabet = 'ACGT'
nuc_score_matrix = '{{5,-4,-4,-4},{-4,5,-4,-4},{-4,-4,5,-4},{-4,-4,-4,5}}'

amino_alphabet = 'ARNDCQEGHILKMFPSTWYV'
blosum62 = """\
{{6,-2,-2,-3,-1,-1,-1,0,-2,-2,-2,-1,-1,-3,-1,2,0,-4,-3,0},\
{-2,8,-1,-2,-5,1,0,-3,0,-5,-3,3,-2,-4,-3,-1,-2,-4,-3,-4},\
{-2,-1,9,2,-4,0,0,-1,1,-5,-5,0,-3,-5,-3,1,0,-6,-3,-4},\
{-3,-2,2,9,-5,0,2,-2,-2,-5,-5,-1,-5,-5,-2,0,-2,-6,-5,-5},\
{-1,-5,-4,-5,13,-4,-5,-4,-4,-2,-2,-5,-2,-4,-4,-1,-1,-3,-4,-1},\
{-1,1,0,0,-4,8,3,-3,1,-4,-3,2,-1,-5,-2,0,-1,-3,-2,-3},\
{-1,0,0,2,-5,3,7,-3,0,-5,-4,1,-3,-5,-2,0,-1,-4,-3,-4},\
{0,-3,-1,-2,-4,-3,-3,8,-3,-6,-5,-2,-4,-5,-3,0,-2,-4,-5,-5},\
{-2,0,1,-2,-4,1,0,-3,11,-5,-4,-1,-2,-2,-3,-1,-3,-4,3,-5},\
{-2,-5,-5,-5,-2,-4,-5,-6,-5,6,2,-4,2,0,-4,-4,-1,-4,-2,4},\
{-2,-3,-5,-5,-2,-3,-4,-5,-4,2,6,-4,3,1,-4,-4,-2,-2,-2,1},\
{-1,3,0,-1,-5,2,1,-2,-1,-4,-4,7,-2,-5,-2,0,-1,-4,-3,-3},\
{-1,-2,-3,-5,-2,-1,-3,-4,-2,2,3,-2,8,0,-4,-2,-1,-2,-1,1},\
{-3,-4,-5,-5,-4,-5,-5,-5,-2,0,1,-5,0,9,-5,-4,-3,1,4,-1},\
{-1,-3,-3,-2,-4,-2,-2,-3,-3,-4,-4,-2,-4,-5,11,-1,-2,-6,-4,-4},\
{2,-1,1,0,-1,0,0,0,-1,-4,-4,0,-2,-4,-1,6,2,-4,-3,-2},\
{0,-2,0,-2,-1,-1,-1,-2,-3,-1,-2,-1,-1,-3,-2,2,7,-4,-2,0},\
{-4,-4,-6,-6,-3,-3,-4,-4,-4,-4,-2,-4,-2,1,-6,-4,-4,16,3,-4},\
{-3,-3,-3,-5,-4,-2,-3,-5,3,-2,-2,-3,-1,4,-4,-3,-2,3,10,-2},\
{0,-4,-4,-5,-1,-3,-4,-5,-5,4,1,-3,1,-1,-4,-2,0,-4,-2,6}};
"""


class HyGo():
    def __init__(self, cwd=os.getcwd(), nthreads=1, alphabet='ACGT',
                 gap_open=20, gap_open2=20, gap_extend=10, gap_extend2=10,
                 no_terminal_penalty=1):
        self.__instance = HyPhy._THyPhy(cwd, nthreads)
        self.call('alignOptions = {};')

        # default settings
        self.set_alphabet()
        self.set_matrix()
        self.set_gap_open()
        self.set_gap_open(20, is_first=False)
        self.set_affine()
        self.set_gap_extend()
        self.set_gap_extend(10, is_first=False)
        self.set_terminal()


    def _get_version(self):
        self.call('GetString(v, HYPHY_VERSION, 1);')
        return self.get('v', as_string=True)

    def set_alphabet(self, alphabet='ACGT'):
        self.call('alignOptions["SEQ_ALIGN_CHARACTER_MAP"]="%s";' % (alphabet,))

    def set_matrix(self, matrix=nuc_score_matrix):
        self.call('alignOptions["SEQ_ALIGN_SCORE_MATRIX"]=%s;' % (matrix,))

    def set_gap_open(self, gap_open=20, is_first=True):
        self.call('alignOptions["SEQ_ALIGN_GAP_OPEN%s"]=%d;' % (
            '' if is_first else '2',
            gap_open
        ))

    def set_affine(self, affine=1):
        self.call('alignOptions["SEQ_ALIGN_AFFINE"]=%d;' % (affine,))

    def set_gap_extend(self, gap_extend=10, is_first=True):
        self.call('alignOptions["SEQ_ALIGN_GAP_EXTEND%s"]=%d;' % (
            '' if is_first else '2',
            gap_extend
        ))

    def set_terminal(self, use_terminal_penalty=False):
        self.call('alignOptions["SEQ_ALIGN_NO_TP"]=%d;' % (
            0 if use_terminal_penalty else 1
        ))

    def call(self, msg, flush=False):
        """
        Executes batch language
        :param msg: HBL command
        :param flush: if False, preserves execution state of the system
        :return: any stdout from HBL execution
        """
        p = self.__instance.ExecuteBF(msg, False)
        return p.sData

    def get(self, expr, as_string=False):
        """
        Retrieve parameter value by name
        :param expr: HBL expression
        :return:
        """
        val = None
        self.call('val=%s;' % (expr,))
        self.call('fprintf(stdout, val);')
        if as_string:
            exec 'val="""%s"""' % (self.stdout(),)
        else:
            exec 'val=%s' % (self.stdout(),)  #FIXME: only works for scalars
        return val

    def stdout(self):
        return self.__instance.GetStdout().sData

    def stderr(self):
        return self.__instance.GetErrors().sData

    def warnings(self):
        return self.__instance.GetWarnings().sData

    def align(self, ref, query):
        """
        Pairwise Gotoh alignment
        :param ref:
        :param query:
        :return:
        """
        self.call('inStr={{"%s","%s"}};' % (ref, query))
        self.call('AlignSequences(aligned, inStr, alignOptions);')

        aligned = None
        exec 'aligned = ' + self.call('return aligned;', False)
        aligned = aligned['0']

        return {
            'score': aligned.get('0', None),
            'ref': aligned.get('1', None),
            'query': aligned.get('2', None)
        }

