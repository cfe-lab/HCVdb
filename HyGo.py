import HyPhy
import os

class HyGo():
    def __init__(self, cwd=os.getcwd(), nthreads=1, alphabet='ACGT',
                 matrix='{{5,-4,-4,-4},{-4,5,-4,-4},{-4,-4,5,-4},{-4,-4,-4,5}}',
                 gap_open=20, gap_open2=20, gap_extend=10, gap_extend2=10,
                 no_terminal_penalty=1):
        self.__instance = HyPhy._THyPhy(cwd, nthreads)

        self.call('alignOptions = {};')
        self.call('alignOptions["SEQ_ALIGN_CHARACTER_MAP"]="%s";' % (alphabet,))
        self.call('alignOptions["SEQ_ALIGN_SCORE_MATRIX"]=%s;' % (matrix,))
        self.call('alignOptions["SEQ_ALIGN_GAP_OPEN"]=%d;' % (gap_open,))
        self.call('alignOptions["SEQ_ALIGN_GAP_OPEN2"]=%d;' % (gap_open2,))
        self.call('alignOptions["SEQ_ALIGN_AFFINE"]=1;')
        self.call('alignOptions["SEQ_ALIGN_GAP_EXTEND"]=%d;' % (gap_extend,))
        self.call('alignOptions["SEQ_ALIGN_GAP_EXTEND2"]=%d;' % (gap_extend2,))
        self.call('alignOptions["SEQ_ALIGN_NO_TP"]=%d;' % (no_terminal_penalty,))


    def _get_version(self):
        self.call('GetString(v, HYPHY_VERSION, 1);')
        return self.get('v', as_string=True)

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

