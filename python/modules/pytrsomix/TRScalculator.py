from Bio import SeqIO
from _trs_omix import lib, ffi
from io import StringIO
import pandas as pd
import parasail
from enum import Enum


class TRS_cols(Enum):
    SEQ_COLUMN = ">SEQ"
    GENOME_COLUMN = "GENOME"


class TRScalculator:
    def __init__(self, sequence=b"sequence.fasta", trs=b"trs.txt", interiors=b"interiors.txt", tmin=100, tmax=3000, mode=0, seq_id=0):
        self.sequence = sequence
        self.trs = trs
        self.interiors = interiors
        self.tmin = tmin
        self.tmax = tmax
        self.mode = mode
        self.seq_id = seq_id
        self.result = None

    def get_fasta_record_id(self):
        pass
        
    
    def calculate(self):

        fasta_file_content = list(SeqIO.parse(self.sequence, format="fasta"))
        if len(fasta_file_content) != 1:
            exit(-1)

        pGfn = ffi.cast("char *", ffi.from_buffer(self.sequence))
        pTfn = ffi.cast("char *", ffi.from_buffer(self.trs))
        pIfn = ffi.cast("char *", ffi.from_buffer(self.interiors))

        tmin = ffi.cast("long long", self.tmin)
        tmax = ffi.cast("long long", self.tmax)
        mode = ffi.cast("int", self.mode)

        self.result = lib.PerformTRSCalculation(pGfn, pTfn, pIfn, tmin, tmax, mode)
        self.result = ffi.string(self.result).decode("utf-8")
        self.result = StringIO(self.result)
        self.result = pd.read_csv(self.result, sep=";")

        genom_name = fasta_file_content[0].name
        self.result[TRS_cols.GENOME_COLUMN.value] = genom_name
        
    @property
    def Result(self):
        return self.result

from Bio import pairwise2
from Bio.Seq import Seq
class SeqAnalyzer():
    def __init__(self, seqs: list):
        self.seqs = seqs
        self.seqs_combined = pd.concat(self.seqs, axis=0).reset_index()

    def calculate_all_alignments_biopython(self, idx):
        objective_seq = Seq(self.seqs_combined.loc[idx, TRS_cols.SEQ_COLUMN.value])
        remaining_seq = self.seqs_combined.drop([idx], axis=0)[TRS_cols.SEQ_COLUMN.value]
        algns = []
        for seq in remaining_seq[:100]:
            print(seq)
            seq_ = Seq(seq)
            a = pairwise2.align.globalxx(objective_seq, seq_)
            algns.append(a)
        return algns

    def calculate_all_alignments_parasail(self, idx):
        objective_seq = self.seqs_combined.loc[idx, TRS_cols.SEQ_COLUMN.value]
        remaining_seq = self.seqs_combined.drop([idx], axis=0)[TRS_cols.SEQ_COLUMN.value]
        algns = []
        for seq in remaining_seq:
            a = parasail.nw_scan_16(objective_seq, seq, 1, 1, parasail.blosum75)
            algns.append(a)
        return algns
    
    def calculate_single_alignment(self, idx, jdx):
        seq1 = Seq(self.seqs_combined.loc[idx, TRS_cols.SEQ_COLUMN.value])
        seq2 = Seq(self.seqs_combined.loc[jdx, TRS_cols.SEQ_COLUMN.value])
        alignment = pairwise2.align.globalxx(seq1, seq2)
        return alignment

    def get_best_score(self, alignments):
        scores = [item[2] for item in alignments]
        return max(scores)

    @property
    def Nseq(self):
        return self.seqs_combined.shape[0]

    @property
    def Combined(self):
        return self.seqs_combined 
