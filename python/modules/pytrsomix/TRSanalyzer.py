from Bio import pairwise2
from Bio.Seq import Seq
import pandas as pd
import parasail
from enum import Enum

class TRS_cols(Enum):
    SEQ_COLUMN = ">SEQ"
    GENOME_COLUMN = "GENOME"

class SeqAnalyzer():
    def __init__(self, seqs: list):
        self.seqs = seqs
        if len(seqs) > 0:
            self.seqs_combined = pd.concat(self.seqs, axis=0).reset_index(drop=True)
        else:
            self.seqs_combined = pd.DataFrame()

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

    def calculate_all_alignments_nw(self, idx):
        objective_seq = self.seqs_combined.loc[idx, TRS_cols.SEQ_COLUMN.value]
        remaining_seq = self.seqs_combined[TRS_cols.SEQ_COLUMN.value]
        algns = {}
        for seq, idx in zip(remaining_seq, remaining_seq.index):
            a = parasail.nw_scan_16(objective_seq, seq, 1, 1, parasail.blosum75)
            algns[idx] = a
        return algns
    
    def calculate_all_alignments_sw(self, idx):
        objective_seq = self.seqs_combined.loc[idx, TRS_cols.SEQ_COLUMN.value]
        remaining_seq = self.seqs_combined[TRS_cols.SEQ_COLUMN.value]
        algns = {}
        for seq, idx in zip(remaining_seq, remaining_seq.index):
            a = parasail.sw_scan_16(objective_seq, seq, 1, 1, parasail.blosum75)
            algns[idx] = a
        return algns

    def calculate_single_alignment_nw(self, idx, jdx):
        seq_idx = self.seqs_combined.loc[idx, TRS_cols.SEQ_COLUMN.value]
        seq_jdx = self.seqs_combined.loc[jdx, TRS_cols.SEQ_COLUMN.value]
        a = parasail.nw_scan_16(seq_idx, seq_jdx, 1, 1, parasail.blosum75)
        return a
    
    def calculate_single_alignment_sw(self, idx, jdx):
        seq_idx = self.seqs_combined.loc[idx, TRS_cols.SEQ_COLUMN.value]
        seq_jdx = self.seqs_combined.loc[jdx, TRS_cols.SEQ_COLUMN.value]
        a = parasail.sw_scan_16(seq_idx, seq_jdx, 1, 1, parasail.blosum75)
        return a
    
    def calculate_single_alignment_biopython(self, idx, jdx):
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
    
    @Combined.setter
    def Combined(self, value):
        self.seqs_combined = value

class AlignmentAnalyzer():
    def __init__(self, algns: dict):
        self.algns = algns

    def get_sorted_scores(self):
        return pd.DataFrame([(idx, val.score) for idx, val in self.algns.items()], columns=["index", "score"]).set_index("index")