import cffi
import numpy as np
import pandas as pd
from io import StringIO

with open("TRS-omix.c", "r") as f:
    lines = f.read()

ffi = cffi.FFI()
ffi.set_source("_trs_omix_python", lines)
ffi.cdef(""" 
void PrintError(int e);
typedef struct {...;} NPt;
typedef struct {...;} NLt;
long long int ImportGenome(char *fn, NLt *nlt);
char * PerformTRSCalculation(char *gfn, char* tfn, char* ifn, long tmin, long tmax, int mode);
""")

ffi.compile(verbose=0)