import cffi

# #with open("src/TRS-omix.c", "r") as f:
# #    lines = f.read()

# ffibuilder = cffi.FFI()
# ffibuilder.set_source("_trs_omix", "", sources=["src/TRS-omix.c"])
# # ffibuilder.set_source("_trs_omix_python", "", include_dirs=["src"], sources=["src/TRS-omix.c"])
# ffibuilder.cdef("""
# char * PerformTRSCalculation(char *gfn, char* tfn, char* ifn, long tmin, long tmax, int mode);
# """)


with open("src/TRS-omix.c", "r") as f:
    lines = f.read()

ffibuilder = cffi.FFI()
ffibuilder.set_source("_trs_omix", lines)
ffibuilder.cdef(""" 
char * PerformTRSCalculation(char *gfn, char* tfn, char* ifn, long tmin, long tmax, int mode);
""")

ffibuilder.compile(verbose=0)



if __name__ == "__main__":    # not when running with setuptools
    ffibuilder.compile(verbose=True)