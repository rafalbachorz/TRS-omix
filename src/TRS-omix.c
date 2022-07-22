#include <stdio.h>
#include <stdlib.h>
/*
+---------------------------------------------------------------------------+
| NPt – data structure serving to store the position of s single nucleotide |
| in the genome and the sign representing this nucleotide.                  |
+---------------------------------------------------------------------------+
*/
typedef struct _NPt {

    long long no;
    char sgn;

    struct _NPt *next;
    struct _NPt *prev;

} NPt;
/*
+---------------------------------------------------------------------------+
| NLt – data structure serving to store the data on the genome with the     |
| use of data structure NPt. in the genome and the sign representing this   |
| nucleotide.                                                               |
+---------------------------------------------------------------------------+
*/
typedef struct _NLt {

    NPt *head;
    NPt *tail;

} NLt;
/*
+---------------------------------------------------------------------------+
| MPt – data structure serving to store the number of the TRS motifs, the   |
| number of the class of TRS motifs and a single sign representing the      |
| nucleotide of the TRS motif from the file trs.txt.                        |
+---------------------------------------------------------------------------+
*/
typedef struct _MPt {

    int no;
    int cl;
    char n;

    struct _MPt *r;
    struct _MPt *b;

} MPt;
/*
+---------------------------------------------------------------------------+
| MLt – data structure serving to store the TRS motifs contained in the file|
| trs.txt with the use of data structure MPt.                               |
+---------------------------------------------------------------------------+
*/
typedef struct _MLt {

    MPt *rt;
    MPt *r;
    MPt *b;

} MLt;
/*
+---------------------------------------------------------------------------+
| VPt – data structure serving to store the left and the right positions    |
| of the flanking sequence which are found, as well as TRS constituting the |
| base for the flanking sequence.                                           |
+---------------------------------------------------------------------------+
*/
typedef struct _VPt {

    int inf;

    int no;
    int cl;

    long long int ps;
    long long int pe;

    NLt *nlt;

    struct _VPt *next;
    struct _VPt *prev;

} VPt;
/*
+---------------------------------------------------------------------------+
| VLt – data structure serving to store the position of the flanking        |
| sequence.                                                                 |
+---------------------------------------------------------------------------+
*/
typedef struct _VLt {

    VPt *head;
    VPt *tail;

} VLt;
/*
+---------------------------------------------------------------------------+
| InitNLt() – this function initializes the linked list of NLt type.        |
+---------------------------------------------------------------------------+
*/
void InitNLt(NLt *nlt) {

    nlt->head = NULL;
    nlt->tail = NULL;
}
/*
+---------------------------------------------------------------------------+
| InitMLt() – this function initializes the linked list of MLt type.        |
+---------------------------------------------------------------------------+
*/
void InitMLt(MLt *l) {

    l->rt = NULL;
    l->r =  NULL;
    l->b =  NULL;
}
/*
+---------------------------------------------------------------------------+
| InitVLt() – this function initializes the linked list of VLt type.        |
+---------------------------------------------------------------------------+
*/
void InitVLt(VLt *vlt) {

    vlt->head = NULL;
    vlt->tail = NULL;
}
/*
+---------------------------------------------------------------------------+
| FreeNLt() – this function frees memory assigned to the linked list of NLt |
| type.                                                                     |
+---------------------------------------------------------------------------+
*/
void FreeNLt(NLt *nlt) {

    if (nlt->head != NULL) {

        NPt *ctnpt = nlt->head, *tmp;

        while (ctnpt != NULL) {
            tmp = ctnpt->next;
            free(ctnpt);
            ctnpt = tmp;
        }
    }
}
/*
+---------------------------------------------------------------------------+
| void FreeMLt() – this function frees the memory assigned to the linked    |
| list of MLt type.                                                         |
+---------------------------------------------------------------------------+
*/
void FreeMLt(MLt *mlt) {

    if (mlt->rt != NULL) {

        MPt *ctmpt1 = mlt->rt, *ctmpt2, *tmp;

        while (ctmpt1 != NULL) {

            ctmpt2 = ctmpt1;

            while (ctmpt2 != NULL) {
                tmp = ctmpt2->r;
                free(ctmpt2);
                ctmpt2 = tmp;
            }

            ctmpt1 = ctmpt1->b;
        }
    }
}
/*
+---------------------------------------------------------------------------+
| void FreeVLt() – this function frees the memory assigned to the linked    |
| list of VLt type.                                                         |
+---------------------------------------------------------------------------+
*/
void FreeVLt(VLt *vlt) {

    if (vlt->head != NULL) {

        VPt *ctvpt = vlt->head, *tmp;

        while (ctvpt != NULL) {
            tmp = ctvpt->next;
            FreeNLt(ctvpt->nlt);
            free(ctvpt);
            ctvpt = tmp;
        }
    }
}
/*
+---------------------------------------------------------------------------+
| AddPtToNLt() – this function adds data of NPt type to the linked list of  |
| NLt type.                                                                 |
+---------------------------------------------------------------------------+
*/
int AddPtToNLt(NLt *nlt, long long int no, char sgn) {

    NPt *npt = (NPt*)malloc(sizeof(NPt));

    int ret = 1;

    if (npt == NULL) ret = -1;
    else {
        npt->no = no;
        npt->sgn = sgn;

        npt->next = NULL;
        npt->prev = NULL;

        if (nlt->head == NULL) {
            nlt->head = npt;
            nlt->tail = npt;
        }
        else {
            nlt->tail->next = npt;
            npt->prev = nlt->tail;
            nlt->tail = nlt->tail->next;
        }
    }
    return ret;
}
/*
+---------------------------------------------------------------------------+
| AddPtToMLt() – this function adds data of MPt type to the linked list of  |
| MLt type.                                                                 |
+---------------------------------------------------------------------------+
*/
int AddPtToMLt(int k, MLt *mlt, int no, int cl, char sgn) {

    MPt *mpt = (MPt*)malloc(sizeof(MPt));

    int ret = 1;

    if (mpt == NULL) ret = -1;
    else {
        mpt->no = no;
        mpt->cl = cl;
        mpt->n  = sgn;
        mpt->r  = NULL;
        mpt->b  = NULL;

        if (mlt->rt == NULL) {
            mlt->rt = mpt;
            mlt->r = mpt;
            mlt->b = mpt;
        }
        else {
            if (k == 1) {
                mlt->r->r = mpt;
                mlt->r = mpt;
            }
            else
                if (k == 2) {
                    mlt->b->b = mpt;
                    mlt->b = mpt;
                    mlt->r = mpt;
                }
        }
    }

    return ret;
}
/*
+---------------------------------------------------------------------------+
| AddPtToVLt() – this function adds data of VPt type to the linked list of  |
| VLt type.                                                                 |
+---------------------------------------------------------------------------+
*/
int AddPtToVLt(VLt *vlt, int inf, int no, int cl, long long int ps, long long int pe, NLt *nlt) {

    VPt *vpt = (VPt*)malloc(sizeof(VPt));

    int ret = 1;

    if (vpt == NULL) ret = -1;
    else {
        vpt->inf = inf;

        vpt->no = no;
        vpt->cl = cl;
        vpt->ps = ps;
        vpt->pe = pe;

        vpt->nlt = nlt;

        vpt->next = NULL;
        vpt->prev = NULL;

        if (vlt->head == NULL) {
            vlt->head = vpt;
            vlt->tail = vpt;
        }
        else {
            vlt->tail->next = vpt;
            vpt->prev = vlt->tail;
            vlt->tail = vlt->tail->next;
        }
    }

    return ret;
}
/*
+---------------------------------------------------------------------------+
| PrtNLtToFile() – this function records the content of the linked list in  |
| a file.                                                                   |
+---------------------------------------------------------------------------+
*/
int PrtNLtToFile(FILE *f, NPt *npt) {

    NPt *ctnpt = npt;

    int ret = 1;

    while (ctnpt != NULL) {

        if (fprintf(f,"%c",ctnpt->sgn) < 0) {
            ret = -1;
            break;
        }

        ctnpt = ctnpt->next;
    }

    return ret;
}
/*
+---------------------------------------------------------------------------+
| ImportGenome() – this function imports the genome from the file           |
| sequence.fasta and records it in the linked list of NLt type.             |
+---------------------------------------------------------------------------+
*/
long long int ImportGenome(char *fn, NLt *nlt) {

    long long int ret = 0;

    char sgn;

    FILE *f;

    if ((f = fopen(fn,"r")) != NULL) {

        while ((sgn = fgetc(f)) != EOF) {

            if (sgn == '>') {

                while ((sgn = fgetc(f)) != EOF) if (sgn == '\n') break;

                if (ferror(f) != 0)  {
                    fclose(f);
                    ret = -110;
                }
            }
            else {

                if (sgn != '\n' && sgn != '\t' && sgn != '\v' && sgn != ' '  && sgn != '\0' && sgn != '\b' && sgn != '\r' && sgn != EOF  && sgn != 'a'  && sgn !='\f')
                {
                    ret++;

                    if (AddPtToNLt(nlt,ret,sgn) < 0) {
                        fclose(f);
                        ret = -120;
                    }
                }
            }
        }

        if (ferror(f) == 0) fclose(f);
        else {
            fclose(f);
            ret = -110;
        }
    }
    else ret = -100;

    if ((ret != -100) && (nlt->head == NULL)) ret = -130;

    return ret;
}
/*
+---------------------------------------------------------------------------+
| CopyNLt() – this function creates a copy of the linked list of NLt type.  |
+---------------------------------------------------------------------------+
*/
int CopyNLt(NLt *rnlt, NLt *nlt) {

    NPt *ctnpt = nlt->head;

    int ret = 1;

    if (ctnpt != NULL) {

        while (ctnpt != NULL) {

            if (AddPtToNLt(rnlt,ctnpt->no,ctnpt->sgn) < 0) {
                ret = -150;
                break;
            }

            ctnpt = ctnpt->next;
        }
    }
    else ret = -140;

    return ret;
}
/*
+---------------------------------------------------------------------------+
| JoinNLtWithNLt() – this function combines two linked lists of NLt type in |
| such a way that the end of one linked list is added another linked list.  |
+---------------------------------------------------------------------------+
*/
void JoinNLtWithNLt(NLt *nlt1, NLt *nlt2) {

    if (nlt2->head != NULL) {

        if (nlt1->head == NULL) {
            nlt1->head = nlt2->head;
            nlt1->tail = nlt2->tail;
        }
        else {
            nlt1->tail->next = nlt2->head;
            nlt2->head->prev = nlt1->tail;
            nlt1->tail = nlt2->tail;
        }
    }
}
/*
+---------------------------------------------------------------------------+
| ImportTRSToMLt() – this function imports TRS motifs from the file         |
| trs.txt and records them in the linked list of MLt type.                  |
+---------------------------------------------------------------------------+
*/
int ImportTRSToMLt(char *fn, MLt *mlt) {

    int p = 0, pl = 0, no = 0, cl = 1, reset = 0;

    char sgn, l;

    FILE *f;

    if ((f = fopen(fn,"r")) != NULL) {

        if ((l = fgetc(f)) != EOF) {

            if (l != '#') pl = -170;
        }
        else pl = -170;

        if (pl >= 0) {

            while ((sgn = fgetc(f)) != EOF) {

                if (sgn == 'A' ||  sgn == 'C' || sgn == 'G' || sgn == 'T' || sgn == '#' || sgn == '\n') {

                    if (reset == 1) {
                        reset = 0;
                        p = 0;
                    }

                    if ((sgn == 'A' || sgn == 'C' || sgn == 'G' || sgn == 'T') && (l == '#')) {

                        p++;

                        no++;

                        if (AddPtToMLt(2,mlt,no,cl,sgn) < 0) {
                            pl = -180;
                            break;
                        }
                    }

                    if ((sgn == 'A' || sgn == 'C' || sgn == 'G' || sgn == 'T') && (l == 'A' || l == 'C' || l == 'G' || l == 'T')) {

                        p++;

                        if (AddPtToMLt(1,mlt,no,cl,sgn) < 0) {
                            pl = -180;
                            break;
                        }
                    }

                    if ((sgn == '#') && (l == 'A' || l == 'C' || l == 'G' || l == 'T')) {

                        if ((p % 3) != 0) {
                            pl = -190;
                            break;
                        }
                        else {

                            if (pl == 0) {
                                reset = 1;
                                pl = p;
                            }
                            else {

                                if (pl == p) reset = 1;
                                else {
                                    pl = -190;
                                    break;
                                }
                            }
                        }
                    }

                    if ((sgn == '\n') && (l == 'A' || l == 'C' || l == 'G' || l == 'T')) {

                        cl++;

                        if ((p % 3) != 0) {
                            pl = -190;
                            break;
                        }
                        else {

                            if (pl == 0) {
                                reset = 1;
                                pl = p;
                            }
                            else {

                                if (pl == p) reset = 1;
                                else {
                                    pl = -190;
                                    break;
                                }
                            }
                        }
                    }

                    if ((sgn == '\n' && l == '\n') || (sgn == '#' && l == '#') || (sgn == '\n' && l == '#')) {
                        pl = -170;
                        break;
                    }

                    l = sgn;
                }
                else {
                    pl = -170;
                    break;
                }
            }
        }

        if (pl >= 0) {

            if (ferror(f) == 0) {

                if ((p % 3) != 0) {
                    pl = -190;
                    fclose(f);
                }
                else
                {
                    if (pl == p) fclose(f);
                    else {
                        pl = -190;
                        fclose(f);
                    }
                }
            }
            else {
                pl = -200;
                fclose(f);
            }
        }
    }
    else pl = -160;

    return pl;
}
/*
+---------------------------------------------------------------------------+
|  FirstToLastInListNLt() – this function shifts the first element of the   |
|  linked list NLt to the end of the linked list and varies the position and|
|  sign representing the nucleotide.                                        |
+---------------------------------------------------------------------------+
*/
void FirstToLastInListNLt(NLt *nlt, long long int no, char sgn) {

    NPt *first = nlt->head->next;

    nlt->head->next->prev = NULL;
    nlt->head->no = no;
    nlt->head->sgn = sgn;
    nlt->head->next = NULL;
    nlt->head->prev = nlt->tail;
    nlt->tail->next = nlt->head;
    nlt->tail = nlt->head;
    nlt->head = first;
}
/*
+---------------------------------------------------------------------------+
| UpLe() – this function replace a small letter representing nucleotide     |
| into a capital letter representing nucleotide.                            |
+---------------------------------------------------------------------------+
*/
char UpLe(char sgn)
{
    if(sgn>=97) sgn = sgn-32;

    return sgn;
 }
/*
+---------------------------------------------------------------------------+
|  CompareOneTRS() – this function checks if there occurs one TRS motif.    |
+---------------------------------------------------------------------------+
*/
int CompareOneTRS(MPt *mpt, NPt *npt) {

    MPt *ctmpt = mpt;
    NPt *ctnpt = npt;

    int ret = 1;

    do {

        if (UpLe(ctmpt->n) != UpLe(ctnpt->sgn)) {
            ret = 0;
            break;
        }
        else {
            ctmpt = ctmpt->r;
            ctnpt = ctnpt->next;
        }

    } while (ctnpt != NULL);

    return ret;
}
/*
+---------------------------------------------------------------------------+
| FindFistrTRSPosition() – this function searches the left TRS flanking     |
| sequence.                                                                 |
+---------------------------------------------------------------------------+
*/
VPt* FindFistrTRSPosition(MLt *mlt, NLt *nlt) {

    NPt *npt3 = nlt->head, *npt1 = nlt->head;

    MPt *ctmpt = mlt->rt;
    NPt *ctnpt = nlt->head;

    int ret = 0;

    VPt *vpt = (VPt*)malloc(sizeof(VPt)), *rvpt = NULL;

    if (vpt != NULL) {

        vpt->inf = 0;

        vpt->no = 0;
        vpt->cl = 0;
        vpt->ps = 0;
        vpt->pe = 0;

        vpt->nlt = NULL;

        vpt->next = NULL;
        vpt->prev = NULL;

        while (ctmpt != NULL) {

            if (CompareOneTRS(ctmpt, ctnpt) == 1) {

                NLt  *nlt1 = (NLt*)malloc(sizeof(NLt));

                if (nlt1 != NULL) {

                    InitNLt(nlt1);

                    for(int i = 1; i <= 3; i++) {

                        if (AddPtToNLt(nlt1,npt1->no,npt1->sgn) < 0) {
                            ret = -1;
                            break;
                        }

                        npt1 = npt1->next;
                    }

                    if (ret == 0) {

                        vpt->inf = 1;

                        vpt->no = ctmpt->no;
                        vpt->cl = ctmpt->cl;
                        vpt->ps = ctnpt->no;
                        vpt->pe = -1;

                        vpt->nlt = nlt1;
                    }
                    else {
                        FreeNLt(nlt1);
                        break;
                    }
                }
                else {
                    ret = -1;
                    break;
                }

                break;
            }

            ctmpt = ctmpt->b;
        }

    }
    else ret = -1;

    if (ret == -1) {
        rvpt = vpt;
        vpt = NULL;
        free(rvpt);
    }

    return vpt;

}
/*
+---------------------------------------------------------------------------+
|  CompareLists() – this function checks if the TRS repeats itself.         |
+---------------------------------------------------------------------------+
*/
int CompareLists(NPt *npt1, NPt *npt2) {

    NPt *ctnpt1 = npt1;
    NPt *ctnpt2 = npt2;

    int ret = 1;

    do {

        if (ctnpt2 == NULL) {
            ret = -1;
            break;
        }

        if (UpLe(ctnpt1->sgn) != UpLe(ctnpt2->sgn)) {
            ret = 0;
            break;
        }
        else {
            ctnpt1 = ctnpt1->next;
            ctnpt2 = ctnpt2->next;
        }

    } while (ctnpt1 != NULL);

    return ret;
}
/*
+---------------------------------------------------------------------------+
|  LC_TRSPositionsFindAndSaveToVLt() – this function finds the positions   |
|  of flanking sequences and records them in the linked list VLt together   |
|  with their TRSs which make the base of these flanking sequences for the  |
|  linear case.                                                             |
+---------------------------------------------------------------------------+
*/
int LC_TRSPositionsFindAndSaveToVLt(MLt *mlt, NLt *nlt, VLt *vlt, int len) {

    int i, r, ret = 1, con = 0;

    NPt *ctnpt = nlt->head;

    VPt *vpt;

    NLt *snlt = (NLt*)malloc(sizeof(NLt));

    if (snlt != NULL) {

        InitNLt(snlt);

        for(i = 1; i <= len; i++) {

            if (ctnpt != NULL) {

                if ((AddPtToNLt(snlt,ctnpt->no,ctnpt->sgn)) < 0) {
                    ctnpt = NULL;
                    ret = -220;
                    break;
                }
            }
            else {
                ret = -230;
                break;
            }

            if (i != len) ctnpt = ctnpt->next;
        }

        if ((i - 1) == len) con = 3;
    }
    else {
        ctnpt = NULL;
        ret = -210;
    }

    while (ctnpt != NULL) {

        if (con == 1) {

            for(i = 1; i <= len; i++) {

                if (ctnpt != NULL) FirstToLastInListNLt(snlt,ctnpt->no,ctnpt->sgn);
                else {
                    con = 9;
                    break;
                }

                if (i != len) ctnpt = ctnpt->next;
            }

            if ((i - 1) == len) con = 3;
        }

        if (con == 2) {

            FirstToLastInListNLt(snlt,ctnpt->no,ctnpt->sgn);

            con = 3;
        }

        if (con == 3) {

            if ((vpt = FindFistrTRSPosition(mlt,snlt)) == NULL) {
                ret = -240;
                break;
            }
            else {

                if (vpt->inf == 0) {

                    VPt *rvpt = vpt;
                    vpt = NULL;
                    free(rvpt);

                    ctnpt = ctnpt->next;
                    con = 2;
                }
                else con = 4;
            }
        }

        if (con == 4) {

            if ((r = CompareLists(vpt->nlt->head,ctnpt->next)) == 0) {

                if (AddPtToVLt(vlt,vpt->inf,vpt->no,vpt->cl,vpt->ps,ctnpt->no,vpt->nlt) < 0) {
                    ret = -240;
                    break;
                }
                else {
                    VPt *rvpt = vpt;
                    vpt = NULL;
                    free(rvpt);
                }

                ctnpt = ctnpt->next;

                con = 1;
            }
            else {

                if (r == 1) {

                    if ((ctnpt->next != NULL) && (ctnpt->next->next != NULL) && (ctnpt->next->next->next != NULL)) ctnpt = ctnpt->next->next->next;

                    con = 4;
                }
                else {

                    if (AddPtToVLt(vlt,vpt->inf,vpt->no,vpt->cl,vpt->ps,ctnpt->no,vpt->nlt) < 0) {
                        ret = -240;
                        break;
                    }
                    else {
                        VPt *rvpt = vpt;
                        vpt = NULL;
                        free(rvpt);
                    }

                    con = 9;
                }
            }
        }

        if (con == 9) {
            ret = 1;
            break;
        }
    }

    return ret;
}
/*
+---------------------------------------------------------------------------+
|  CC_TRSPositionsFindAndSaveToVLt() – this function finds the positions of|
|  flanking sequences and records them in the linked list VLt together with |
|  their TRSs which make the base for these flanking sequences for the      |
|  linear circular case.                                                    |
+---------------------------------------------------------------------------+
*/
int CC_TRSPositionsFindAndSaveToVLt(MLt *mlt, NLt *nlt, VLt *vlt, int len) {

    int i, r, k = 0, ret = 1, con = 0, stop = 0;

    NPt *ctnpt = nlt->head;
    VPt *vpt;

    NLt *snlt = (NLt*)malloc(sizeof(NLt));

    if (snlt != NULL) {

        InitNLt(snlt);

        for(i = 1; i <= len; i++) {

            if (ctnpt != NULL) {

                if (AddPtToNLt(snlt,ctnpt->no,ctnpt->sgn) < 0) {
                    ctnpt = NULL;
                    ret = -260;
                    break;
                }

                if (k < ctnpt->no) k = ctnpt->no;
                else stop = 1;
            }
            else {
                ret = -270;
                break;
            }

            if (i != len) ctnpt = ctnpt->next;
        }

        if ((i - 1) == len) con = 3;
    }
    else {
        ctnpt = NULL;
        ret = -250;
    }

    while (ctnpt != NULL) {

        if (stop == 2) con = 9;

        if (con == 1) {

            for(i = 1; i <= len; i++) {

                if (ctnpt != NULL) {

                    FirstToLastInListNLt(snlt,ctnpt->no,ctnpt->sgn);

                    if (k < ctnpt->no) k = ctnpt->no;
                    else stop = 1;
                }
                else {
                    con = 9;
                    break;
                }

                if (i != len) ctnpt = ctnpt->next;
            }

            if ((i - 1) == len) con = 3;
        }

        if (con == 2) {

            FirstToLastInListNLt(snlt,ctnpt->no,ctnpt->sgn);

            if (k < ctnpt->no) k = ctnpt->no;
            else stop = 1;

            con = 3;
        }

        if (con == 3) {

            if ((vpt = FindFistrTRSPosition(mlt,snlt)) == NULL) {
                ret = -280;
                break;
            }
            else {

                if (vpt->inf == 0) {

                    VPt *rvpt = vpt;
                    vpt = NULL;
                    free(rvpt);

                    ctnpt = ctnpt->next;
                    con = 2;
                }
                else con = 4;
            }
        }

        if (con == 4)
        {
            if ((r = CompareLists(vpt->nlt->head,ctnpt->next)) == 0) {

                if (AddPtToVLt(vlt,vpt->inf,vpt->no,vpt->cl,vpt->ps,ctnpt->no,vpt->nlt) < 0) {
                    ret = -280;
                    break;
                }

                if (stop != 1) con = 1;
                else {
                    if (vpt->ps < ctnpt->no) {
                        con = 1;
                        stop = 2;
                    }
                    else con = 1;
                }

                VPt *rvpt = vpt;
                vpt = NULL;
                free(rvpt);

                ctnpt = ctnpt->next;
            }
            else {

                if (r == 1) {

                    if ((ctnpt->next != NULL) && (ctnpt->next->next != NULL) && (ctnpt->next->next->next != NULL)) {

                        ctnpt = ctnpt->next->next->next;

                        if (k < ctnpt->no) k = ctnpt->no;
                        else stop = 1;
                    }

                    con = 4;
                }
                else {

                    if (AddPtToVLt(vlt,vpt->inf,vpt->no,vpt->cl,vpt->ps,ctnpt->no,vpt->nlt) < 0) {
                        ret = -280;
                        break;
                    }
                    else {
                        VPt *rvpt = vpt;
                        vpt = NULL;
                        free(rvpt);
                    }

                    con = 9;
                }
            }
        }

        if (con == 9) {
            ret = 1;
            break;
        }
    }

    return ret;
}
/*
+---------------------------------------------------------------------------+
|  LC_InteriorsFindAndSaveToFile() – this function generates the file       |
|  interiors.txt for the linear case with conditions.                       |
+---------------------------------------------------------------------------+
*/
int LC_InteriorsFindAndSaveToFile(char *ifn, VLt *vlt, NLt *nlt, long long int min, long long int max) {

    NPt *ctnpt = nlt->head;
    VPt *ctvpt = vlt->head;

    int ret = 1;

    FILE *f;

    if ((f = fopen(ifn,"w")) != NULL) {

        if (fprintf(f,"L-NoClass;L-No;LFS;Len(LFS);L-POS(LFS);R-POS(LFS);R-NoClass;R-No;RFS;Len(RFS);L-POS(RFS);R-POS(RFS);>SEQ;Len(SEQ)\n") < 0) {
            ret = -300;
            ctvpt->next = NULL;
        }

        while (ctvpt->next != NULL) {

            if (((ctvpt->next->ps - ctvpt->pe) > min) && ((ctvpt->next->ps - ctvpt->pe) < (max +2))) {

                if ((ctvpt->cl > 0)&&(ctvpt->cl < 10)) {

                    if ((ctvpt->no > 0)&&(ctvpt->no < 10)) {

                        if (fprintf(f,"0%d;0%d;",ctvpt->cl,ctvpt->no) < 0) {
                            ret = -300;
                            break;
                        }
                    }
                    else {

                        if (fprintf(f,"0%d;%d;",ctvpt->cl,ctvpt->no) < 0) {
                            ret = -300;
                            break;
                        }
                    }
                }
                else {

                    if ((ctvpt->no > 0)&&(ctvpt->no < 10)) {

                        if (fprintf(f,"%d;0%d;",ctvpt->cl,ctvpt->no) < 0) {
                            ret = -300;
                            break;
                        }
                    }
                    else {

                        if (fprintf(f,"%d;%d;",ctvpt->cl,ctvpt->no) < 0) {
                            ret = -300;
                            break;
                        }
                    }
                }

                for(int i = 1; i <= (abs(ctvpt->pe - ctvpt->ps + 1)/3); i++) {

                    if (PrtNLtToFile(f,ctvpt->nlt->head) < 0) {
                        ret = -300;
                        break;
                    }
                }

                if (ret == -300) break;

                if ((ctvpt->next->cl > 0)&&(ctvpt->next->cl < 10))
                {
                    if ((ctvpt->next->no > 0)&&(ctvpt->next->no < 10)) {

                        if (fprintf(f,";%d;%d;%d;0%d;0%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                            ret = -300;
                            break;
                        }
                    }
                    else {
                        if (fprintf(f,";%d;%d;%d;0%d;%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                            ret = -300;
                            break;
                        }
                    }
                }
                else {

                    if ((ctvpt->next->no > 0)&&(ctvpt->next->no < 10)) {

                        if (fprintf(f,";%d;%d;%d;%d;0%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps, ctvpt->pe, ctvpt->next->cl,ctvpt->next->no) < 0) {
                            ret = -300;
                            break;
                        }
                    }
                    else {
                        if (fprintf(f,";%d;%d;%d;%d;%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                            ret = -300;
                            break;
                        }
                    }
                }

                for(int i = 1; i <= (abs(ctvpt->next->pe - ctvpt->next->ps + 1)/3); i++) {

                    if (PrtNLtToFile(f,ctvpt->next->nlt->head) < 0) {
                        ret = -300;
                        break;
                    }
                }

                if (ret == -300) break;

                if (fprintf(f,";%d;%d;%d;>",abs(ctvpt->next->pe - ctvpt->next->ps + 1),ctvpt->next->ps,ctvpt->next->pe) < 0) {
                    ret = -300;
                    break;
                }

                while (ctvpt->pe != ctnpt->no) {
                    ctnpt = ctnpt->next;
                }

                ctnpt = ctnpt->next;

                int li = 0;

                while (ctnpt->no != ctvpt->next->ps) {

                    if (fprintf(f,"%c",ctnpt->sgn) < 0) {
                        ret = -300;
                        break;
                    }

                    li++;

                    ctnpt = ctnpt->next;
                }

                if (ret == -300) break;

                if (fprintf(f,";%d\n",li) < 0) {
                    ret = -300;
                    break;
                }
            }

            ctvpt = ctvpt->next;
        }

        fclose(f);
    }
    else ret = -290;

    return ret;
}
/*
+---------------------------------------------------------------------------+
|  CC_InteriorsFindAndSaveToFile() – this function generates the file       |
|  interiors.txt for the circular case with conditions.                     |
+---------------------------------------------------------------------------+
*/
int CC_InteriorsFindAndSaveToFile(char *ifn, VLt *vlt, NLt *nlt, long long int size, long long int min, long long int max) {

    NPt *ctnpt = nlt->head;
    VPt *ctvpt = vlt->head;

    int ret = 1;

    FILE *f;

    if ((f = fopen(ifn,"w")) != NULL) {

        if (fprintf(f,"L-NoClass;L-No;LFS;Len(LFS);L-POS(LFS);R-POS(LFS);R-NoClass;R-No;RFS;Len(RFS);L-POS(RFS);R-POS(RFS);>SEQ;Len(SEQ)\n") < 0) {
            ret = -320;
            ctvpt->next = NULL;
        }

        while (ctvpt->next != NULL) {

            if (ctvpt->next->ps > ctvpt->pe) {

                if (((ctvpt->next->ps - ctvpt->pe) > min) && ((ctvpt->next->ps - ctvpt->pe) < (max +2))) {

                    if ((ctvpt->cl > 0)&&(ctvpt->cl < 10)) {

                        if ((ctvpt->no > 0)&&(ctvpt->no < 10)) {

                            if (fprintf(f,"0%d;0%d;",ctvpt->cl,ctvpt->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                        else {
                            if (fprintf(f,"0%d;%d;",ctvpt->cl,ctvpt->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                    }
                    else {

                        if ((ctvpt->no > 0)&&(ctvpt->no < 10)) {
                            if (fprintf(f,"%d;0%d;",ctvpt->cl,ctvpt->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                        else {
                            if (fprintf(f,"%d;%d;",ctvpt->cl,ctvpt->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                    }

                    for(int i = 1; i <= (abs(ctvpt->pe - ctvpt->ps + 1)/3); i++) {

                        if (PrtNLtToFile(f,ctvpt->nlt->head) < 0) {
                            ret = -320;
                            break;
                        }
                    }

                    if (ret == -320) break;

                    if ((ctvpt->next->cl > 0)&&(ctvpt->next->cl < 10)) {

                        if ((ctvpt->next->no > 0)&&(ctvpt->next->no < 10)) {

                            if (fprintf(f,";%d;%d;%d;0%d;0%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                        else {

                            if (fprintf(f,";%d;%d;%d;0%d;%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                    }
                    else {

                        if ((ctvpt->next->no > 0)&&(ctvpt->next->no < 10)) {

                            if (fprintf(f,";%d;%d;%d;%d;0%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                        else {

                            if (fprintf(f,";%d;%d;%d;%d;%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                    }

                    for(int i = 1; i <= (abs(ctvpt->next->pe - ctvpt->next->ps + 1)/3); i++) {

                        if (PrtNLtToFile(f,ctvpt->next->nlt->head) < 0) {
                            ret = -320;
                            break;
                        }
                    }

                    if (ret == -320) break;

                    if (fprintf(f,";%d;%d;%d;>",abs(ctvpt->next->pe - ctvpt->next->ps + 1),ctvpt->next->ps,ctvpt->next->pe) < 0) {
                        ret = -320;
                        break;
                    }

                    while (ctvpt->pe != ctnpt->no){
                        ctnpt = ctnpt->next;
                    }

                    ctnpt = ctnpt->next;

                    int li = 0;

                    while (ctnpt->no != ctvpt->next->ps) {

                        if (fprintf(f,"%c",ctnpt->sgn) < 0) {
                            ret = -320;
                            break;
                        }

                        li++;

                        ctnpt = ctnpt->next;
                    }

                    if (ret == -320) break;

                    if (fprintf(f,";%d\n",li) < 0) {
                        ret = -320;
                        break;
                    }
                }

                ctvpt = ctvpt->next;
            }
            else {

                if ((((size - ctvpt->pe) + ctvpt->next->ps) > min) && (((size - ctvpt->pe) + ctvpt->next->ps) < (max +2))) {

                    if ((ctvpt->cl > 0)&&(ctvpt->cl < 10)) {

                        if ((ctvpt->no > 0)&&(ctvpt->no < 10)) {

                            if (fprintf(f,"0%d;0%d;",ctvpt->cl,ctvpt->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                        else {

                            if (fprintf(f,"0%d;%d;",ctvpt->cl,ctvpt->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                    }
                    else {

                        if ((ctvpt->no > 0)&&(ctvpt->no < 10)) {

                            if (fprintf(f,"%d;0%d;",ctvpt->cl,ctvpt->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                        else {

                            if (fprintf(f,"%d;%d;",ctvpt->cl,ctvpt->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                    }

                    for(int i = 1; i <= (abs(ctvpt->pe - ctvpt->ps + 1)/3); i++) {

                        if (PrtNLtToFile(f,ctvpt->nlt->head) < 0) {
                            ret = -320;
                            break;
                        }
                    }

                    if (ret == -320) break;

                    if ((ctvpt->next->cl > 0)&&(ctvpt->next->cl < 10)) {

                        if ((ctvpt->next->no > 0)&&(ctvpt->next->no < 10)) {

                            if (fprintf(f,";%d;%d;%d;0%d;0%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                        else {
                            if (fprintf(f,";%d;%d;%d;0%d;%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                    }
                    else {

                        if ((ctvpt->next->no > 0)&&(ctvpt->next->no < 10)) {

                            if (fprintf(f,";%d;%d;%d;%d;0%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                        else {

                            if (fprintf(f,";%d;%d;%d;%d;%d;",abs(ctvpt->pe - ctvpt->ps + 1),ctvpt->ps,ctvpt->pe,ctvpt->next->cl,ctvpt->next->no) < 0) {
                                ret = -320;
                                break;
                            }
                        }
                    }

                    for(int i = 1; i <= (abs(ctvpt->next->pe - ctvpt->next->ps + 1)/3); i++) {

                        if (PrtNLtToFile(f,ctvpt->next->nlt->head) < 0) {
                            ret = -320;
                            break;
                        }
                    }

                    if (ret == -320) break;

                    if (fprintf(f,";%d;%d;%d;>",abs(ctvpt->next->pe - ctvpt->next->ps + 1),ctvpt->next->ps,ctvpt->next->pe) < 0) {
                        ret = -320;
                        break;
                    }

                    while (ctvpt->pe != ctnpt->no) {
                        ctnpt = ctnpt->next;
                    }

                    ctnpt = ctnpt->next;

                    int li = 0;

                    while (ctnpt->no != ctvpt->next->ps) {

                        if (fprintf(f,"%c",ctnpt->sgn) < 0) {
                            ret = -320;
                            break;
                        }

                        li++;

                        ctnpt = ctnpt->next;
                    }

                    if (ret == -320) break;

                    if (fprintf(f,";%d\n",li) < 0) {
                        ret = -320;
                        break;
                    }
                }

                ctvpt = ctvpt->next;
            }
        }

        fclose(f);
    }
    else ret = -310;

    return ret;
}
/*
+---------------------------------------------------------------------------+
|  PrintError() – this function returns information on the basic errors     |
|  served.                                                                  |
+---------------------------------------------------------------------------+
*/
void PrintError(int e) {

    printf("\nError: ");

    if (e == -100) printf("%d The file sequence.fasta is missing, or there is an error with opening sequence.fasta.\n\n", e);
    if (e == -110) printf("%d The data (character) reading error from file sequence.fasta.\n\n", e);
    if (e == -120) printf("%d The data saving error to the linked list NLt.\n\n", e);
    if (e == -130) printf("%d The file sequence.fasta is nucleotides missing.\n\n", e);
    if (e == -140) printf("%d Initialization error of an additional linked list NLt.\n\n", e);
    if (e == -150) printf("%d The computer memory (RAM) is full or an error while writing a genome nucleotide into memory.\n\n", e);
    if (e == -160) printf("%d The file trs.txt is missing, or there is an error with opening trs.txt.\n\n", e);
    if (e == -170) printf("%d The file trs.txt is not formatted correctly. Check record in the trs.txt file.\n\n", e);
    if (e == -180) printf("%d The computer memory (RAM) is full or incorrect record in the trs.txt file.\n\n", e);
    if (e == -190) printf("%d Check the record of TRS motif in the in the trs.txt file \n\n", e);
    if (e == -200) printf("%d Incorrect read data (character) from the trs.txt file.\n\n", e);
    if (e == -210) printf("%d The computer memory (RAM) is full or incorrect record of the TRS motif position in the RAM.\n\n", e);
    if (e == -220) printf("%d The data saving error to the linked list NLt during searching for a TRS motif.\n\n", e);
    if (e == -230) printf("%d The genome taken from the sequence.fasta file contains less nucleotides than the TRS motif.\n\n", e);
    if (e == -240) printf("%d The computer memory (RAM) is full or an error during writing the flanking sequence of the linear case.\n\n", e);
    if (e == -250) printf("%d The computer memory (RAM) is full or an error during writing of the TRS motif position in the RAM of the circular case.\n\n", e);
    if (e == -260) printf("%d The data saving error to the linked list NLt or during searching for a TRS motif of the circular case.\n\n", e);
    if (e == -270) printf("%d The genome taken from the sequence.fasta file contains less nucleotides than the TRS motif of the circular case.\n\n", e);
    if (e == -280) printf("%d The computer memory (RAM) is full or an error during writing the flanking sequence of the circular case.\n\n", e);
    if (e == -290) printf("%d The file interiors.txt creating error of the linear case.\n\n", e);
    if (e == -300) printf("%d The file interiors.txt saving error of the linear case.\n\n", e);
    if (e == -310) printf("%d The file interiors.txt creating error of the circular case.\n\n", e);
    if (e == -320) printf("%d The file interiors.txt saving error of the circular case.\n\n", e);

}
/*
+---------------------------------------------------------------------------+
| ExitN() – this function free the memory assigned to the linked lists      |
| of NLt types after the occurrence of one of the basic errors served.      |                                                              |
+---------------------------------------------------------------------------+
*/
void ExitN(NLt *nlt, int e) {

    FreeNLt(nlt);
    PrintError(e);
    int t; scanf("%i",&t);
    exit(1);
}
/*
+---------------------------------------------------------------------------+
| ExitNN() – this function free the memory assigned to the linked lists     |
| of NLt types after the occurrence of one of the basic errors served.      |                                                              |
+---------------------------------------------------------------------------+
*/
void ExitNN(NLt *nlt1, NLt *nlt2, int e) {

    FreeNLt(nlt1);
    FreeNLt(nlt2);
    PrintError(e);
    int t; scanf("%i",&t);
    exit(1);
}
/*
+---------------------------------------------------------------------------+
| ExitNM() – this function free the memory assigned to the linked lists     |
| of MLt types after the occurrence of one of the basic errors served.      |                                                              |
+---------------------------------------------------------------------------+
*/
void ExitNM(NLt *nlt, MLt *mlt, int e) {

    FreeNLt(nlt);
    FreeMLt(mlt);
    PrintError(e);
    int t; scanf("%i",&t);
    exit(1);
}
/*
+---------------------------------------------------------------------------+
| ExitNMV() – this function free the memory assigned to the linked lists     |
| of VLt types after the occurrence of one of the basic errors served.      |                                                              |
+---------------------------------------------------------------------------+
*/
void ExitNMV(NLt *nlt, MLt *mlt, VLt *vlt, int e) {

    FreeNLt(nlt);
    FreeMLt(mlt);
    FreeVLt(vlt);
    PrintError(e);
    int t; scanf("%i",&t);
    exit(1);
}
int main()
{
    //char *gfn = "sequence.fasta", *tfn = "trs.txt", *ifn = "interiors.txt";
    char *tfn = "trs.txt", *ifn = "interiors.txt";

    int  e = 0, opt = 0, len = 0;

    long long int tsize = 0;

    float tmin, tmax;
    char gfn[256];

    NLt gl; InitNLt(&gl);
    NLt gl2; InitNLt(&gl2);
    MLt ml; InitMLt(&ml);
    VLt vl; InitVLt(&vl);

    printf("\n");
    printf("1. Analysis of the linear case with conditions \n\n");
    printf("2. Analysis of the circular case with conditions \n\n");
    printf("Select option: ");

    // debug part
    // tmin = 3.0;
    // tmax = 6.0;
    // gfn = "sequence.fasta";

    // if ((tsize = ImportGenome(gfn,&gl)) < 0) ExitN(&gl,(int)tsize);
    // if ((len = ImportTRSToMLt(tfn,&ml)) < 0) ExitNM(&gl,&ml,len);
    // if ((e = LC_TRSPositionsFindAndSaveToVLt(&ml,&gl,&vl,len)) < 0) ExitNMV(&gl,&ml,&vl,e);
    // if ((e = LC_InteriorsFindAndSaveToFile(ifn,&vl,&gl,tmin,tmax)) < 0) ExitNMV(&gl,&ml,&vl,e);
    // FreeNLt(&gl);
    // FreeMLt(&ml);
    // FreeVLt(&vl);
    // return 0;


    scanf("%i",&opt);

        switch (opt)
        {
            case 0: break;


            case 1: printf("\n\n Set the minimum (Min) and the maximum (Max) lenght of interior:\n");
                    printf("\n  Min="); scanf("%f",&tmin);
                    printf("\n  Max="); scanf("%f",&tmax);
                    printf("\n  Genom file="); scanf("%s", &gfn);
                    //printf("\n The fasta file: %s", &gfn);

                    if (tmin < 1) {
                        printf("\n Error: Min < 1\n\n");
                        exit (1);
                    }

                    if (tmin > tmax) {
                        printf("\n Error: Min > Max\n\n");
                        exit (1);
                    }

                    printf("\n\nSTART\n");
                    
                    if ((tsize = ImportGenome(gfn,&gl)) < 0) {
                        printf("\n issue in importing the genome...");
                        ExitN(&gl,(int)tsize);
                    }

                    if ((len = ImportTRSToMLt(tfn,&ml)) < 0) {
                        printf("\n issue in importing the TRS file...");
                        ExitNM(&gl,&ml,len);
                    }
                    if ((e = LC_TRSPositionsFindAndSaveToVLt(&ml,&gl,&vl,len)) < 0) {
                        printf("\n issue in finding the sequence...");
                        ExitNMV(&gl,&ml,&vl,e);
                    }

                    if ((e = LC_InteriorsFindAndSaveToFile(ifn,&vl,&gl,tmin,tmax)) < 0) {
                        printf("\n issue in importing the saving the resulting file...");
                        ExitNMV(&gl,&ml,&vl,e);
                    }

                    printf("END\n");

                    break;


            case 2: printf("\n\n Set the minimum (Min) and the maximum (Max) lenght of interior:\n");
                    printf("\n  Min="); scanf("%f",&tmin);
                    printf("\n  Max="); scanf("%f",&tmax);
                    printf("\n  Genom file="); scanf("%s",&gfn);

                    if (tmin < 1) {
                        printf("\n Error: Min < 1\n\n");
                        exit (1);
                    }

                    if (tmin > tmax) {
                        printf("\n Error: Min > Max\n\n");
                        exit (1);
                    }

                    printf("\n\nSTART\n");

                    if ((tsize = ImportGenome(gfn,&gl)) < 0) ExitN(&gl,(int)tsize);

                    if ((e = CopyNLt(&gl2, &gl)) < 0) ExitNN(&gl,&gl2,e); JoinNLtWithNLt(&gl,&gl2);

                    if ((len = ImportTRSToMLt(tfn,&ml)) < 0) ExitNM(&gl,&ml,len);

                    if ((e = CC_TRSPositionsFindAndSaveToVLt(&ml,&gl,&vl,len)) < 0) ExitNMV(&gl,&ml,&vl,e);

                    if ((e = CC_InteriorsFindAndSaveToFile(ifn,&vl,&gl,tsize,tmin,tmax)) < 0) ExitNMV(&gl,&ml,&vl,e);

                    printf("END\n\n");

                    break;
        }

        FreeNLt(&gl);
        FreeMLt(&ml);
        FreeVLt(&vl);

        return 0;
}
