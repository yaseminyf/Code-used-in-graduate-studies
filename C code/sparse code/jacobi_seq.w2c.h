/* Include builtin types and operators */
#include <whirl2c.h>

/* File-level types */
struct split_ld {
  _IEEE64 msd;
  _IEEE64 lsd;
};

union predef_ld_union {
  _QUAD ldpart;
  struct split_ld split;
};

union predef_ldv_union {
  _QUAD ldpart0;
  struct split_ld split0;
};

struct __file_s {
  _INT64 _cnt;
  _UINT8 * _ptr;
  _UINT8 * _base;
  _UINT8 _reserved[2LL];
  _UINT16 _flag;
  _UINT32 _file;
};

struct VECTOR {
  _INT32 dim;
  _UINT8 fill[4LL];
  _IEEE64 * ve;
};

struct MATRIX {
  _UINT32 rowdefined;
  _INT32 m;
  _INT32 n;
  _INT32 nonzeros;
  _INT32 * Iind;
  _INT32 * Jind;
  _INT32 * rowptr;
  _INT32 * colind;
  _IEEE64 * Valuecol;
  _IEEE64 * Valuerow;
};

struct anonymous53 {
  _UINT32 _class[257LL];
  _INT16 _lower[257LL];
  _INT16 _upper[257LL];
  _UINT8 _cswidth[7LL];
  _UINT8 _fill[1LL];
};

struct _csinfo {
  _UINT8 _eucwidth[3LL];
  _UINT8 _scrwidth[3LL];
  _UINT8 _mb_cur_max;
};

struct _euc_func {
  _INT32 _is_euc;
  _UINT8 fill0[4LL];
  _INT32(*_isset2)(_INT32);
  _INT32(*_isset3)(_INT32);
  _INT32(*_iscodeset)(_INT32, _INT32);
};

struct _collate_res {
  _UINT8 _coll_as_cmp;
};

struct __attr {
  struct anonymous53 * _ctype_tbl;
  struct _csinfo _csinfo0;
  _UINT8 fill2[1LL];
  struct _euc_func _euc_func0;
  struct _collate_res _collate_res0;
  _UINT8 fill1[7LL];
  void * _fill0[15LL];
};

struct comtype4 {
  _INT32 MSGLVB;
  _INT32 IERRB;
  _INT32 MAXSB;
  _INT32 MCOLS;
  _INT32 MSEQNS;
  _INT32 MDEQNS;
  _INT32 MSCONS;
  _INT32 MDCONS;
};

struct comtype5 {
  _INT32 STAGE;
  _INT32 IOUNIT;
  _INT32 MXUSED;
  _INT32 MXREQD;
  _INT32 NCOLS;
  _INT32 NSEQNS;
  _INT32 NDEQNS;
  _INT32 NSCONS;
  _INT32 NDCONS;
  _INT32 NZEQNS;
  _INT32 NZCONS;
  _INT32 NZMAX;
  _INT32 NEDGES;
  _INT32 METHOD;
  _INT32 NOFNZ;
  _INT32 NOFSUB;
  _INT32 ICPAD[28LL];
  _INT32 RKDEF1;
  _INT32 RKDEF2;
  _INT32 NSHORT;
  _INT32 NLONG;
  _INT32 NSP1;
  _INT32 NLP1;
};

struct comtype6 {
  _INT32 PERM;
  _INT32 INVP;
  _INT32 RHS;
  _INT32 DIAG;
  _INT32 XRNZ;
  _INT32 RNZ;
  _INT32 XNZSUB;
  _INT32 NZSUB;
  _INT32 DSEQNS;
  _INT32 DSBEQN;
  _INT32 DSCONS;
  _INT32 DSBCON;
  _INT32 ROWMSK;
  _INT32 XZROWS;
  _INT32 IMPAD[35LL];
  _INT32 LLFREE;
};

struct comtype1 {
  _INT32 IPRNTE;
  _INT32 IPRNTS;
  _INT32 MAXINT;
  _IEEE32 RATIOS;
  _IEEE32 RATIOL;
  _IEEE32 MCHEPS;
  _IEEE32 TIME;
};

struct comtype2 {
  _INT32 mdeg;
  _INT32 ehead;
  _INT32 tag;
  _INT32 mdnode;
};

