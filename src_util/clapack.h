#ifndef __CLAPACK_H
#define __CLAPACK_H

/** \file clapack.h
\brief Contains generic headers for CLAPACK.
\ingroup generic_headers
*/

#ifdef __cplusplus
extern "C" {
#endif
/* Subroutine */ int cbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, f2c_real *d__, f2c_real *e, f2c_complex *vt, integer *ldvt, 
	f2c_complex *u, integer *ldu, f2c_complex *c__, integer *ldc, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int cgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	integer *kl, integer *ku, f2c_complex *ab, integer *ldab, f2c_real *d__, 
	f2c_real *e, f2c_complex *q, integer *ldq, f2c_complex *pt, integer *ldpt, 
	f2c_complex *c__, integer *ldc, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	f2c_complex *ab, integer *ldab, integer *ipiv, f2c_real *anorm, f2c_real *rcond, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	f2c_complex *ab, integer *ldab, f2c_real *r__, f2c_real *c__, f2c_real *rowcnd, f2c_real 
	*colcnd, f2c_real *amax, integer *info);

/* Subroutine */ int cgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, f2c_complex *ab, integer *ldab, f2c_complex *afb, integer *
	ldafb, integer *ipiv, f2c_complex *b, integer *ldb, f2c_complex *x, integer *
	ldx, f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, integer *
	info);

/* Subroutine */ int cgbsv_(integer *n, integer *kl, integer *ku, integer *
	nrhs, f2c_complex *ab, integer *ldab, integer *ipiv, f2c_complex *b, integer *
	ldb, integer *info);

/* Subroutine */ int cgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	integer *ku, integer *nrhs, f2c_complex *ab, integer *ldab, f2c_complex *afb,
	integer *ldafb, integer *ipiv, char *equed, f2c_real *r__, f2c_real *c__, 
	f2c_complex *b, integer *ldb, f2c_complex *x, integer *ldx, f2c_real *rcond, f2c_real 
	*ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	f2c_complex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int cgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	f2c_complex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int cgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, f2c_complex *ab, integer *ldab, integer *ipiv, f2c_complex 
	*b, integer *ldb, integer *info);

/* Subroutine */ int cgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, f2c_real *scale, integer *m, f2c_complex *v, integer *ldv, 
	integer *info);

/* Subroutine */ int cgebal_(char *job, integer *n, f2c_complex *a, integer *lda, 
	integer *ilo, integer *ihi, f2c_real *scale, integer *info);

/* Subroutine */ int cgebd2_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_real *d__, f2c_real *e, f2c_complex *tauq, f2c_complex *taup, f2c_complex *work, 
	integer *info);

/* Subroutine */ int cgebrd_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_real *d__, f2c_real *e, f2c_complex *tauq, f2c_complex *taup, f2c_complex *work, 
	integer *lwork, integer *info);

/* Subroutine */ int cgecon_(char *norm, integer *n, f2c_complex *a, integer *lda,
	f2c_real *anorm, f2c_real *rcond, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cgeequ_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_real *r__, f2c_real *c__, f2c_real *rowcnd, f2c_real *colcnd, f2c_real *amax, 
	integer *info);

/* Subroutine */ int cgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	f2c_complex *a, integer *lda, integer *sdim, f2c_complex *w, f2c_complex *vs, 
	integer *ldvs, f2c_complex *work, integer *lwork, f2c_real *rwork, logical *
	bwork, integer *info);

/* Subroutine */ int cgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, f2c_complex *a, integer *lda, integer *sdim, f2c_complex *
	w, f2c_complex *vs, integer *ldvs, f2c_real *rconde, f2c_real *rcondv, f2c_complex *
	work, integer *lwork, f2c_real *rwork, logical *bwork, integer *info);

/* Subroutine */ int cgeev_(char *jobvl, char *jobvr, integer *n, f2c_complex *a, 
	integer *lda, f2c_complex *w, f2c_complex *vl, integer *ldvl, f2c_complex *vr, 
	integer *ldvr, f2c_complex *work, integer *lwork, f2c_real *rwork, integer *
	info);

/* Subroutine */ int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, f2c_complex *a, integer *lda, f2c_complex *w, f2c_complex *vl, 
	integer *ldvl, f2c_complex *vr, integer *ldvr, integer *ilo, integer *ihi,
	f2c_real *scale, f2c_real *abnrm, f2c_real *rconde, f2c_real *rcondv, f2c_complex *work, 
	integer *lwork, f2c_real *rwork, integer *info);

/* Subroutine */ int cgegs_(char *jobvsl, char *jobvsr, integer *n, f2c_complex *
	a, integer *lda, f2c_complex *b, integer *ldb, f2c_complex *alpha, f2c_complex *
	beta, f2c_complex *vsl, integer *ldvsl, f2c_complex *vsr, integer *ldvsr, 
	f2c_complex *work, integer *lwork, f2c_real *rwork, integer *info);

/* Subroutine */ int cgegv_(char *jobvl, char *jobvr, integer *n, f2c_complex *a, 
	integer *lda, f2c_complex *b, integer *ldb, f2c_complex *alpha, f2c_complex *beta,
	f2c_complex *vl, integer *ldvl, f2c_complex *vr, integer *ldvr, f2c_complex *
	work, integer *lwork, f2c_real *rwork, integer *info);

/* Subroutine */ int cgehd2_(integer *n, integer *ilo, integer *ihi, f2c_complex *
	a, integer *lda, f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cgehrd_(integer *n, integer *ilo, integer *ihi, f2c_complex *
	a, integer *lda, f2c_complex *tau, f2c_complex *work, integer *lwork, integer 
	*info);

/* Subroutine */ int cgelq2_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cgelqf_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cgels_(char *trans, integer *m, integer *n, integer *
	nrhs, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, f2c_complex *
	work, integer *lwork, integer *info);

/* Subroutine */ int cgelsx_(integer *m, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *b, integer *ldb, integer *jpvt, f2c_real *rcond,
	integer *rank, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cgelsy_(integer *m, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *b, integer *ldb, integer *jpvt, f2c_real *rcond,
	integer *rank, f2c_complex *work, integer *lwork, f2c_real *rwork, integer *
	info);

/* Subroutine */ int cgeql2_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cgeqlf_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cgeqp3_(integer *m, integer *n, f2c_complex *a, integer *lda,
	integer *jpvt, f2c_complex *tau, f2c_complex *work, integer *lwork, f2c_real *
	rwork, integer *info);

/* Subroutine */ int cgeqpf_(integer *m, integer *n, f2c_complex *a, integer *lda,
	integer *jpvt, f2c_complex *tau, f2c_complex *work, f2c_real *rwork, integer *
	info);

/* Subroutine */ int cgeqr2_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cgeqrf_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cgerfs_(char *trans, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *af, integer *ldaf, integer *ipiv, f2c_complex *
	b, integer *ldb, f2c_complex *x, integer *ldx, f2c_real *ferr, f2c_real *berr, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cgerq2_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cgerqf_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cgesc2_(integer *n, f2c_complex *a, integer *lda, f2c_complex *
	rhs, integer *ipiv, integer *jpiv, f2c_real *scale);

/* Subroutine */ int cgesv_(integer *n, integer *nrhs, f2c_complex *a, integer *
	lda, integer *ipiv, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, f2c_complex *a, integer *lda, f2c_complex *af, integer *ldaf, integer *
	ipiv, char *equed, f2c_real *r__, f2c_real *c__, f2c_complex *b, integer *ldb, 
	f2c_complex *x, integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cgetc2_(integer *n, f2c_complex *a, integer *lda, integer *
	ipiv, integer *jpiv, integer *info);

/* Subroutine */ int cgetf2_(integer *m, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, integer *info);

/* Subroutine */ int cgetrf_(integer *m, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, integer *info);

/* Subroutine */ int cgetri_(integer *n, f2c_complex *a, integer *lda, integer *
	ipiv, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cgetrs_(char *trans, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, integer *ipiv, f2c_complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int cggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, f2c_real *lscale, f2c_real *rscale, integer *m, f2c_complex *v, 
	integer *ldv, integer *info);

/* Subroutine */ int cggbal_(char *job, integer *n, f2c_complex *a, integer *lda, 
	f2c_complex *b, integer *ldb, integer *ilo, integer *ihi, f2c_real *lscale, 
	f2c_real *rscale, f2c_real *work, integer *info);

/* Subroutine */ int cgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, integer *n, f2c_complex *a, integer *lda, f2c_complex *b, integer *
	ldb, integer *sdim, f2c_complex *alpha, f2c_complex *beta, f2c_complex *vsl, 
	integer *ldvsl, f2c_complex *vsr, integer *ldvsr, f2c_complex *work, integer *
	lwork, f2c_real *rwork, logical *bwork, integer *info);

/* Subroutine */ int cggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, integer *n, f2c_complex *a, integer *lda, f2c_complex *b,
	integer *ldb, integer *sdim, f2c_complex *alpha, f2c_complex *beta, f2c_complex *
	vsl, integer *ldvsl, f2c_complex *vsr, integer *ldvsr, f2c_real *rconde, f2c_real 
	*rcondv, f2c_complex *work, integer *lwork, f2c_real *rwork, integer *iwork, 
	integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int cggev_(char *jobvl, char *jobvr, integer *n, f2c_complex *a, 
	integer *lda, f2c_complex *b, integer *ldb, f2c_complex *alpha, f2c_complex *beta,
	f2c_complex *vl, integer *ldvl, f2c_complex *vr, integer *ldvr, f2c_complex *
	work, integer *lwork, f2c_real *rwork, integer *info);

/* Subroutine */ int cggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb,
	f2c_complex *alpha, f2c_complex *beta, f2c_complex *vl, integer *ldvl, f2c_complex *
	vr, integer *ldvr, integer *ilo, integer *ihi, f2c_real *lscale, f2c_real *
	rscale, f2c_real *abnrm, f2c_real *bbnrm, f2c_real *rconde, f2c_real *rcondv, f2c_complex 
	*work, integer *lwork, f2c_real *rwork, integer *iwork, logical *bwork, 
	integer *info);

/* Subroutine */ int cggglm_(integer *n, integer *m, integer *p, f2c_complex *a, 
	integer *lda, f2c_complex *b, integer *ldb, f2c_complex *d__, f2c_complex *x, 
	f2c_complex *y, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb,
	f2c_complex *q, integer *ldq, f2c_complex *z__, integer *ldz, integer *info);

/* Subroutine */ int cgglse_(integer *m, integer *n, integer *p, f2c_complex *a, 
	integer *lda, f2c_complex *b, integer *ldb, f2c_complex *c__, f2c_complex *d__, 
	f2c_complex *x, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cggqrf_(integer *n, integer *m, integer *p, f2c_complex *a, 
	integer *lda, f2c_complex *taua, f2c_complex *b, integer *ldb, f2c_complex *taub, 
	f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cggrqf_(integer *m, integer *p, integer *n, f2c_complex *a, 
	integer *lda, f2c_complex *taua, f2c_complex *b, integer *ldb, f2c_complex *taub, 
	f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, f2c_complex *a, integer *
	lda, f2c_complex *b, integer *ldb, f2c_real *alpha, f2c_real *beta, f2c_complex *u, 
	integer *ldu, f2c_complex *v, integer *ldv, f2c_complex *q, integer *ldq, 
	f2c_complex *work, f2c_real *rwork, integer *iwork, integer *info);

/* Subroutine */ int cggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, f2c_complex *a, integer *lda, f2c_complex *b, integer 
	*ldb, f2c_real *tola, f2c_real *tolb, integer *k, integer *l, f2c_complex *u, 
	integer *ldu, f2c_complex *v, integer *ldv, f2c_complex *q, integer *ldq, 
	integer *iwork, f2c_real *rwork, f2c_complex *tau, f2c_complex *work, integer *
	info);

/* Subroutine */ int cgtcon_(char *norm, integer *n, f2c_complex *dl, f2c_complex *
	d__, f2c_complex *du, f2c_complex *du2, integer *ipiv, f2c_real *anorm, f2c_real *
	rcond, f2c_complex *work, integer *info);

/* Subroutine */ int cgtrfs_(char *trans, integer *n, integer *nrhs, f2c_complex *
	dl, f2c_complex *d__, f2c_complex *du, f2c_complex *dlf, f2c_complex *df, f2c_complex *
	duf, f2c_complex *du2, integer *ipiv, f2c_complex *b, integer *ldb, f2c_complex *
	x, integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int cgtsv_(integer *n, integer *nrhs, f2c_complex *dl, f2c_complex *
	d__, f2c_complex *du, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, f2c_complex *dl, f2c_complex *d__, f2c_complex *du, f2c_complex *dlf, f2c_complex *
	df, f2c_complex *duf, f2c_complex *du2, integer *ipiv, f2c_complex *b, integer *
	ldb, f2c_complex *x, integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cgttrf_(integer *n, f2c_complex *dl, f2c_complex *d__, f2c_complex *
	du, f2c_complex *du2, integer *ipiv, integer *info);

/* Subroutine */ int cgttrs_(char *trans, integer *n, integer *nrhs, f2c_complex *
	dl, f2c_complex *d__, f2c_complex *du, f2c_complex *du2, integer *ipiv, f2c_complex *
	b, integer *ldb, integer *info);

/* Subroutine */ int cgtts2_(integer *itrans, integer *n, integer *nrhs, 
	f2c_complex *dl, f2c_complex *d__, f2c_complex *du, f2c_complex *du2, integer *ipiv, 
	f2c_complex *b, integer *ldb);

/* Subroutine */ int chbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	f2c_complex *ab, integer *ldab, f2c_real *w, f2c_complex *z__, integer *ldz, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int chbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	f2c_complex *ab, integer *ldab, f2c_real *w, f2c_complex *z__, integer *ldz, 
	f2c_complex *work, integer *lwork, f2c_real *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int chbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, f2c_complex *ab, integer *ldab, f2c_complex *q, integer *ldq, 
	f2c_real *vl, f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, integer *
	m, f2c_real *w, f2c_complex *z__, integer *ldz, f2c_complex *work, f2c_real *rwork, 
	integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int chbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, f2c_complex *ab, integer *ldab, f2c_complex *bb, integer *ldbb, 
	f2c_complex *x, integer *ldx, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int chbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, f2c_complex *ab, integer *ldab, f2c_complex *bb, integer *ldbb, 
	f2c_real *w, f2c_complex *z__, integer *ldz, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int chbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, f2c_complex *ab, integer *ldab, f2c_complex *bb, 
	integer *ldbb, f2c_complex *q, integer *ldq, f2c_real *vl, f2c_real *vu, integer *
	il, integer *iu, f2c_real *abstol, integer *m, f2c_real *w, f2c_complex *z__, 
	integer *ldz, f2c_complex *work, f2c_real *rwork, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int chbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	f2c_complex *ab, integer *ldab, f2c_real *d__, f2c_real *e, f2c_complex *q, integer *
	ldq, f2c_complex *work, integer *info);

/* Subroutine */ int checon_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, f2c_real *anorm, f2c_real *rcond, f2c_complex *work, integer *
	info);

/* Subroutine */ int cheev_(char *jobz, char *uplo, integer *n, f2c_complex *a, 
	integer *lda, f2c_real *w, f2c_complex *work, integer *lwork, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int cheevd_(char *jobz, char *uplo, integer *n, f2c_complex *a, 
	integer *lda, f2c_real *w, f2c_complex *work, integer *lwork, f2c_real *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int cheevr_(char *jobz, char *range, char *uplo, integer *n, 
	f2c_complex *a, integer *lda, f2c_real *vl, f2c_real *vu, integer *il, integer *
	iu, f2c_real *abstol, integer *m, f2c_real *w, f2c_complex *z__, integer *ldz, 
	integer *isuppz, f2c_complex *work, integer *lwork, f2c_real *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int cheevx_(char *jobz, char *range, char *uplo, integer *n, 
	f2c_complex *a, integer *lda, f2c_real *vl, f2c_real *vu, integer *il, integer *
	iu, f2c_real *abstol, integer *m, f2c_real *w, f2c_complex *z__, integer *ldz, 
	f2c_complex *work, integer *lwork, f2c_real *rwork, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int chegs2_(integer *itype, char *uplo, integer *n, f2c_complex *
	a, integer *lda, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int chegst_(integer *itype, char *uplo, integer *n, f2c_complex *
	a, integer *lda, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int chegv_(integer *itype, char *jobz, char *uplo, integer *
	n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, f2c_real *w, 
	f2c_complex *work, integer *lwork, f2c_real *rwork, integer *info);

/* Subroutine */ int chegvd_(integer *itype, char *jobz, char *uplo, integer *
	n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, f2c_real *w, 
	f2c_complex *work, integer *lwork, f2c_real *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int chegvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, 
	f2c_real *vl, f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, integer *
	m, f2c_real *w, f2c_complex *z__, integer *ldz, f2c_complex *work, integer *lwork,
	f2c_real *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int cherfs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *af, integer *ldaf, integer *ipiv, f2c_complex *
	b, integer *ldb, f2c_complex *x, integer *ldx, f2c_real *ferr, f2c_real *berr, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int chesv_(char *uplo, integer *n, integer *nrhs, f2c_complex *a,
	integer *lda, integer *ipiv, f2c_complex *b, integer *ldb, f2c_complex *work,
	integer *lwork, integer *info);

/* Subroutine */ int chesvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_complex *a, integer *lda, f2c_complex *af, integer *ldaf, integer *
	ipiv, f2c_complex *b, integer *ldb, f2c_complex *x, integer *ldx, f2c_real *rcond,
	f2c_real *ferr, f2c_real *berr, f2c_complex *work, integer *lwork, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int chetf2_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, integer *info);

/* Subroutine */ int chetrd_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	f2c_real *d__, f2c_real *e, f2c_complex *tau, f2c_complex *work, integer *lwork, 
	integer *info);

/* Subroutine */ int chetrf_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int chetri_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, f2c_complex *work, integer *info);

/* Subroutine */ int chetrs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, integer *ipiv, f2c_complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int chgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, f2c_complex *a, integer *lda, f2c_complex *b, 
	integer *ldb, f2c_complex *alpha, f2c_complex *beta, f2c_complex *q, integer *ldq,
	f2c_complex *z__, integer *ldz, f2c_complex *work, integer *lwork, f2c_real *
	rwork, integer *info);

/* Subroutine */ int chpcon_(char *uplo, integer *n, f2c_complex *ap, integer *
	ipiv, f2c_real *anorm, f2c_real *rcond, f2c_complex *work, integer *info);

/* Subroutine */ int chpev_(char *jobz, char *uplo, integer *n, f2c_complex *ap, 
	f2c_real *w, f2c_complex *z__, integer *ldz, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int chpevd_(char *jobz, char *uplo, integer *n, f2c_complex *ap, 
	f2c_real *w, f2c_complex *z__, integer *ldz, f2c_complex *work, integer *lwork, 
	f2c_real *rwork, integer *lrwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int chpevx_(char *jobz, char *range, char *uplo, integer *n, 
	f2c_complex *ap, f2c_real *vl, f2c_real *vu, integer *il, integer *iu, f2c_real *
	abstol, integer *m, f2c_real *w, f2c_complex *z__, integer *ldz, f2c_complex *
	work, f2c_real *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int chpgst_(integer *itype, char *uplo, integer *n, f2c_complex *
	ap, f2c_complex *bp, integer *info);

/* Subroutine */ int chpgv_(integer *itype, char *jobz, char *uplo, integer *
	n, f2c_complex *ap, f2c_complex *bp, f2c_real *w, f2c_complex *z__, integer *ldz, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int chpgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, f2c_complex *ap, f2c_complex *bp, f2c_real *w, f2c_complex *z__, integer *ldz, 
	f2c_complex *work, integer *lwork, f2c_real *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int chpgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, f2c_complex *ap, f2c_complex *bp, f2c_real *vl, f2c_real *vu, 
	integer *il, integer *iu, f2c_real *abstol, integer *m, f2c_real *w, f2c_complex *
	z__, integer *ldz, f2c_complex *work, f2c_real *rwork, integer *iwork, 
	integer *ifail, integer *info);

/* Subroutine */ int chprfs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, f2c_complex *afp, integer *ipiv, f2c_complex *b, integer *ldb, f2c_complex *x,
	integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int chpsv_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, integer *ipiv, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int chpsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_complex *ap, f2c_complex *afp, integer *ipiv, f2c_complex *b, integer *
	ldb, f2c_complex *x, integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int chptrd_(char *uplo, integer *n, f2c_complex *ap, f2c_real *d__, 
	f2c_real *e, f2c_complex *tau, integer *info);

/* Subroutine */ int chptrf_(char *uplo, integer *n, f2c_complex *ap, integer *
	ipiv, integer *info);

/* Subroutine */ int chptri_(char *uplo, integer *n, f2c_complex *ap, integer *
	ipiv, f2c_complex *work, integer *info);

/* Subroutine */ int chptrs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, integer *ipiv, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int chsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, f2c_complex *h__, integer *ldh, f2c_complex *w, f2c_complex *
	vl, integer *ldvl, f2c_complex *vr, integer *ldvr, integer *mm, integer *
	m, f2c_complex *work, f2c_real *rwork, integer *ifaill, integer *ifailr, 
	integer *info);

/* Subroutine */ int chseqr_(char *job, char *compz, integer *n, integer *ilo,
	integer *ihi, f2c_complex *h__, integer *ldh, f2c_complex *w, f2c_complex *z__, 
	integer *ldz, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int clabrd_(integer *m, integer *n, integer *nb, f2c_complex *a, 
	integer *lda, f2c_real *d__, f2c_real *e, f2c_complex *tauq, f2c_complex *taup, 
	f2c_complex *x, integer *ldx, f2c_complex *y, integer *ldy);

/* Subroutine */ int clacgv_(integer *n, f2c_complex *x, integer *incx);

/* Subroutine */ int clacon_(integer *n, f2c_complex *v, f2c_complex *x, f2c_real *est, 
	integer *kase);

/* Subroutine */ int clacp2_(char *uplo, integer *m, integer *n, f2c_real *a, 
	integer *lda, f2c_complex *b, integer *ldb);

/* Subroutine */ int clacpy_(char *uplo, integer *m, integer *n, f2c_complex *a, 
	integer *lda, f2c_complex *b, integer *ldb);

/* Subroutine */ int clacrm_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_real *b, integer *ldb, f2c_complex *c__, integer *ldc, f2c_real *rwork);

/* Subroutine */ int clacrt_(integer *n, f2c_complex *cx, integer *incx, f2c_complex *
	cy, integer *incy, f2c_complex *c__, f2c_complex *s);

/* Subroutine */ int claed0_(integer *qsiz, integer *n, f2c_real *d__, f2c_real *e, 
	f2c_complex *q, integer *ldq, f2c_complex *qstore, integer *ldqs, f2c_real *rwork,
	integer *iwork, integer *info);

/* Subroutine */ int claed7_(integer *n, integer *cutpnt, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, f2c_real *d__, f2c_complex *
	q, integer *ldq, f2c_real *rho, integer *indxq, f2c_real *qstore, integer *
	qptr, integer *prmptr, integer *perm, integer *givptr, integer *
	givcol, f2c_real *givnum, f2c_complex *work, f2c_real *rwork, integer *iwork, 
	integer *info);

/* Subroutine */ int claed8_(integer *k, integer *n, integer *qsiz, f2c_complex *
	q, integer *ldq, f2c_real *d__, f2c_real *rho, integer *cutpnt, f2c_real *z__, 
	f2c_real *dlamda, f2c_complex *q2, integer *ldq2, f2c_real *w, integer *indxp, 
	integer *indx, integer *indxq, integer *perm, integer *givptr, 
	integer *givcol, f2c_real *givnum, integer *info);

/* Subroutine */ int claein_(logical *rightv, logical *noinit, integer *n, 
	f2c_complex *h__, integer *ldh, f2c_complex *w, f2c_complex *v, f2c_complex *b, 
	integer *ldb, f2c_real *rwork, f2c_real *eps3, f2c_real *smlnum, integer *info);

/* Subroutine */ int claesy_(f2c_complex *a, f2c_complex *b, f2c_complex *c__, f2c_complex *
	rt1, f2c_complex *rt2, f2c_complex *evscal, f2c_complex *cs1, f2c_complex *sn1);

/* Subroutine */ int claev2_(f2c_complex *a, f2c_complex *b, f2c_complex *c__, f2c_real *rt1, 
	f2c_real *rt2, f2c_real *cs1, f2c_complex *sn1);

/* Subroutine */ int clags2_(logical *upper, f2c_real *a1, f2c_complex *a2, f2c_real *a3, 
	f2c_real *b1, f2c_complex *b2, f2c_real *b3, f2c_real *csu, f2c_complex *snu, f2c_real *csv, 
	f2c_complex *snv, f2c_real *csq, f2c_complex *snq);

/* Subroutine */ int clagtm_(char *trans, integer *n, integer *nrhs, f2c_real *
	alpha, f2c_complex *dl, f2c_complex *d__, f2c_complex *du, f2c_complex *x, integer *
	ldx, f2c_real *beta, f2c_complex *b, integer *ldb);

/* Subroutine */ int clahef_(char *uplo, integer *n, integer *nb, integer *kb,
	f2c_complex *a, integer *lda, integer *ipiv, f2c_complex *w, integer *ldw, 
	integer *info);

/* Subroutine */ int clahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, f2c_complex *h__, integer *ldh, f2c_complex *w, 
	integer *iloz, integer *ihiz, f2c_complex *z__, integer *ldz, integer *
	info);

/* Subroutine */ int clahrd_(integer *n, integer *k, integer *nb, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *t, integer *ldt, f2c_complex *y, 
	integer *ldy);

/* Subroutine */ int claic1_(integer *job, integer *j, f2c_complex *x, f2c_real *sest,
	f2c_complex *w, f2c_complex *gamma, f2c_real *sestpr, f2c_complex *s, f2c_complex *c__);

/* Subroutine */ int clals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, f2c_complex *b, integer *ldb, f2c_complex *bx, 
	integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, f2c_real *givnum, integer *ldgnum, f2c_real *poles, f2c_real *
	difl, f2c_real *difr, f2c_real *z__, integer *k, f2c_real *c__, f2c_real *s, f2c_real *
	rwork, integer *info);

/* Subroutine */ int clalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, f2c_complex *b, integer *ldb, f2c_complex *bx, integer *ldbx, 
	f2c_real *u, integer *ldu, f2c_real *vt, integer *k, f2c_real *difl, f2c_real *difr, 
	f2c_real *z__, f2c_real *poles, integer *givptr, integer *givcol, integer *
	ldgcol, integer *perm, f2c_real *givnum, f2c_real *c__, f2c_real *s, f2c_real *rwork, 
	integer *iwork, integer *info);

/* Subroutine */ int clapll_(integer *n, f2c_complex *x, integer *incx, f2c_complex *
	y, integer *incy, f2c_real *ssmin);

/* Subroutine */ int clapmt_(logical *forwrd, integer *m, integer *n, f2c_complex 
	*x, integer *ldx, integer *k);

/* Subroutine */ int claqgb_(integer *m, integer *n, integer *kl, integer *ku,
	f2c_complex *ab, integer *ldab, f2c_real *r__, f2c_real *c__, f2c_real *rowcnd, f2c_real 
	*colcnd, f2c_real *amax, char *equed);

/* Subroutine */ int claqge_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_real *r__, f2c_real *c__, f2c_real *rowcnd, f2c_real *colcnd, f2c_real *amax, char *
	equed);

/* Subroutine */ int claqhb_(char *uplo, integer *n, integer *kd, f2c_complex *ab,
	integer *ldab, f2c_real *s, f2c_real *scond, f2c_real *amax, char *equed);

/* Subroutine */ int claqhe_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	f2c_real *s, f2c_real *scond, f2c_real *amax, char *equed);

/* Subroutine */ int claqhp_(char *uplo, integer *n, f2c_complex *ap, f2c_real *s, 
	f2c_real *scond, f2c_real *amax, char *equed);

/* Subroutine */ int claqp2_(integer *m, integer *n, integer *offset, f2c_complex 
	*a, integer *lda, integer *jpvt, f2c_complex *tau, f2c_real *vn1, f2c_real *vn2, 
	f2c_complex *work);

/* Subroutine */ int claqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, f2c_complex *a, integer *lda, integer *jpvt, f2c_complex *
	tau, f2c_real *vn1, f2c_real *vn2, f2c_complex *auxv, f2c_complex *f, integer *ldf);

/* Subroutine */ int claqsb_(char *uplo, integer *n, integer *kd, f2c_complex *ab,
	integer *ldab, f2c_real *s, f2c_real *scond, f2c_real *amax, char *equed);

/* Subroutine */ int claqsp_(char *uplo, integer *n, f2c_complex *ap, f2c_real *s, 
	f2c_real *scond, f2c_real *amax, char *equed);

/* Subroutine */ int claqsy_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	f2c_real *s, f2c_real *scond, f2c_real *amax, char *equed);

/* Subroutine */ int clar1v_(integer *n, integer *b1, integer *bn, f2c_real *
	sigma, f2c_real *d__, f2c_real *l, f2c_real *ld, f2c_real *lld, f2c_real *gersch, f2c_complex 
	*z__, f2c_real *ztz, f2c_real *mingma, integer *r__, integer *isuppz, f2c_real *
	work);

/* Subroutine */ int clar2v_(integer *n, f2c_complex *x, f2c_complex *y, f2c_complex *z__,
	integer *incx, f2c_real *c__, f2c_complex *s, integer *incc);

/* Subroutine */ int clarcm_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_complex *b, integer *ldb, f2c_complex *c__, integer *ldc, f2c_real *rwork);

/* Subroutine */ int clarf_(char *side, integer *m, integer *n, f2c_complex *v, 
	integer *incv, f2c_complex *tau, f2c_complex *c__, integer *ldc, f2c_complex *
	work);

/* Subroutine */ int clarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, f2c_complex *v, integer *ldv, 
	f2c_complex *t, integer *ldt, f2c_complex *c__, integer *ldc, f2c_complex *work, 
	integer *ldwork);

/* Subroutine */ int clarfg_(integer *n, f2c_complex *alpha, f2c_complex *x, integer *
	incx, f2c_complex *tau);

/* Subroutine */ int clarft_(char *direct, char *storev, integer *n, integer *
	k, f2c_complex *v, integer *ldv, f2c_complex *tau, f2c_complex *t, integer *ldt);

/* Subroutine */ int clarfx_(char *side, integer *m, integer *n, f2c_complex *v, 
	f2c_complex *tau, f2c_complex *c__, integer *ldc, f2c_complex *work);

/* Subroutine */ int clargv_(integer *n, f2c_complex *x, integer *incx, f2c_complex *
	y, integer *incy, f2c_real *c__, integer *incc);

/* Subroutine */ int clarnv_(integer *idist, integer *iseed, integer *n, 
	f2c_complex *x);

/* Subroutine */ int clarrv_(integer *n, f2c_real *d__, f2c_real *l, integer *isplit, 
	integer *m, f2c_real *w, integer *iblock, f2c_real *gersch, f2c_real *tol, 
	f2c_complex *z__, integer *ldz, integer *isuppz, f2c_real *work, integer *
	iwork, integer *info);

/* Subroutine */ int clartg_(f2c_complex *f, f2c_complex *g, f2c_real *cs, f2c_complex *sn, 
	f2c_complex *r__);

/* Subroutine */ int clartv_(integer *n, f2c_complex *x, integer *incx, f2c_complex *
	y, integer *incy, f2c_real *c__, f2c_complex *s, integer *incc);

/* Subroutine */ int clarz_(char *side, integer *m, integer *n, integer *l, 
	f2c_complex *v, integer *incv, f2c_complex *tau, f2c_complex *c__, integer *ldc, 
	f2c_complex *work);

/* Subroutine */ int clarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, f2c_complex *v, 
	integer *ldv, f2c_complex *t, integer *ldt, f2c_complex *c__, integer *ldc, 
	f2c_complex *work, integer *ldwork);

/* Subroutine */ int clarzt_(char *direct, char *storev, integer *n, integer *
	k, f2c_complex *v, integer *ldv, f2c_complex *tau, f2c_complex *t, integer *ldt);

/* Subroutine */ int clascl_(char *type__, integer *kl, integer *ku, f2c_real *
	cfrom, f2c_real *cto, integer *m, integer *n, f2c_complex *a, integer *lda, 
	integer *info);

/* Subroutine */ int claset_(char *uplo, integer *m, integer *n, f2c_complex *
	alpha, f2c_complex *beta, f2c_complex *a, integer *lda);

/* Subroutine */ int clasr_(char *side, char *pivot, char *direct, integer *m,
	integer *n, f2c_real *c__, f2c_real *s, f2c_complex *a, integer *lda);

/* Subroutine */ int classq_(integer *n, f2c_complex *x, integer *incx, f2c_real *
	scale, f2c_real *sumsq);

/* Subroutine */ int claswp_(integer *n, f2c_complex *a, integer *lda, integer *
	k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int clasyf_(char *uplo, integer *n, integer *nb, integer *kb,
	f2c_complex *a, integer *lda, integer *ipiv, f2c_complex *w, integer *ldw, 
	integer *info);

/* Subroutine */ int clatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, f2c_complex *ab, integer *ldab, f2c_complex *
	x, f2c_real *scale, f2c_real *cnorm, integer *info);

/* Subroutine */ int clatdf_(integer *ijob, integer *n, f2c_complex *z__, integer 
	*ldz, f2c_complex *rhs, f2c_real *rdsum, f2c_real *rdscal, integer *ipiv, integer 
	*jpiv);

/* Subroutine */ int clatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, f2c_complex *ap, f2c_complex *x, f2c_real *scale, f2c_real *cnorm,
	integer *info);

/* Subroutine */ int clatrd_(char *uplo, integer *n, integer *nb, f2c_complex *a, 
	integer *lda, f2c_real *e, f2c_complex *tau, f2c_complex *w, integer *ldw);

/* Subroutine */ int clatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, f2c_complex *a, integer *lda, f2c_complex *x, f2c_real *scale,
	f2c_real *cnorm, integer *info);

/* Subroutine */ int clatrz_(integer *m, integer *n, integer *l, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work);

/* Subroutine */ int clatzm_(char *side, integer *m, integer *n, f2c_complex *v, 
	integer *incv, f2c_complex *tau, f2c_complex *c1, f2c_complex *c2, integer *ldc, 
	f2c_complex *work);

/* Subroutine */ int clauu2_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *info);

/* Subroutine */ int clauum_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *info);

/* Subroutine */ int cpbcon_(char *uplo, integer *n, integer *kd, f2c_complex *ab,
	integer *ldab, f2c_real *anorm, f2c_real *rcond, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int cpbequ_(char *uplo, integer *n, integer *kd, f2c_complex *ab,
	integer *ldab, f2c_real *s, f2c_real *scond, f2c_real *amax, integer *info);

/* Subroutine */ int cpbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, f2c_complex *ab, integer *ldab, f2c_complex *afb, integer *ldafb, 
	f2c_complex *b, integer *ldb, f2c_complex *x, integer *ldx, f2c_real *ferr, f2c_real *
	berr, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cpbstf_(char *uplo, integer *n, integer *kd, f2c_complex *ab,
	integer *ldab, integer *info);

/* Subroutine */ int cpbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, f2c_complex *ab, integer *ldab, f2c_complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int cpbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, f2c_complex *ab, integer *ldab, f2c_complex *afb, integer *
	ldafb, char *equed, f2c_real *s, f2c_complex *b, integer *ldb, f2c_complex *x, 
	integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_complex *work, 
	f2c_real *rwork, integer *info);

/* Subroutine */ int cpbtf2_(char *uplo, integer *n, integer *kd, f2c_complex *ab,
	integer *ldab, integer *info);

/* Subroutine */ int cpbtrf_(char *uplo, integer *n, integer *kd, f2c_complex *ab,
	integer *ldab, integer *info);

/* Subroutine */ int cpbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, f2c_complex *ab, integer *ldab, f2c_complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int cpocon_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	f2c_real *anorm, f2c_real *rcond, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cpoequ_(integer *n, f2c_complex *a, integer *lda, f2c_real *s, 
	f2c_real *scond, f2c_real *amax, integer *info);

/* Subroutine */ int cporfs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *af, integer *ldaf, f2c_complex *b, integer *ldb,
	f2c_complex *x, integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_complex *work, 
	f2c_real *rwork, integer *info);

/* Subroutine */ int cposv_(char *uplo, integer *n, integer *nrhs, f2c_complex *a,
	integer *lda, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_complex *a, integer *lda, f2c_complex *af, integer *ldaf, char *
	equed, f2c_real *s, f2c_complex *b, integer *ldb, f2c_complex *x, integer *ldx, 
	f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int cpotf2_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *info);

/* Subroutine */ int cpotrf_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *info);

/* Subroutine */ int cpotri_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *info);

/* Subroutine */ int cpotrs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cppcon_(char *uplo, integer *n, f2c_complex *ap, f2c_real *anorm,
	f2c_real *rcond, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cppequ_(char *uplo, integer *n, f2c_complex *ap, f2c_real *s, 
	f2c_real *scond, f2c_real *amax, integer *info);

/* Subroutine */ int cpprfs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, f2c_complex *afp, f2c_complex *b, integer *ldb, f2c_complex *x, integer *ldx, 
	f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cppsv_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_complex *ap, f2c_complex *afp, char *equed, f2c_real *s, f2c_complex *b, 
	integer *ldb, f2c_complex *x, integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real 
	*berr, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int cpptrf_(char *uplo, integer *n, f2c_complex *ap, integer *
	info);

/* Subroutine */ int cpptri_(char *uplo, integer *n, f2c_complex *ap, integer *
	info);

/* Subroutine */ int cpptrs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cptcon_(integer *n, f2c_real *d__, f2c_complex *e, f2c_real *anorm, 
	f2c_real *rcond, f2c_real *rwork, integer *info);

/* Subroutine */ int cptrfs_(char *uplo, integer *n, integer *nrhs, f2c_real *d__,
	f2c_complex *e, f2c_real *df, f2c_complex *ef, f2c_complex *b, integer *ldb, f2c_complex 
	*x, integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int cptsv_(integer *n, integer *nrhs, f2c_real *d__, f2c_complex *e, 
	f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cptsvx_(char *fact, integer *n, integer *nrhs, f2c_real *d__,
	f2c_complex *e, f2c_real *df, f2c_complex *ef, f2c_complex *b, integer *ldb, f2c_complex 
	*x, integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_complex *work, 
	f2c_real *rwork, integer *info);

/* Subroutine */ int cpttrf_(integer *n, f2c_real *d__, f2c_complex *e, integer *info);

/* Subroutine */ int cpttrs_(char *uplo, integer *n, integer *nrhs, f2c_real *d__,
	f2c_complex *e, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cptts2_(integer *iuplo, integer *n, integer *nrhs, f2c_real *
	d__, f2c_complex *e, f2c_complex *b, integer *ldb);

/* Subroutine */ int crot_(integer *n, f2c_complex *cx, integer *incx, f2c_complex *
	cy, integer *incy, f2c_real *c__, f2c_complex *s);

/* Subroutine */ int cspcon_(char *uplo, integer *n, f2c_complex *ap, integer *
	ipiv, f2c_real *anorm, f2c_real *rcond, f2c_complex *work, integer *info);

/* Subroutine */ int cspmv_(char *uplo, integer *n, f2c_complex *alpha, f2c_complex *
	ap, f2c_complex *x, integer *incx, f2c_complex *beta, f2c_complex *y, integer *
	incy);

/* Subroutine */ int cspr_(char *uplo, integer *n, f2c_complex *alpha, f2c_complex *x,
	integer *incx, f2c_complex *ap);

/* Subroutine */ int csprfs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, f2c_complex *afp, integer *ipiv, f2c_complex *b, integer *ldb, f2c_complex *x,
	integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int cspsv_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, integer *ipiv, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int cspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_complex *ap, f2c_complex *afp, integer *ipiv, f2c_complex *b, integer *
	ldb, f2c_complex *x, integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int csptrf_(char *uplo, integer *n, f2c_complex *ap, integer *
	ipiv, integer *info);

/* Subroutine */ int csptri_(char *uplo, integer *n, f2c_complex *ap, integer *
	ipiv, f2c_complex *work, integer *info);

/* Subroutine */ int csptrs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	ap, integer *ipiv, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int csrot_(integer *n, f2c_complex *cx, integer *incx, f2c_complex *
	cy, integer *incy, f2c_real *c__, f2c_real *s);

/* Subroutine */ int csrscl_(integer *n, f2c_real *sa, f2c_complex *sx, integer *incx);

/* Subroutine */ int cstedc_(char *compz, integer *n, f2c_real *d__, f2c_real *e, 
	f2c_complex *z__, integer *ldz, f2c_complex *work, integer *lwork, f2c_real *
	rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
	info);

/* Subroutine */ int cstein_(integer *n, f2c_real *d__, f2c_real *e, integer *m, f2c_real 
	*w, integer *iblock, integer *isplit, f2c_complex *z__, integer *ldz, 
	f2c_real *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int csteqr_(char *compz, integer *n, f2c_real *d__, f2c_real *e, 
	f2c_complex *z__, integer *ldz, f2c_real *work, integer *info);

/* Subroutine */ int csycon_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, f2c_real *anorm, f2c_real *rcond, f2c_complex *work, integer *
	info);

/* Subroutine */ int csymv_(char *uplo, integer *n, f2c_complex *alpha, f2c_complex *
	a, integer *lda, f2c_complex *x, integer *incx, f2c_complex *beta, f2c_complex *y,
	integer *incy);

/* Subroutine */ int csyr_(char *uplo, integer *n, f2c_complex *alpha, f2c_complex *x,
	integer *incx, f2c_complex *a, integer *lda);

/* Subroutine */ int csyrfs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *af, integer *ldaf, integer *ipiv, f2c_complex *
	b, integer *ldb, f2c_complex *x, integer *ldx, f2c_real *ferr, f2c_real *berr, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int csysv_(char *uplo, integer *n, integer *nrhs, f2c_complex *a,
	integer *lda, integer *ipiv, f2c_complex *b, integer *ldb, f2c_complex *work,
	integer *lwork, integer *info);

/* Subroutine */ int csysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_complex *a, integer *lda, f2c_complex *af, integer *ldaf, integer *
	ipiv, f2c_complex *b, integer *ldb, f2c_complex *x, integer *ldx, f2c_real *rcond,
	f2c_real *ferr, f2c_real *berr, f2c_complex *work, integer *lwork, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int csytf2_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, integer *info);

/* Subroutine */ int csytrf_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int csytri_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	integer *ipiv, f2c_complex *work, integer *info);

/* Subroutine */ int csytrs_(char *uplo, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, integer *ipiv, f2c_complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int ctbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, f2c_complex *ab, integer *ldab, f2c_real *rcond, f2c_complex *work, 
	f2c_real *rwork, integer *info);

/* Subroutine */ int ctbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, f2c_complex *ab, integer *ldab, f2c_complex *b, 
	integer *ldb, f2c_complex *x, integer *ldx, f2c_real *ferr, f2c_real *berr, 
	f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int ctbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, f2c_complex *ab, integer *ldab, f2c_complex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int ctgevc_(char *side, char *howmny, logical *select, 
	integer *n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, 
	f2c_complex *vl, integer *ldvl, f2c_complex *vr, integer *ldvr, integer *mm, 
	integer *m, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int ctgex2_(logical *wantq, logical *wantz, integer *n, 
	f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, f2c_complex *q, 
	integer *ldq, f2c_complex *z__, integer *ldz, integer *j1, integer *info);

/* Subroutine */ int ctgexc_(logical *wantq, logical *wantz, integer *n, 
	f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, f2c_complex *q, 
	integer *ldq, f2c_complex *z__, integer *ldz, integer *ifst, integer *
	ilst, integer *info);

/* Subroutine */ int ctgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, f2c_complex *a, integer *lda, f2c_complex *b, 
	integer *ldb, f2c_complex *alpha, f2c_complex *beta, f2c_complex *q, integer *ldq,
	f2c_complex *z__, integer *ldz, integer *m, f2c_real *pl, f2c_real *pr, f2c_real *
	dif, f2c_complex *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int ctgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, f2c_complex *a, integer *
	lda, f2c_complex *b, integer *ldb, f2c_real *tola, f2c_real *tolb, f2c_real *alpha, 
	f2c_real *beta, f2c_complex *u, integer *ldu, f2c_complex *v, integer *ldv, 
	f2c_complex *q, integer *ldq, f2c_complex *work, integer *ncycle, integer *
	info);

/* Subroutine */ int ctgsna_(char *job, char *howmny, logical *select, 
	integer *n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, 
	f2c_complex *vl, integer *ldvl, f2c_complex *vr, integer *ldvr, f2c_real *s, f2c_real 
	*dif, integer *mm, integer *m, f2c_complex *work, integer *lwork, integer 
	*iwork, integer *info);

/* Subroutine */ int ctgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, f2c_complex *c__, 
	integer *ldc, f2c_complex *d__, integer *ldd, f2c_complex *e, integer *lde, 
	f2c_complex *f, integer *ldf, f2c_real *scale, f2c_real *rdsum, f2c_real *rdscal, 
	integer *info);

/* Subroutine */ int ctgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, f2c_complex *c__, 
	integer *ldc, f2c_complex *d__, integer *ldd, f2c_complex *e, integer *lde, 
	f2c_complex *f, integer *ldf, f2c_real *scale, f2c_real *dif, f2c_complex *work, 
	integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int ctpcon_(char *norm, char *uplo, char *diag, integer *n, 
	f2c_complex *ap, f2c_real *rcond, f2c_complex *work, f2c_real *rwork, integer *info);

/* Subroutine */ int ctprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, f2c_complex *ap, f2c_complex *b, integer *ldb, f2c_complex *x, 
	integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int ctptri_(char *uplo, char *diag, integer *n, f2c_complex *ap, 
	integer *info);

/* Subroutine */ int ctptrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, f2c_complex *ap, f2c_complex *b, integer *ldb, integer *info);

/* Subroutine */ int ctrcon_(char *norm, char *uplo, char *diag, integer *n, 
	f2c_complex *a, integer *lda, f2c_real *rcond, f2c_complex *work, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int ctrevc_(char *side, char *howmny, logical *select, 
	integer *n, f2c_complex *t, integer *ldt, f2c_complex *vl, integer *ldvl, 
	f2c_complex *vr, integer *ldvr, integer *mm, integer *m, f2c_complex *work, 
	f2c_real *rwork, integer *info);

/* Subroutine */ int ctrexc_(char *compq, integer *n, f2c_complex *t, integer *
	ldt, f2c_complex *q, integer *ldq, integer *ifst, integer *ilst, integer *
	info);

/* Subroutine */ int ctrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, 
	f2c_complex *x, integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_complex *work, f2c_real 
	*rwork, integer *info);

/* Subroutine */ int ctrsen_(char *job, char *compq, logical *select, integer 
	*n, f2c_complex *t, integer *ldt, f2c_complex *q, integer *ldq, f2c_complex *w, 
	integer *m, f2c_real *s, f2c_real *sep, f2c_complex *work, integer *lwork, 
	integer *info);

/* Subroutine */ int ctrsna_(char *job, char *howmny, logical *select, 
	integer *n, f2c_complex *t, integer *ldt, f2c_complex *vl, integer *ldvl, 
	f2c_complex *vr, integer *ldvr, f2c_real *s, f2c_real *sep, integer *mm, integer *
	m, f2c_complex *work, integer *ldwork, f2c_real *rwork, integer *info);

/* Subroutine */ int ctrsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, 
	f2c_complex *c__, integer *ldc, f2c_real *scale, integer *info);

/* Subroutine */ int ctrti2_(char *uplo, char *diag, integer *n, f2c_complex *a, 
	integer *lda, integer *info);

/* Subroutine */ int ctrtri_(char *uplo, char *diag, integer *n, f2c_complex *a, 
	integer *lda, integer *info);

/* Subroutine */ int ctrtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, f2c_complex *a, integer *lda, f2c_complex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int ctzrqf_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, integer *info);

/* Subroutine */ int ctzrzf_(integer *m, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cung2l_(integer *m, integer *n, integer *k, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cung2r_(integer *m, integer *n, integer *k, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cungbr_(char *vect, integer *m, integer *n, integer *k, 
	f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *work, integer *lwork,
	integer *info);

/* Subroutine */ int cunghr_(integer *n, integer *ilo, integer *ihi, f2c_complex *
	a, integer *lda, f2c_complex *tau, f2c_complex *work, integer *lwork, integer 
	*info);

/* Subroutine */ int cungl2_(integer *m, integer *n, integer *k, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cunglq_(integer *m, integer *n, integer *k, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cungql_(integer *m, integer *n, integer *k, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cungqr_(integer *m, integer *n, integer *k, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cungr2_(integer *m, integer *n, integer *k, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work, integer *info);

/* Subroutine */ int cungrq_(integer *m, integer *n, integer *k, f2c_complex *a, 
	integer *lda, f2c_complex *tau, f2c_complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cungtr_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	f2c_complex *tau, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunm2l_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *info);

/* Subroutine */ int cunm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *info);

/* Subroutine */ int cunmbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, 
	f2c_complex *c__, integer *ldc, f2c_complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cunmhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, f2c_complex *a, integer *lda, f2c_complex *tau, 
	f2c_complex *c__, integer *ldc, f2c_complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cunml2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *info);

/* Subroutine */ int cunmlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunmql_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunmqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunmr2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *info);

/* Subroutine */ int cunmr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, f2c_complex *a, integer *lda, f2c_complex *tau, 
	f2c_complex *c__, integer *ldc, f2c_complex *work, integer *info);

/* Subroutine */ int cunmrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunmrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, f2c_complex *a, integer *lda, f2c_complex *tau, 
	f2c_complex *c__, integer *ldc, f2c_complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cunmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, f2c_complex *a, integer *lda, f2c_complex *tau, f2c_complex *c__, 
	integer *ldc, f2c_complex *work, integer *lwork, integer *info);

/* Subroutine */ int cupgtr_(char *uplo, integer *n, f2c_complex *ap, f2c_complex *
	tau, f2c_complex *q, integer *ldq, f2c_complex *work, integer *info);

/* Subroutine */ int cupmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, f2c_complex *ap, f2c_complex *tau, f2c_complex *c__, integer *ldc, 
	f2c_complex *work, integer *info);

/* Subroutine */ int dbdsdc_(char *uplo, char *compq, integer *n, doublereal *
	d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, 
	integer *ldvt, doublereal *q, integer *iq, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ int dbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt, 
	integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *
	ldc, doublereal *work, integer *info);

/* Subroutine */ int ddisna_(char *job, integer *m, integer *n, doublereal *
	d__, doublereal *sep, integer *info);

/* Subroutine */ int dgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *
	d__, doublereal *e, doublereal *q, integer *ldq, doublereal *pt, 
	integer *ldpt, doublereal *c__, integer *ldc, doublereal *work, 
	integer *info);

/* Subroutine */ int dgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	doublereal *ab, integer *ldab, integer *ipiv, doublereal *anorm, 
	doublereal *rcond, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info);

/* Subroutine */ int dgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, 
	integer *ldafb, integer *ipiv, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dgbsv_(integer *n, integer *kl, integer *ku, integer *
	nrhs, doublereal *ab, integer *ldab, integer *ipiv, doublereal *b, 
	integer *ldb, integer *info);

/* Subroutine */ int dgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	integer *ku, integer *nrhs, doublereal *ab, integer *ldab, 
	doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	doublereal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int dgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	doublereal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int dgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv, 
	doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *
	ldv, integer *info);

/* Subroutine */ int dgebal_(char *job, integer *n, doublereal *a, integer *
	lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);

/* Subroutine */ int dgebd2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	taup, doublereal *work, integer *info);

/* Subroutine */ int dgebrd_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	taup, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dgecon_(char *norm, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ int dgeequ_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, integer *info);

/* Subroutine */ int dgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	doublereal *a, integer *lda, integer *sdim, doublereal *wr, 
	doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work, 
	integer *lwork, logical *bwork, integer *info);

/* Subroutine */ int dgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, doublereal *a, integer *lda, integer *sdim, 
	doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, 
	doublereal *rconde, doublereal *rcondv, doublereal *work, integer *
	lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, 
	integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, 
	integer *lwork, integer *info);

/* Subroutine */ int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublereal *a, integer *lda, doublereal *wr, 
	doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, 
	integer *ldvr, integer *ilo, integer *ihi, doublereal *scale, 
	doublereal *abnrm, doublereal *rconde, doublereal *rcondv, doublereal 
	*work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int dgegs_(char *jobvsl, char *jobvsr, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	alphar, doublereal *alphai, doublereal *beta, doublereal *vsl, 
	integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work, 
	integer *lwork, integer *info);

/* Subroutine */ int dgegv_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dgehd2_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *info);

/* Subroutine */ int dgehrd_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info);

/* Subroutine */ int dgelq2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dgelqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dgels_(char *trans, integer *m, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dgelsd_(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	integer *iwork, integer *info);

/* Subroutine */ int dgelss_(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ int dgelsx_(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	info);

/* Subroutine */ int dgelsy_(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	lwork, integer *info);

/* Subroutine */ int dgeql2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dgeqlf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dgeqp3_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ int dgeqpf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dgeqr2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dgeqrf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dgerfs_(char *trans, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	integer *info);

/* Subroutine */ int dgerq2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dgerqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dgesc2_(integer *n, doublereal *a, integer *lda, 
	doublereal *rhs, integer *ipiv, integer *jpiv, doublereal *scale);

/* Subroutine */ int dgesdd_(char *jobz, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *s, doublereal *u, integer *ldu, 
	doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	integer *iwork, integer *info);

/* Subroutine */ int dgesv_(integer *n, integer *nrhs, doublereal *a, integer 
	*lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ int dgetc2_(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, integer *jpiv, integer *info);

/* Subroutine */ int dgetf2_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);

/* Subroutine */ int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);

/* Subroutine */ int dgetri_(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dgetrs_(char *trans, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info);

/* Subroutine */ int dggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublereal *v, integer *ldv, integer *info);

/* Subroutine */ int dggbal_(char *job, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, integer *ilo, integer *ihi, 
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info);

/* Subroutine */ int dgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, integer *sdim, doublereal *alphar, doublereal *alphai, 
	doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr, 
	integer *ldvsr, doublereal *work, integer *lwork, logical *bwork, 
	integer *info);

/* Subroutine */ int dggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, char *sense, integer *n, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, integer *sdim, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl,
	doublereal *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	rcondv, doublereal *work, integer *lwork, integer *iwork, integer *
	liwork, logical *bwork, integer *info);

/* Subroutine */ int dggev_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, 
	doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	rcondv, doublereal *work, integer *lwork, integer *iwork, logical *
	bwork, integer *info);

/* Subroutine */ int dggglm_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *d__, 
	doublereal *x, doublereal *y, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *q, integer *ldq, doublereal *z__, integer *
	ldz, integer *info);

/* Subroutine */ int dgglse_(integer *m, integer *n, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	doublereal *d__, doublereal *x, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dggqrf_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *taua, doublereal *b, integer *ldb, 
	doublereal *taub, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dggrqf_(integer *m, integer *p, integer *n, doublereal *
	a, integer *lda, doublereal *taua, doublereal *b, integer *ldb, 
	doublereal *taub, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer 
	*ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork, 
	integer *info);

/* Subroutine */ int dggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer 
	*l, doublereal *u, integer *ldu, doublereal *v, integer *ldv, 
	doublereal *q, integer *ldq, integer *iwork, doublereal *tau, 
	doublereal *work, integer *info);

/* Subroutine */ int dgtcon_(char *norm, integer *n, doublereal *dl, 
	doublereal *d__, doublereal *du, doublereal *du2, integer *ipiv, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ int dgtrfs_(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf, 
	doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);

/* Subroutine */ int dgtsv_(integer *n, integer *nrhs, doublereal *dl, 
	doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer 
	*info);

/* Subroutine */ int dgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *
	dlf, doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ int dgttrf_(integer *n, doublereal *dl, doublereal *d__, 
	doublereal *du, doublereal *du2, integer *ipiv, integer *info);

/* Subroutine */ int dgttrs_(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, 
	integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dgtts2_(integer *itrans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, 
	integer *ipiv, doublereal *b, integer *ldb);

/* Subroutine */ int dhgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dhsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, doublereal *h__, integer *ldh, doublereal *wr, 
	doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, 
	integer *ldvr, integer *mm, integer *m, doublereal *work, integer *
	ifaill, integer *ifailr, integer *info);

/* Subroutine */ int dhseqr_(char *job, char *compz, integer *n, integer *ilo,
	integer *ihi, doublereal *h__, integer *ldh, doublereal *wr, 
	doublereal *wi, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *info);

/* Subroutine */ int dlabad_(doublereal *small, doublereal *large);

/* Subroutine */ int dlabrd_(integer *m, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq, 
	doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer 
	*ldy);

/* Subroutine */ int dlacon_(integer *n, doublereal *v, doublereal *x, 
	integer *isgn, doublereal *est, integer *kase);

/* Subroutine */ int dlacpy_(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb);

/* Subroutine */ int dladiv_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *p, doublereal *q);

/* Subroutine */ int dlae2_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *rt1, doublereal *rt2);

/* Subroutine */ int dlaebz_(integer *ijob, integer *nitmax, integer *n, 
	integer *mmax, integer *minp, integer *nbmin, doublereal *abstol, 
	doublereal *reltol, doublereal *pivmin, doublereal *d__, doublereal *
	e, doublereal *e2, integer *nval, doublereal *ab, doublereal *c__, 
	integer *mout, integer *nab, doublereal *work, integer *iwork, 
	integer *info);

/* Subroutine */ int dlaed0_(integer *icompq, integer *qsiz, integer *n, 
	doublereal *d__, doublereal *e, doublereal *q, integer *ldq, 
	doublereal *qstore, integer *ldqs, doublereal *work, integer *iwork, 
	integer *info);

/* Subroutine */ int dlaed1_(integer *n, doublereal *d__, doublereal *q, 
	integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt, 
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dlaed2_(integer *k, integer *n, integer *n1, doublereal *
	d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, 
	doublereal *z__, doublereal *dlamda, doublereal *w, doublereal *q2, 
	integer *indx, integer *indxc, integer *indxp, integer *coltyp, 
	integer *info);

/* Subroutine */ int dlaed3_(integer *k, integer *n, integer *n1, doublereal *
	d__, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda,
	doublereal *q2, integer *indx, integer *ctot, doublereal *w, 
	doublereal *s, integer *info);

/* Subroutine */ int dlaed4_(integer *n, integer *i__, doublereal *d__, 
	doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dlam,
	integer *info);

/* Subroutine */ int dlaed5_(integer *i__, doublereal *d__, doublereal *z__, 
	doublereal *delta, doublereal *rho, doublereal *dlam);

/* Subroutine */ int dlaed6_(integer *kniter, logical *orgati, doublereal *
	rho, doublereal *d__, doublereal *z__, doublereal *finit, doublereal *
	tau, integer *info);

/* Subroutine */ int dlaed7_(integer *icompq, integer *n, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__, 
	doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer 
	*cutpnt, doublereal *qstore, integer *qptr, integer *prmptr, integer *
	perm, integer *givptr, integer *givcol, doublereal *givnum, 
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dlaed8_(integer *icompq, integer *k, integer *n, integer 
	*qsiz, doublereal *d__, doublereal *q, integer *ldq, integer *indxq, 
	doublereal *rho, integer *cutpnt, doublereal *z__, doublereal *dlamda,
	doublereal *q2, integer *ldq2, doublereal *w, integer *perm, integer 
	*givptr, integer *givcol, doublereal *givnum, integer *indxp, integer 
	*indx, integer *info);

/* Subroutine */ int dlaed9_(integer *k, integer *kstart, integer *kstop, 
	integer *n, doublereal *d__, doublereal *q, integer *ldq, doublereal *
	rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds, 
	integer *info);

/* Subroutine */ int dlaeda_(integer *n, integer *tlvls, integer *curlvl, 
	integer *curpbm, integer *prmptr, integer *perm, integer *givptr, 
	integer *givcol, doublereal *givnum, doublereal *q, integer *qptr, 
	doublereal *z__, doublereal *ztemp, integer *info);

/* Subroutine */ int dlaein_(logical *rightv, logical *noinit, integer *n, 
	doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi, 
	doublereal *vr, doublereal *vi, doublereal *b, integer *ldb, 
	doublereal *work, doublereal *eps3, doublereal *smlnum, doublereal *
	bignum, integer *info);

/* Subroutine */ int dlaev2_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1);

/* Subroutine */ int dlaexc_(logical *wantq, integer *n, doublereal *t, 
	integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1, 
	integer *n2, doublereal *work, integer *info);

/* Subroutine */ int dlag2_(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *safmin, doublereal *scale1, doublereal *
	scale2, doublereal *wr1, doublereal *wr2, doublereal *wi);

/* Subroutine */ int dlags2_(logical *upper, doublereal *a1, doublereal *a2, 
	doublereal *a3, doublereal *b1, doublereal *b2, doublereal *b3, 
	doublereal *csu, doublereal *snu, doublereal *csv, doublereal *snv, 
	doublereal *csq, doublereal *snq);

/* Subroutine */ int dlagtf_(integer *n, doublereal *a, doublereal *lambda, 
	doublereal *b, doublereal *c__, doublereal *tol, doublereal *d__, 
	integer *in, integer *info);

/* Subroutine */ int dlagtm_(char *trans, integer *n, integer *nrhs, 
	doublereal *alpha, doublereal *dl, doublereal *d__, doublereal *du, 
	doublereal *x, integer *ldx, doublereal *beta, doublereal *b, integer 
	*ldb);

/* Subroutine */ int dlagts_(integer *job, integer *n, doublereal *a, 
	doublereal *b, doublereal *c__, doublereal *d__, integer *in, 
	doublereal *y, doublereal *tol, integer *info);

/* Subroutine */ int dlagv2_(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *
	snr);

/* Subroutine */ int dlahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal 
	*wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__, 
	integer *ldz, integer *info);

/* Subroutine */ int dlahrd_(integer *n, integer *k, integer *nb, doublereal *
	a, integer *lda, doublereal *tau, doublereal *t, integer *ldt, 
	doublereal *y, integer *ldy);

/* Subroutine */ int dlaic1_(integer *job, integer *j, doublereal *x, 
	doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
	sestpr, doublereal *s, doublereal *c__);

/* Subroutine */ int dlaln2_(logical *ltrans, integer *na, integer *nw, 
	doublereal *smin, doublereal *ca, doublereal *a, integer *lda, 
	doublereal *d1, doublereal *d2, doublereal *b, integer *ldb, 
	doublereal *wr, doublereal *wi, doublereal *x, integer *ldx, 
	doublereal *scale, doublereal *xnorm, integer *info);

/* Subroutine */ int dlals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, doublereal *b, integer *ldb, doublereal 
	*bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *
	poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *
	k, doublereal *c__, doublereal *s, doublereal *work, integer *info);

/* Subroutine */ int dlalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, doublereal *b, integer *ldb, doublereal *bx, integer *
	ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *k, 
	doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	work, integer *iwork, integer *info);

/* Subroutine */ int dlalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, doublereal *d__, doublereal *e, doublereal *b, integer *ldb, 
	doublereal *rcond, integer *rank, doublereal *work, integer *iwork, 
	integer *info);

/* Subroutine */ int dlamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1);

/* Subroutine */ int dlamc2_(integer *beta, integer *t, logical *rnd, 
	doublereal *eps, integer *emin, doublereal *rmin, integer *emax, 
	doublereal *rmax);

/* Subroutine */ int dlamc4_(integer *emin, doublereal *start, integer *base);

/* Subroutine */ int dlamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, doublereal *rmax);

/* Subroutine */ int dlamrg_(integer *n1, integer *n2, doublereal *a, integer 
	*dtrd1, integer *dtrd2, integer *index);

/* Subroutine */ int dlanv2_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r,
	doublereal *rt2i, doublereal *cs, doublereal *sn);

/* Subroutine */ int dlapll_(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *ssmin);

/* Subroutine */ int dlapmt_(logical *forwrd, integer *m, integer *n, 
	doublereal *x, integer *ldx, integer *k);

/* Subroutine */ int dlaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int dlaqge_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, char *equed);

/* Subroutine */ int dlaqp2_(integer *m, integer *n, integer *offset, 
	doublereal *a, integer *lda, integer *jpvt, doublereal *tau, 
	doublereal *vn1, doublereal *vn2, doublereal *work);

/* Subroutine */ int dlaqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, doublereal *a, integer *lda, integer *jpvt, 
	doublereal *tau, doublereal *vn1, doublereal *vn2, doublereal *auxv, 
	doublereal *f, integer *ldf);

/* Subroutine */ int dlaqsb_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	char *equed);

/* Subroutine */ int dlaqsp_(char *uplo, integer *n, doublereal *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int dlaqsy_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int dlaqtr_(logical *ltran, logical *lf2c_real, integer *n, 
	doublereal *t, integer *ldt, doublereal *b, doublereal *w, doublereal 
	*scale, doublereal *x, doublereal *work, integer *info);

/* Subroutine */ int dlar1v_(integer *n, integer *b1, integer *bn, doublereal 
	*sigma, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	lld, doublereal *gersch, doublereal *z__, doublereal *ztz, doublereal 
	*mingma, integer *r__, integer *isuppz, doublereal *work);

/* Subroutine */ int dlar2v_(integer *n, doublereal *x, doublereal *y, 
	doublereal *z__, integer *incx, doublereal *c__, doublereal *s, 
	integer *incc);

/* Subroutine */ int dlarf_(char *side, integer *m, integer *n, doublereal *v,
	integer *incv, doublereal *tau, doublereal *c__, integer *ldc, 
	doublereal *work);

/* Subroutine */ int dlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, 
	doublereal *work, integer *ldwork);

/* Subroutine */ int dlarfg_(integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *tau);

/* Subroutine */ int dlarft_(char *direct, char *storev, integer *n, integer *
	k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, 
	integer *ldt);

/* Subroutine */ int dlarfx_(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work);

/* Subroutine */ int dlargv_(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *c__, integer *incc);

/* Subroutine */ int dlarnv_(integer *idist, integer *iseed, integer *n, 
	doublereal *x);

/* Subroutine */ int dlarrb_(integer *n, doublereal *d__, doublereal *l, 
	doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast, 
	doublereal *sigma, doublereal *reltol, doublereal *w, doublereal *
	wgap, doublereal *werr, doublereal *work, integer *iwork, integer *
	info);

/* Subroutine */ int dlarre_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *tol, integer *nsplit, integer *isplit, integer *m, 
	doublereal *w, doublereal *woff, doublereal *gersch, doublereal *work,
	integer *info);

/* Subroutine */ int dlarrf_(integer *n, doublereal *d__, doublereal *l, 
	doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast, 
	doublereal *w, doublereal *dplus, doublereal *lplus, doublereal *work,
	integer *iwork, integer *info);

/* Subroutine */ int dlarrv_(integer *n, doublereal *d__, doublereal *l, 
	integer *isplit, integer *m, doublereal *w, integer *iblock, 
	doublereal *gersch, doublereal *tol, doublereal *z__, integer *ldz, 
	integer *isuppz, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dlartg_(doublereal *f, doublereal *g, doublereal *cs, 
	doublereal *sn, doublereal *r__);

/* Subroutine */ int dlartv_(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *c__, doublereal *s, integer 
	*incc);

/* Subroutine */ int dlaruv_(integer *iseed, integer *n, doublereal *x);

/* Subroutine */ int dlarz_(char *side, integer *m, integer *n, integer *l, 
	doublereal *v, integer *incv, doublereal *tau, doublereal *c__, 
	integer *ldc, doublereal *work);

/* Subroutine */ int dlarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublereal *v,
	integer *ldv, doublereal *t, integer *ldt, doublereal *c__, integer *
	ldc, doublereal *work, integer *ldwork);

/* Subroutine */ int dlarzt_(char *direct, char *storev, integer *n, integer *
	k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, 
	integer *ldt);

/* Subroutine */ int dlas2_(doublereal *f, doublereal *g, doublereal *h__, 
	doublereal *ssmin, doublereal *ssmax);

/* Subroutine */ int dlascl_(char *type__, integer *kl, integer *ku, 
	doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	doublereal *a, integer *lda, integer *info);

/* Subroutine */ int dlasd0_(integer *n, integer *sqre, doublereal *d__, 
	doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *
	ldvt, integer *smlsiz, integer *iwork, doublereal *work, integer *
	info);

/* Subroutine */ int dlasd1_(integer *nl, integer *nr, integer *sqre, 
	doublereal *d__, doublereal *alpha, doublereal *beta, doublereal *u, 
	integer *ldu, doublereal *vt, integer *ldvt, integer *idxq, integer *
	iwork, doublereal *work, integer *info);

/* Subroutine */ int dlasd2_(integer *nl, integer *nr, integer *sqre, integer 
	*k, doublereal *d__, doublereal *z__, doublereal *alpha, doublereal *
	beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, 
	doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2, 
	integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
	idxq, integer *coltyp, integer *info);

/* Subroutine */ int dlasd3_(integer *nl, integer *nr, integer *sqre, integer 
	*k, doublereal *d__, doublereal *q, integer *ldq, doublereal *dsigma, 
	doublereal *u, integer *ldu, doublereal *u2, integer *ldu2, 
	doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2, 
	integer *idxc, integer *ctot, doublereal *z__, integer *info);

/* Subroutine */ int dlasd4_(integer *n, integer *i__, doublereal *d__, 
	doublereal *z__, doublereal *delta, doublereal *rho, doublereal *
	sigma, doublereal *work, integer *info);

/* Subroutine */ int dlasd5_(integer *i__, doublereal *d__, doublereal *z__, 
	doublereal *delta, doublereal *rho, doublereal *dsigma, doublereal *
	work);

/* Subroutine */ int dlasd6_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, doublereal *d__, doublereal *vf, doublereal *vl, 
	doublereal *alpha, doublereal *beta, integer *idxq, integer *perm, 
	integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *
	difr, doublereal *z__, integer *k, doublereal *c__, doublereal *s, 
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dlasd7_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *k, doublereal *d__, doublereal *z__, 
	doublereal *zw, doublereal *vf, doublereal *vfw, doublereal *vl, 
	doublereal *vlw, doublereal *alpha, doublereal *beta, doublereal *
	dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, 
	integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	integer *ldgnum, doublereal *c__, doublereal *s, integer *info);

/* Subroutine */ int dlasd8_(integer *icompq, integer *k, doublereal *d__, 
	doublereal *z__, doublereal *vf, doublereal *vl, doublereal *difl, 
	doublereal *difr, integer *lddifr, doublereal *dsigma, doublereal *
	work, integer *info);

/* Subroutine */ int dlasd9_(integer *icompq, integer *ldu, integer *k, 
	doublereal *d__, doublereal *z__, doublereal *vf, doublereal *vl, 
	doublereal *difl, doublereal *difr, doublereal *dsigma, doublereal *
	work, integer *info);

/* Subroutine */ int dlasda_(integer *icompq, integer *smlsiz, integer *n, 
	integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer 
	*ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, 
	doublereal *z__, doublereal *poles, integer *givptr, integer *givcol, 
	integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__, 
	doublereal *s, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dlasdq_(char *uplo, integer *sqre, integer *n, integer *
	ncvt, integer *nru, integer *ncc, doublereal *d__, doublereal *e, 
	doublereal *vt, integer *ldvt, doublereal *u, integer *ldu, 
	doublereal *c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int dlasdt_(integer *n, integer *lvl, integer *nd, integer *
	inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ int dlaset_(char *uplo, integer *m, integer *n, doublereal *
	alpha, doublereal *beta, doublereal *a, integer *lda);

/* Subroutine */ int dlasq1_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *work, integer *info);

/* Subroutine */ int dlasq2_(integer *n, doublereal *z__, integer *info);

/* Subroutine */ int dlasq3_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *dmin__, doublereal *sigma, doublereal *desig,
	doublereal *qmax, integer *nfail, integer *iter, integer *ndiv, 
	logical *ieee);

/* Subroutine */ int dlasq4_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, integer *n0in, doublereal *dmin__, doublereal *dmin1, 
	doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, 
	doublereal *tau, integer *ttype);

/* Subroutine */ int dlasq5_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *tau, doublereal *dmin__, doublereal *dmin1, 
	doublereal *dmin2, doublereal *dn, doublereal *dnm1, doublereal *dnm2,
	logical *ieee);

/* Subroutine */ int dlasq6_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2,
	doublereal *dn, doublereal *dnm1, doublereal *dnm2);

/* Subroutine */ int dlasr_(char *side, char *pivot, char *direct, integer *m,
	integer *n, doublereal *c__, doublereal *s, doublereal *a, integer *
	lda);

/* Subroutine */ int dlasrt_(char *id, integer *n, doublereal *d__, integer *
	info);

/* Subroutine */ int dlassq_(integer *n, doublereal *x, integer *incx, 
	doublereal *scale, doublereal *sumsq);

/* Subroutine */ int dlasv2_(doublereal *f, doublereal *g, doublereal *h__, 
	doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *
	csr, doublereal *snl, doublereal *csl);

/* Subroutine */ int dlaswp_(integer *n, doublereal *a, integer *lda, integer 
	*k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int dlasy2_(logical *ltranl, logical *ltranr, integer *isgn, 
	integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
	tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale, 
	doublereal *x, integer *ldx, doublereal *xnorm, integer *info);

/* Subroutine */ int dlasyf_(char *uplo, integer *n, integer *nb, integer *kb,
	doublereal *a, integer *lda, integer *ipiv, doublereal *w, integer *
	ldw, integer *info);

/* Subroutine */ int dlatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublereal *ab, integer *ldab, 
	doublereal *x, doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int dlatdf_(integer *ijob, integer *n, doublereal *z__, 
	integer *ldz, doublereal *rhs, doublereal *rdsum, doublereal *rdscal, 
	integer *ipiv, integer *jpiv);

/* Subroutine */ int dlatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *ap, doublereal *x, doublereal *scale, 
	doublereal *cnorm, integer *info);

/* Subroutine */ int dlatrd_(char *uplo, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *e, doublereal *tau, doublereal *w, 
	integer *ldw);

/* Subroutine */ int dlatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *a, integer *lda, doublereal *x, 
	doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int dlatrz_(integer *m, integer *n, integer *l, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work);

/* Subroutine */ int dlatzm_(char *side, integer *m, integer *n, doublereal *
	v, integer *incv, doublereal *tau, doublereal *c1, doublereal *c2, 
	integer *ldc, doublereal *work);

/* Subroutine */ int dlauu2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ int dlauum_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ int dopgtr_(char *uplo, integer *n, doublereal *ap, 
	doublereal *tau, doublereal *q, integer *ldq, doublereal *work, 
	integer *info);

/* Subroutine */ int dopmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublereal *ap, doublereal *tau, doublereal *c__, integer 
	*ldc, doublereal *work, integer *info);

/* Subroutine */ int dorg2l_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dorg2r_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dorgbr_(char *vect, integer *m, integer *n, integer *k, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info);

/* Subroutine */ int dorghr_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info);

/* Subroutine */ int dorgl2_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dorglq_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dorgql_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dorgqr_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dorgr2_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int dorgrq_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dorgtr_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dorm2l_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int dorm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int dormbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dormhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dorml2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int dormlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dormql_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dormqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dormr2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int dormr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int dormrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dormrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dormtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dpbcon_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublereal *
	work, integer *iwork, integer *info);

/* Subroutine */ int dpbequ_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	integer *info);

/* Subroutine */ int dpbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);

/* Subroutine */ int dpbstf_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);

/* Subroutine */ int dpbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb, 
	integer *info);

/* Subroutine */ int dpbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, 
	integer *ldafb, char *equed, doublereal *s, doublereal *b, integer *
	ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr,
	doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dpbtf2_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);

/* Subroutine */ int dpbtrf_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);

/* Subroutine */ int dpbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb, 
	integer *info);

/* Subroutine */ int dpocon_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ int dpoequ_(integer *n, doublereal *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int dporfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);

/* Subroutine */ int dposv_(char *uplo, integer *n, integer *nrhs, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *
	x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *
	berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dpotf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ int dpotrf_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ int dpotri_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ int dpotrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ int dppcon_(char *uplo, integer *n, doublereal *ap, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ int dppequ_(char *uplo, integer *n, doublereal *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int dpprfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, doublereal *afp, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dppsv_(char *uplo, integer *n, integer *nrhs, doublereal 
	*ap, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *ap, doublereal *afp, char *equed, doublereal *s, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ int dpptrf_(char *uplo, integer *n, doublereal *ap, integer *
	info);

/* Subroutine */ int dpptri_(char *uplo, integer *n, doublereal *ap, integer *
	info);

/* Subroutine */ int dpptrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dptcon_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *info);

/* Subroutine */ int dpteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int dptrfs_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer 
	*ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	doublereal *work, integer *info);

/* Subroutine */ int dptsv_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dptsvx_(char *fact, integer *n, integer *nrhs, 
	doublereal *d__, doublereal *e, doublereal *df, doublereal *ef, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	info);

/* Subroutine */ int dpttrf_(integer *n, doublereal *d__, doublereal *e, 
	integer *info);

/* Subroutine */ int dpttrs_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dptts2_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb);

/* Subroutine */ int drscl_(integer *n, doublereal *sa, doublereal *sx, 
	integer *incx);

/* Subroutine */ int dsbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *info);

/* Subroutine */ int dsbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int dsbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, doublereal *ab, integer *ldab, doublereal *q, integer *
	ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu, 
	doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *iwork, integer *ifail, 
	integer *info);

/* Subroutine */ int dsbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *x, integer *ldx, doublereal *work, integer *info);

/* Subroutine */ int dsbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int dsbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int dsbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, doublereal *ab, integer *ldab, doublereal *
	bb, integer *ldbb, doublereal *q, integer *ldq, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	*m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int dsbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *d__, doublereal *e, 
	doublereal *q, integer *ldq, doublereal *work, integer *info);

/* Subroutine */ int dspcon_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer 
	*iwork, integer *info);

/* Subroutine */ int dspev_(char *jobz, char *uplo, integer *n, doublereal *
	ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int dspevd_(char *jobz, char *uplo, integer *n, doublereal *
	ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int dspevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *ap, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *iwork, integer *ifail, 
	integer *info);

/* Subroutine */ int dspgst_(integer *itype, char *uplo, integer *n, 
	doublereal *ap, doublereal *bp, integer *info);

/* Subroutine */ int dspgv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *info);

/* Subroutine */ int dspgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int dspgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublereal *ap, doublereal *bp, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	*m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int dsprfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b, 
	integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dspsv_(char *uplo, integer *n, integer *nrhs, doublereal 
	*ap, integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int dspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b, 
	integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, 
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	integer *info);

/* Subroutine */ int dsptrd_(char *uplo, integer *n, doublereal *ap, 
	doublereal *d__, doublereal *e, doublereal *tau, integer *info);

/* Subroutine */ int dsptrf_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, integer *info);

/* Subroutine */ int dsptri_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *work, integer *info);

/* Subroutine */ int dsptrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, integer *ipiv, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ int dstebz_(char *range, char *order, integer *n, doublereal 
	*vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
	doublereal *d__, doublereal *e, integer *m, integer *nsplit, 
	doublereal *w, integer *iblock, integer *isplit, doublereal *work, 
	integer *iwork, integer *info);

/* Subroutine */ int dstedc_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int dstegr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int dstein_(integer *n, doublereal *d__, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info);

/* Subroutine */ int dsteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int dsterf_(integer *n, doublereal *d__, doublereal *e, 
	integer *info);

/* Subroutine */ int dstev_(char *jobz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int dstevd_(char *jobz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int dstevr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int dstevx_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info);

/* Subroutine */ int dsycon_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *
	work, integer *iwork, integer *info);

/* Subroutine */ int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a,
	integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *info);

/* Subroutine */ int dsyevd_(char *jobz, char *uplo, integer *n, doublereal *
	a, integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int dsyevr_(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
	integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int dsygs2_(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ int dsygst_(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ int dsygv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *w, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dsygvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *w, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int dsygvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublereal *a, integer *lda, doublereal *b, integer 
	*ldb, doublereal *vl, doublereal *vu, integer *il, integer *iu, 
	doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *ifail, integer *info);

/* Subroutine */ int dsyrfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	integer *info);

/* Subroutine */ int dsysv_(char *uplo, integer *n, integer *nrhs, doublereal 
	*a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, 
	doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dsysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *
	ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int dsytd2_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info);

/* Subroutine */ int dsytf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);

/* Subroutine */ int dsytrd_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, doublereal *
	work, integer *lwork, integer *info);

/* Subroutine */ int dsytrf_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dsytri_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *work, integer *info);

/* Subroutine */ int dsytrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info);

/* Subroutine */ int dtbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, doublereal *ab, integer *ldab, doublereal *rcond, 
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dtbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal 
	*b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dtbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal 
	*b, integer *ldb, integer *info);

/* Subroutine */ int dtgevc_(char *side, char *howmny, logical *select, 
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer 
	*mm, integer *m, doublereal *work, integer *info);

/* Subroutine */ int dtgex2_(logical *wantq, logical *wantz, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *j1, integer *
	n1, integer *n2, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dtgexc_(logical *wantq, logical *wantz, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *ifst, 
	integer *ilst, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int dtgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	integer *m, doublereal *pl, doublereal *pr, doublereal *dif, 
	doublereal *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int dtgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *tola, 
	doublereal *tolb, doublereal *alpha, doublereal *beta, doublereal *u, 
	integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *
	ldq, doublereal *work, integer *ncycle, integer *info);

/* Subroutine */ int dtgsna_(char *job, char *howmny, logical *select, 
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, 
	doublereal *s, doublereal *dif, integer *mm, integer *m, doublereal *
	work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int dtgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *rdsum, doublereal *rdscal, integer *iwork, integer 
	*pq, integer *info);

/* Subroutine */ int dtgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *dif, doublereal *work, integer *lwork, integer *
	iwork, integer *info);

/* Subroutine */ int dtpcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublereal *ap, doublereal *rcond, doublereal *work, integer *iwork, 
	integer *info);

/* Subroutine */ int dtprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dtptri_(char *uplo, char *diag, integer *n, doublereal *
	ap, integer *info);

/* Subroutine */ int dtptrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ int dtrcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublereal *a, integer *lda, doublereal *rcond, doublereal *work, 
	integer *iwork, integer *info);

/* Subroutine */ int dtrevc_(char *side, char *howmny, logical *select, 
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, 
	doublereal *work, integer *info);

/* Subroutine */ int dtrexc_(char *compq, integer *n, doublereal *t, integer *
	ldt, doublereal *q, integer *ldq, integer *ifst, integer *ilst, 
	doublereal *work, integer *info);

/* Subroutine */ int dtrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int dtrsen_(char *job, char *compq, logical *select, integer 
	*n, doublereal *t, integer *ldt, doublereal *q, integer *ldq, 
	doublereal *wr, doublereal *wi, integer *m, doublereal *s, doublereal 
	*sep, doublereal *work, integer *lwork, integer *iwork, integer *
	liwork, integer *info);

/* Subroutine */ int dtrsna_(char *job, char *howmny, logical *select, 
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, doublereal *s, doublereal *sep, 
	integer *mm, integer *m, doublereal *work, integer *ldwork, integer *
	iwork, integer *info);

/* Subroutine */ int dtrsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *scale, integer *info);

/* Subroutine */ int dtrti2_(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info);

/* Subroutine */ int dtrtri_(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info);

/* Subroutine */ int dtrtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, integer *info);

/* Subroutine */ int dtzrqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, integer *info);

/* Subroutine */ int dtzrzf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

integer icmax1_(integer *n, f2c_complex *cx, integer *incx);

integer ieeeck_(integer *ispec, f2c_real *zero, f2c_real *one);

integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen 
	opts_len);

integer izmax1_(integer *n, doublecomplex *cx, integer *incx);

/* Subroutine */ int sbdsdc_(char *uplo, char *compq, integer *n, f2c_real *d__, 
	f2c_real *e, f2c_real *u, integer *ldu, f2c_real *vt, integer *ldvt, f2c_real *q, 
	integer *iq, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, f2c_real *d__, f2c_real *e, f2c_real *vt, integer *ldvt, f2c_real *
	u, integer *ldu, f2c_real *c__, integer *ldc, f2c_real *work, integer *info);

/* Subroutine */ int sdisna_(char *job, integer *m, integer *n, f2c_real *d__, 
	f2c_real *sep, integer *info);

/* Subroutine */ int sgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	integer *kl, integer *ku, f2c_real *ab, integer *ldab, f2c_real *d__, f2c_real *
	e, f2c_real *q, integer *ldq, f2c_real *pt, integer *ldpt, f2c_real *c__, integer 
	*ldc, f2c_real *work, integer *info);

/* Subroutine */ int sgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	f2c_real *ab, integer *ldab, integer *ipiv, f2c_real *anorm, f2c_real *rcond, 
	f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	f2c_real *ab, integer *ldab, f2c_real *r__, f2c_real *c__, f2c_real *rowcnd, f2c_real *
	colcnd, f2c_real *amax, integer *info);

/* Subroutine */ int sgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, f2c_real *ab, integer *ldab, f2c_real *afb, integer *ldafb,
	integer *ipiv, f2c_real *b, integer *ldb, f2c_real *x, integer *ldx, f2c_real *
	ferr, f2c_real *berr, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sgbsv_(integer *n, integer *kl, integer *ku, integer *
	nrhs, f2c_real *ab, integer *ldab, integer *ipiv, f2c_real *b, integer *ldb, 
	integer *info);

/* Subroutine */ int sgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	integer *ku, integer *nrhs, f2c_real *ab, integer *ldab, f2c_real *afb, 
	integer *ldafb, integer *ipiv, char *equed, f2c_real *r__, f2c_real *c__, 
	f2c_real *b, integer *ldb, f2c_real *x, integer *ldx, f2c_real *rcond, f2c_real *ferr,
	f2c_real *berr, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	f2c_real *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int sgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	f2c_real *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int sgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, f2c_real *ab, integer *ldab, integer *ipiv, f2c_real *b, 
	integer *ldb, integer *info);

/* Subroutine */ int sgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, f2c_real *scale, integer *m, f2c_real *v, integer *ldv, integer 
	*info);

/* Subroutine */ int sgebal_(char *job, integer *n, f2c_real *a, integer *lda, 
	integer *ilo, integer *ihi, f2c_real *scale, integer *info);

/* Subroutine */ int sgebd2_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *d__, f2c_real *e, f2c_real *tauq, f2c_real *taup, f2c_real *work, integer *info);

/* Subroutine */ int sgebrd_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *d__, f2c_real *e, f2c_real *tauq, f2c_real *taup, f2c_real *work, integer *
	lwork, integer *info);

/* Subroutine */ int sgecon_(char *norm, integer *n, f2c_real *a, integer *lda, 
	f2c_real *anorm, f2c_real *rcond, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sgeequ_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *r__, f2c_real *c__, f2c_real *rowcnd, f2c_real *colcnd, f2c_real *amax, integer 
	*info);

/* Subroutine */ int sgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	f2c_real *a, integer *lda, integer *sdim, f2c_real *wr, f2c_real *wi, f2c_real *vs, 
	integer *ldvs, f2c_real *work, integer *lwork, logical *bwork, integer *
	info);

/* Subroutine */ int sgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, f2c_real *a, integer *lda, integer *sdim, f2c_real *wr, 
	f2c_real *wi, f2c_real *vs, integer *ldvs, f2c_real *rconde, f2c_real *rcondv, f2c_real *
	work, integer *lwork, integer *iwork, integer *liwork, logical *bwork,
	integer *info);

/* Subroutine */ int sgeev_(char *jobvl, char *jobvr, integer *n, f2c_real *a, 
	integer *lda, f2c_real *wr, f2c_real *wi, f2c_real *vl, integer *ldvl, f2c_real *vr, 
	integer *ldvr, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, f2c_real *a, integer *lda, f2c_real *wr, f2c_real *wi, f2c_real *
	vl, integer *ldvl, f2c_real *vr, integer *ldvr, integer *ilo, integer *
	ihi, f2c_real *scale, f2c_real *abnrm, f2c_real *rconde, f2c_real *rcondv, f2c_real *work,
	integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int sgegs_(char *jobvsl, char *jobvsr, integer *n, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, f2c_real *alphar, f2c_real *alphai, f2c_real 
	*beta, f2c_real *vsl, integer *ldvsl, f2c_real *vsr, integer *ldvsr, f2c_real *
	work, integer *lwork, integer *info);

/* Subroutine */ int sgegv_(char *jobvl, char *jobvr, integer *n, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, f2c_real *alphar, f2c_real *alphai, f2c_real 
	*beta, f2c_real *vl, integer *ldvl, f2c_real *vr, integer *ldvr, f2c_real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int sgehd2_(integer *n, integer *ilo, integer *ihi, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sgehrd_(integer *n, integer *ilo, integer *ihi, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgelq2_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sgelqf_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgels_(char *trans, integer *m, integer *n, integer *
	nrhs, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int sgelsd_(integer *m, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, f2c_real *s, f2c_real *rcond, integer *
	rank, f2c_real *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int sgelss_(integer *m, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, f2c_real *s, f2c_real *rcond, integer *
	rank, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgelsx_(integer *m, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, integer *jpvt, f2c_real *rcond, 
	integer *rank, f2c_real *work, integer *info);

/* Subroutine */ int sgelsy_(integer *m, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, integer *jpvt, f2c_real *rcond, 
	integer *rank, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgeql2_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sgeqlf_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgeqp3_(integer *m, integer *n, f2c_real *a, integer *lda, 
	integer *jpvt, f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgeqpf_(integer *m, integer *n, f2c_real *a, integer *lda, 
	integer *jpvt, f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sgeqr2_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sgeqrf_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgerfs_(char *trans, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *af, integer *ldaf, integer *ipiv, f2c_real *b, 
	integer *ldb, f2c_real *x, integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_real *
	work, integer *iwork, integer *info);

/* Subroutine */ int sgerq2_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sgerqf_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgesc2_(integer *n, f2c_real *a, integer *lda, f2c_real *rhs, 
	integer *ipiv, integer *jpiv, f2c_real *scale);

/* Subroutine */ int sgesdd_(char *jobz, integer *m, integer *n, f2c_real *a, 
	integer *lda, f2c_real *s, f2c_real *u, integer *ldu, f2c_real *vt, integer *ldvt,
	f2c_real *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int sgesv_(integer *n, integer *nrhs, f2c_real *a, integer *lda, 
	integer *ipiv, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	f2c_real *a, integer *lda, f2c_real *s, f2c_real *u, integer *ldu, f2c_real *vt, 
	integer *ldvt, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, f2c_real *a, integer *lda, f2c_real *af, integer *ldaf, integer *ipiv, 
	char *equed, f2c_real *r__, f2c_real *c__, f2c_real *b, integer *ldb, f2c_real *x, 
	integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int sgetc2_(integer *n, f2c_real *a, integer *lda, integer *ipiv,
	integer *jpiv, integer *info);

/* Subroutine */ int sgetf2_(integer *m, integer *n, f2c_real *a, integer *lda, 
	integer *ipiv, integer *info);

/* Subroutine */ int sgetrf_(integer *m, integer *n, f2c_real *a, integer *lda, 
	integer *ipiv, integer *info);

/* Subroutine */ int sgetri_(integer *n, f2c_real *a, integer *lda, integer *ipiv,
	f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgetrs_(char *trans, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, integer *ipiv, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, f2c_real *lscale, f2c_real *rscale, integer *m, f2c_real *v, 
	integer *ldv, integer *info);

/* Subroutine */ int sggbal_(char *job, integer *n, f2c_real *a, integer *lda, 
	f2c_real *b, integer *ldb, integer *ilo, integer *ihi, f2c_real *lscale, f2c_real 
	*rscale, f2c_real *work, integer *info);

/* Subroutine */ int sgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, integer *n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, 
	integer *sdim, f2c_real *alphar, f2c_real *alphai, f2c_real *beta, f2c_real *vsl, 
	integer *ldvsl, f2c_real *vsr, integer *ldvsr, f2c_real *work, integer *lwork,
	logical *bwork, integer *info);

/* Subroutine */ int sggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, integer *n, f2c_real *a, integer *lda, f2c_real *b, 
	integer *ldb, integer *sdim, f2c_real *alphar, f2c_real *alphai, f2c_real *beta, 
	f2c_real *vsl, integer *ldvsl, f2c_real *vsr, integer *ldvsr, f2c_real *rconde, 
	f2c_real *rcondv, f2c_real *work, integer *lwork, integer *iwork, integer *
	liwork, logical *bwork, integer *info);

/* Subroutine */ int sggev_(char *jobvl, char *jobvr, integer *n, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, f2c_real *alphar, f2c_real *alphai, f2c_real 
	*beta, f2c_real *vl, integer *ldvl, f2c_real *vr, integer *ldvr, f2c_real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int sggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real 
	*alphar, f2c_real *alphai, f2c_real *beta, f2c_real *vl, integer *ldvl, f2c_real *vr, 
	integer *ldvr, integer *ilo, integer *ihi, f2c_real *lscale, f2c_real *rscale,
	f2c_real *abnrm, f2c_real *bbnrm, f2c_real *rconde, f2c_real *rcondv, f2c_real *work, 
	integer *lwork, integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int sggglm_(integer *n, integer *m, integer *p, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, f2c_real *d__, f2c_real *x, f2c_real *y, 
	f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real 
	*q, integer *ldq, f2c_real *z__, integer *ldz, integer *info);

/* Subroutine */ int sgglse_(integer *m, integer *n, integer *p, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, f2c_real *c__, f2c_real *d__, f2c_real *x, 
	f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sggqrf_(integer *n, integer *m, integer *p, f2c_real *a, 
	integer *lda, f2c_real *taua, f2c_real *b, integer *ldb, f2c_real *taub, f2c_real *
	work, integer *lwork, integer *info);

/* Subroutine */ int sggrqf_(integer *m, integer *p, integer *n, f2c_real *a, 
	integer *lda, f2c_real *taua, f2c_real *b, integer *ldb, f2c_real *taub, f2c_real *
	work, integer *lwork, integer *info);

/* Subroutine */ int sggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, f2c_real *a, integer *lda,
	f2c_real *b, integer *ldb, f2c_real *alpha, f2c_real *beta, f2c_real *u, integer *
	ldu, f2c_real *v, integer *ldv, f2c_real *q, integer *ldq, f2c_real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int sggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, 
	f2c_real *tola, f2c_real *tolb, integer *k, integer *l, f2c_real *u, integer *ldu,
	f2c_real *v, integer *ldv, f2c_real *q, integer *ldq, integer *iwork, f2c_real *
	tau, f2c_real *work, integer *info);

/* Subroutine */ int sgtcon_(char *norm, integer *n, f2c_real *dl, f2c_real *d__, 
	f2c_real *du, f2c_real *du2, integer *ipiv, f2c_real *anorm, f2c_real *rcond, f2c_real *
	work, integer *iwork, integer *info);

/* Subroutine */ int sgtrfs_(char *trans, integer *n, integer *nrhs, f2c_real *dl,
	f2c_real *d__, f2c_real *du, f2c_real *dlf, f2c_real *df, f2c_real *duf, f2c_real *du2, 
	integer *ipiv, f2c_real *b, integer *ldb, f2c_real *x, integer *ldx, f2c_real *
	ferr, f2c_real *berr, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sgtsv_(integer *n, integer *nrhs, f2c_real *dl, f2c_real *d__, 
	f2c_real *du, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, f2c_real *dl, f2c_real *d__, f2c_real *du, f2c_real *dlf, f2c_real *df, f2c_real *duf, 
	f2c_real *du2, integer *ipiv, f2c_real *b, integer *ldb, f2c_real *x, integer *
	ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int sgttrf_(integer *n, f2c_real *dl, f2c_real *d__, f2c_real *du, f2c_real *
	du2, integer *ipiv, integer *info);

/* Subroutine */ int sgttrs_(char *trans, integer *n, integer *nrhs, f2c_real *dl,
	f2c_real *d__, f2c_real *du, f2c_real *du2, integer *ipiv, f2c_real *b, integer *ldb,
	integer *info);

/* Subroutine */ int sgtts2_(integer *itrans, integer *n, integer *nrhs, f2c_real 
	*dl, f2c_real *d__, f2c_real *du, f2c_real *du2, integer *ipiv, f2c_real *b, integer *
	ldb);

/* Subroutine */ int shgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, f2c_real *a, integer *lda, f2c_real *b, integer *
	ldb, f2c_real *alphar, f2c_real *alphai, f2c_real *beta, f2c_real *q, integer *ldq, 
	f2c_real *z__, integer *ldz, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int shsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, f2c_real *h__, integer *ldh, f2c_real *wr, f2c_real *wi, f2c_real 
	*vl, integer *ldvl, f2c_real *vr, integer *ldvr, integer *mm, integer *m, 
	f2c_real *work, integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int shseqr_(char *job, char *compz, integer *n, integer *ilo,
	integer *ihi, f2c_real *h__, integer *ldh, f2c_real *wr, f2c_real *wi, f2c_real *z__,
	integer *ldz, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int slabad_(f2c_real *small, f2c_real *large);

/* Subroutine */ int slabrd_(integer *m, integer *n, integer *nb, f2c_real *a, 
	integer *lda, f2c_real *d__, f2c_real *e, f2c_real *tauq, f2c_real *taup, f2c_real *x, 
	integer *ldx, f2c_real *y, integer *ldy);

/* Subroutine */ int slacon_(integer *n, f2c_real *v, f2c_real *x, integer *isgn, 
	f2c_real *est, integer *kase);

/* Subroutine */ int slacpy_(char *uplo, integer *m, integer *n, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb);

/* Subroutine */ int sladiv_(f2c_real *a, f2c_real *b, f2c_real *c__, f2c_real *d__, f2c_real *p, 
	f2c_real *q);

/* Subroutine */ int slae2_(f2c_real *a, f2c_real *b, f2c_real *c__, f2c_real *rt1, f2c_real *rt2);

/* Subroutine */ int slaebz_(integer *ijob, integer *nitmax, integer *n, 
	integer *mmax, integer *minp, integer *nbmin, f2c_real *abstol, f2c_real *
	reltol, f2c_real *pivmin, f2c_real *d__, f2c_real *e, f2c_real *e2, integer *nval, 
	f2c_real *ab, f2c_real *c__, integer *mout, integer *nab, f2c_real *work, integer 
	*iwork, integer *info);

/* Subroutine */ int slaed0_(integer *icompq, integer *qsiz, integer *n, f2c_real 
	*d__, f2c_real *e, f2c_real *q, integer *ldq, f2c_real *qstore, integer *ldqs, 
	f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int slaed1_(integer *n, f2c_real *d__, f2c_real *q, integer *ldq, 
	integer *indxq, f2c_real *rho, integer *cutpnt, f2c_real *work, integer *
	iwork, integer *info);

/* Subroutine */ int slaed2_(integer *k, integer *n, integer *n1, f2c_real *d__, 
	f2c_real *q, integer *ldq, integer *indxq, f2c_real *rho, f2c_real *z__, f2c_real *
	dlamda, f2c_real *w, f2c_real *q2, integer *indx, integer *indxc, integer *
	indxp, integer *coltyp, integer *info);

/* Subroutine */ int slaed3_(integer *k, integer *n, integer *n1, f2c_real *d__, 
	f2c_real *q, integer *ldq, f2c_real *rho, f2c_real *dlamda, f2c_real *q2, integer *
	indx, integer *ctot, f2c_real *w, f2c_real *s, integer *info);

/* Subroutine */ int slaed4_(integer *n, integer *i__, f2c_real *d__, f2c_real *z__, 
	f2c_real *delta, f2c_real *rho, f2c_real *dlam, integer *info);

/* Subroutine */ int slaed5_(integer *i__, f2c_real *d__, f2c_real *z__, f2c_real *delta, 
	f2c_real *rho, f2c_real *dlam);

/* Subroutine */ int slaed6_(integer *kniter, logical *orgati, f2c_real *rho, 
	f2c_real *d__, f2c_real *z__, f2c_real *finit, f2c_real *tau, integer *info);

/* Subroutine */ int slaed7_(integer *icompq, integer *n, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, f2c_real *d__, f2c_real *q, 
	integer *ldq, integer *indxq, f2c_real *rho, integer *cutpnt, f2c_real *
	qstore, integer *qptr, integer *prmptr, integer *perm, integer *
	givptr, integer *givcol, f2c_real *givnum, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int slaed8_(integer *icompq, integer *k, integer *n, integer 
	*qsiz, f2c_real *d__, f2c_real *q, integer *ldq, integer *indxq, f2c_real *rho, 
	integer *cutpnt, f2c_real *z__, f2c_real *dlamda, f2c_real *q2, integer *ldq2, 
	f2c_real *w, integer *perm, integer *givptr, integer *givcol, f2c_real *
	givnum, integer *indxp, integer *indx, integer *info);

/* Subroutine */ int slaed9_(integer *k, integer *kstart, integer *kstop, 
	integer *n, f2c_real *d__, f2c_real *q, integer *ldq, f2c_real *rho, f2c_real *dlamda,
	f2c_real *w, f2c_real *s, integer *lds, integer *info);

/* Subroutine */ int slaeda_(integer *n, integer *tlvls, integer *curlvl, 
	integer *curpbm, integer *prmptr, integer *perm, integer *givptr, 
	integer *givcol, f2c_real *givnum, f2c_real *q, integer *qptr, f2c_real *z__, 
	f2c_real *ztemp, integer *info);

/* Subroutine */ int slaein_(logical *rightv, logical *noinit, integer *n, 
	f2c_real *h__, integer *ldh, f2c_real *wr, f2c_real *wi, f2c_real *vr, f2c_real *vi, f2c_real 
	*b, integer *ldb, f2c_real *work, f2c_real *eps3, f2c_real *smlnum, f2c_real *bignum, 
	integer *info);

/* Subroutine */ int slaev2_(f2c_real *a, f2c_real *b, f2c_real *c__, f2c_real *rt1, f2c_real *
	rt2, f2c_real *cs1, f2c_real *sn1);

/* Subroutine */ int slaexc_(logical *wantq, integer *n, f2c_real *t, integer *
	ldt, f2c_real *q, integer *ldq, integer *j1, integer *n1, integer *n2, 
	f2c_real *work, integer *info);

/* Subroutine */ int slag2_(f2c_real *a, integer *lda, f2c_real *b, integer *ldb, 
	f2c_real *safmin, f2c_real *scale1, f2c_real *scale2, f2c_real *wr1, f2c_real *wr2, f2c_real *
	wi);

/* Subroutine */ int slags2_(logical *upper, f2c_real *a1, f2c_real *a2, f2c_real *a3, 
	f2c_real *b1, f2c_real *b2, f2c_real *b3, f2c_real *csu, f2c_real *snu, f2c_real *csv, f2c_real *
	snv, f2c_real *csq, f2c_real *snq);

/* Subroutine */ int slagtf_(integer *n, f2c_real *a, f2c_real *lambda, f2c_real *b, f2c_real 
	*c__, f2c_real *tol, f2c_real *d__, integer *in, integer *info);

/* Subroutine */ int slagtm_(char *trans, integer *n, integer *nrhs, f2c_real *
	alpha, f2c_real *dl, f2c_real *d__, f2c_real *du, f2c_real *x, integer *ldx, f2c_real *
	beta, f2c_real *b, integer *ldb);

/* Subroutine */ int slagts_(integer *job, integer *n, f2c_real *a, f2c_real *b, f2c_real 
	*c__, f2c_real *d__, integer *in, f2c_real *y, f2c_real *tol, integer *info);

/* Subroutine */ int slagv2_(f2c_real *a, integer *lda, f2c_real *b, integer *ldb, 
	f2c_real *alphar, f2c_real *alphai, f2c_real *beta, f2c_real *csl, f2c_real *snl, f2c_real *
	csr, f2c_real *snr);

/* Subroutine */ int slahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, f2c_real *h__, integer *ldh, f2c_real *wr, f2c_real *
	wi, integer *iloz, integer *ihiz, f2c_real *z__, integer *ldz, integer *
	info);

/* Subroutine */ int slahrd_(integer *n, integer *k, integer *nb, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *t, integer *ldt, f2c_real *y, integer *ldy);

/* Subroutine */ int slaic1_(integer *job, integer *j, f2c_real *x, f2c_real *sest, 
	f2c_real *w, f2c_real *gamma, f2c_real *sestpr, f2c_real *s, f2c_real *c__);

/* Subroutine */ int slaln2_(logical *ltrans, integer *na, integer *nw, f2c_real *
	smin, f2c_real *ca, f2c_real *a, integer *lda, f2c_real *d1, f2c_real *d2, f2c_real *b, 
	integer *ldb, f2c_real *wr, f2c_real *wi, f2c_real *x, integer *ldx, f2c_real *scale, 
	f2c_real *xnorm, integer *info);

/* Subroutine */ int slals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, f2c_real *b, integer *ldb, f2c_real *bx, 
	integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, f2c_real *givnum, integer *ldgnum, f2c_real *poles, f2c_real *
	difl, f2c_real *difr, f2c_real *z__, integer *k, f2c_real *c__, f2c_real *s, f2c_real *
	work, integer *info);

/* Subroutine */ int slalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, f2c_real *b, integer *ldb, f2c_real *bx, integer *ldbx, f2c_real *
	u, integer *ldu, f2c_real *vt, integer *k, f2c_real *difl, f2c_real *difr, f2c_real *
	z__, f2c_real *poles, integer *givptr, integer *givcol, integer *ldgcol, 
	integer *perm, f2c_real *givnum, f2c_real *c__, f2c_real *s, f2c_real *work, integer *
	iwork, integer *info);

/* Subroutine */ int slalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, f2c_real *d__, f2c_real *e, f2c_real *b, integer *ldb, f2c_real *rcond, 
	integer *rank, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int slamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1);

/* Subroutine */ int slamc2_(integer *beta, integer *t, logical *rnd, f2c_real *
	eps, integer *emin, f2c_real *rmin, integer *emax, f2c_real *rmax);

/* Subroutine */ int slamc4_(integer *emin, f2c_real *start, integer *base);

/* Subroutine */ int slamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, f2c_real *rmax);

/* Subroutine */ int slamrg_(integer *n1, integer *n2, f2c_real *a, integer *
	strd1, integer *strd2, integer *index);

/* Subroutine */ int slanv2_(f2c_real *a, f2c_real *b, f2c_real *c__, f2c_real *d__, f2c_real *
	rt1r, f2c_real *rt1i, f2c_real *rt2r, f2c_real *rt2i, f2c_real *cs, f2c_real *sn);

/* Subroutine */ int slapll_(integer *n, f2c_real *x, integer *incx, f2c_real *y, 
	integer *incy, f2c_real *ssmin);

/* Subroutine */ int slapmt_(logical *forwrd, integer *m, integer *n, f2c_real *x,
	integer *ldx, integer *k);

/* Subroutine */ int slaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	f2c_real *ab, integer *ldab, f2c_real *r__, f2c_real *c__, f2c_real *rowcnd, f2c_real *
	colcnd, f2c_real *amax, char *equed);

/* Subroutine */ int slaqge_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *r__, f2c_real *c__, f2c_real *rowcnd, f2c_real *colcnd, f2c_real *amax, char *
	equed);

/* Subroutine */ int slaqp2_(integer *m, integer *n, integer *offset, f2c_real *a,
	integer *lda, integer *jpvt, f2c_real *tau, f2c_real *vn1, f2c_real *vn2, f2c_real *
	work);

/* Subroutine */ int slaqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, f2c_real *a, integer *lda, integer *jpvt, f2c_real *tau, 
	f2c_real *vn1, f2c_real *vn2, f2c_real *auxv, f2c_real *f, integer *ldf);

/* Subroutine */ int slaqsb_(char *uplo, integer *n, integer *kd, f2c_real *ab, 
	integer *ldab, f2c_real *s, f2c_real *scond, f2c_real *amax, char *equed);

/* Subroutine */ int slaqsp_(char *uplo, integer *n, f2c_real *ap, f2c_real *s, f2c_real *
	scond, f2c_real *amax, char *equed);

/* Subroutine */ int slaqsy_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	f2c_real *s, f2c_real *scond, f2c_real *amax, char *equed);

/* Subroutine */ int slaqtr_(logical *ltran, logical *lf2c_real, integer *n, f2c_real 
	*t, integer *ldt, f2c_real *b, f2c_real *w, f2c_real *scale, f2c_real *x, f2c_real *work, 
	integer *info);

/* Subroutine */ int slar1v_(integer *n, integer *b1, integer *bn, f2c_real *
	sigma, f2c_real *d__, f2c_real *l, f2c_real *ld, f2c_real *lld, f2c_real *gersch, f2c_real *
	z__, f2c_real *ztz, f2c_real *mingma, integer *r__, integer *isuppz, f2c_real *
	work);

/* Subroutine */ int slar2v_(integer *n, f2c_real *x, f2c_real *y, f2c_real *z__, integer 
	*incx, f2c_real *c__, f2c_real *s, integer *incc);

/* Subroutine */ int slarf_(char *side, integer *m, integer *n, f2c_real *v, 
	integer *incv, f2c_real *tau, f2c_real *c__, integer *ldc, f2c_real *work);

/* Subroutine */ int slarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, f2c_real *v, integer *ldv, 
	f2c_real *t, integer *ldt, f2c_real *c__, integer *ldc, f2c_real *work, integer *
	ldwork);

/* Subroutine */ int slarfg_(integer *n, f2c_real *alpha, f2c_real *x, integer *incx, 
	f2c_real *tau);

/* Subroutine */ int slarft_(char *direct, char *storev, integer *n, integer *
	k, f2c_real *v, integer *ldv, f2c_real *tau, f2c_real *t, integer *ldt);

/* Subroutine */ int slarfx_(char *side, integer *m, integer *n, f2c_real *v, 
	f2c_real *tau, f2c_real *c__, integer *ldc, f2c_real *work);

/* Subroutine */ int slargv_(integer *n, f2c_real *x, integer *incx, f2c_real *y, 
	integer *incy, f2c_real *c__, integer *incc);

/* Subroutine */ int slarnv_(integer *idist, integer *iseed, integer *n, f2c_real 
	*x);

/* Subroutine */ int slarrb_(integer *n, f2c_real *d__, f2c_real *l, f2c_real *ld, f2c_real *
	lld, integer *ifirst, integer *ilast, f2c_real *sigma, f2c_real *reltol, f2c_real 
	*w, f2c_real *wgap, f2c_real *werr, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int slarre_(integer *n, f2c_real *d__, f2c_real *e, f2c_real *tol, 
	integer *nsplit, integer *isplit, integer *m, f2c_real *w, f2c_real *woff, 
	f2c_real *gersch, f2c_real *work, integer *info);

/* Subroutine */ int slarrf_(integer *n, f2c_real *d__, f2c_real *l, f2c_real *ld, f2c_real *
	lld, integer *ifirst, integer *ilast, f2c_real *w, f2c_real *dplus, f2c_real *
	lplus, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int slarrv_(integer *n, f2c_real *d__, f2c_real *l, integer *isplit, 
	integer *m, f2c_real *w, integer *iblock, f2c_real *gersch, f2c_real *tol, f2c_real *
	z__, integer *ldz, integer *isuppz, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int slartg_(f2c_real *f, f2c_real *g, f2c_real *cs, f2c_real *sn, f2c_real *r__);

/* Subroutine */ int slartv_(integer *n, f2c_real *x, integer *incx, f2c_real *y, 
	integer *incy, f2c_real *c__, f2c_real *s, integer *incc);

/* Subroutine */ int slaruv_(integer *iseed, integer *n, f2c_real *x);

/* Subroutine */ int slarz_(char *side, integer *m, integer *n, integer *l, 
	f2c_real *v, integer *incv, f2c_real *tau, f2c_real *c__, integer *ldc, f2c_real *
	work);

/* Subroutine */ int slarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, f2c_real *v, 
	integer *ldv, f2c_real *t, integer *ldt, f2c_real *c__, integer *ldc, f2c_real *
	work, integer *ldwork);

/* Subroutine */ int slarzt_(char *direct, char *storev, integer *n, integer *
	k, f2c_real *v, integer *ldv, f2c_real *tau, f2c_real *t, integer *ldt);

/* Subroutine */ int slas2_(f2c_real *f, f2c_real *g, f2c_real *h__, f2c_real *ssmin, f2c_real *
	ssmax);

/* Subroutine */ int slascl_(char *type__, integer *kl, integer *ku, f2c_real *
	cfrom, f2c_real *cto, integer *m, integer *n, f2c_real *a, integer *lda, 
	integer *info);

/* Subroutine */ int slasd0_(integer *n, integer *sqre, f2c_real *d__, f2c_real *e, 
	f2c_real *u, integer *ldu, f2c_real *vt, integer *ldvt, integer *smlsiz, 
	integer *iwork, f2c_real *work, integer *info);

/* Subroutine */ int slasd1_(integer *nl, integer *nr, integer *sqre, f2c_real *
	d__, f2c_real *alpha, f2c_real *beta, f2c_real *u, integer *ldu, f2c_real *vt, 
	integer *ldvt, integer *idxq, integer *iwork, f2c_real *work, integer *
	info);

/* Subroutine */ int slasd2_(integer *nl, integer *nr, integer *sqre, integer 
	*k, f2c_real *d__, f2c_real *z__, f2c_real *alpha, f2c_real *beta, f2c_real *u, integer *
	ldu, f2c_real *vt, integer *ldvt, f2c_real *dsigma, f2c_real *u2, integer *ldu2, 
	f2c_real *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc,
	integer *idxq, integer *coltyp, integer *info);

/* Subroutine */ int slasd3_(integer *nl, integer *nr, integer *sqre, integer 
	*k, f2c_real *d__, f2c_real *q, integer *ldq, f2c_real *dsigma, f2c_real *u, integer *
	ldu, f2c_real *u2, integer *ldu2, f2c_real *vt, integer *ldvt, f2c_real *vt2, 
	integer *ldvt2, integer *idxc, integer *ctot, f2c_real *z__, integer *
	info);

/* Subroutine */ int slasd4_(integer *n, integer *i__, f2c_real *d__, f2c_real *z__, 
	f2c_real *delta, f2c_real *rho, f2c_real *sigma, f2c_real *work, integer *info);

/* Subroutine */ int slasd5_(integer *i__, f2c_real *d__, f2c_real *z__, f2c_real *delta, 
	f2c_real *rho, f2c_real *dsigma, f2c_real *work);

/* Subroutine */ int slasd6_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, f2c_real *d__, f2c_real *vf, f2c_real *vl, f2c_real *alpha, f2c_real *beta,
	integer *idxq, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, f2c_real *givnum, integer *ldgnum, f2c_real *poles, f2c_real *
	difl, f2c_real *difr, f2c_real *z__, integer *k, f2c_real *c__, f2c_real *s, f2c_real *
	work, integer *iwork, integer *info);

/* Subroutine */ int slasd7_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *k, f2c_real *d__, f2c_real *z__, f2c_real *zw, f2c_real *vf, 
	f2c_real *vfw, f2c_real *vl, f2c_real *vlw, f2c_real *alpha, f2c_real *beta, f2c_real *dsigma,
	integer *idx, integer *idxp, integer *idxq, integer *perm, integer *
	givptr, integer *givcol, integer *ldgcol, f2c_real *givnum, integer *
	ldgnum, f2c_real *c__, f2c_real *s, integer *info);

/* Subroutine */ int slasd8_(integer *icompq, integer *k, f2c_real *d__, f2c_real *
	z__, f2c_real *vf, f2c_real *vl, f2c_real *difl, f2c_real *difr, integer *lddifr, 
	f2c_real *dsigma, f2c_real *work, integer *info);

/* Subroutine */ int slasd9_(integer *icompq, integer *ldu, integer *k, f2c_real *
	d__, f2c_real *z__, f2c_real *vf, f2c_real *vl, f2c_real *difl, f2c_real *difr, f2c_real *
	dsigma, f2c_real *work, integer *info);

/* Subroutine */ int slasda_(integer *icompq, integer *smlsiz, integer *n, 
	integer *sqre, f2c_real *d__, f2c_real *e, f2c_real *u, integer *ldu, f2c_real *vt, 
	integer *k, f2c_real *difl, f2c_real *difr, f2c_real *z__, f2c_real *poles, integer *
	givptr, integer *givcol, integer *ldgcol, integer *perm, f2c_real *givnum,
	f2c_real *c__, f2c_real *s, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int slasdq_(char *uplo, integer *sqre, integer *n, integer *
	ncvt, integer *nru, integer *ncc, f2c_real *d__, f2c_real *e, f2c_real *vt, 
	integer *ldvt, f2c_real *u, integer *ldu, f2c_real *c__, integer *ldc, f2c_real *
	work, integer *info);

/* Subroutine */ int slasdt_(integer *n, integer *lvl, integer *nd, integer *
	inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ int slaset_(char *uplo, integer *m, integer *n, f2c_real *alpha, 
	f2c_real *beta, f2c_real *a, integer *lda);

/* Subroutine */ int slasq1_(integer *n, f2c_real *d__, f2c_real *e, f2c_real *work, 
	integer *info);

/* Subroutine */ int slasq2_(integer *n, f2c_real *z__, integer *info);

/* Subroutine */ int slasq3_(integer *i0, integer *n0, f2c_real *z__, integer *pp,
	f2c_real *dmin__, f2c_real *sigma, f2c_real *desig, f2c_real *qmax, integer *nfail, 
	integer *iter, integer *ndiv, logical *ieee);

/* Subroutine */ int slasq4_(integer *i0, integer *n0, f2c_real *z__, integer *pp,
	integer *n0in, f2c_real *dmin__, f2c_real *dmin1, f2c_real *dmin2, f2c_real *dn, 
	f2c_real *dn1, f2c_real *dn2, f2c_real *tau, integer *ttype);

/* Subroutine */ int slasq5_(integer *i0, integer *n0, f2c_real *z__, integer *pp,
	f2c_real *tau, f2c_real *dmin__, f2c_real *dmin1, f2c_real *dmin2, f2c_real *dn, f2c_real *
	dnm1, f2c_real *dnm2, logical *ieee);

/* Subroutine */ int slasq6_(integer *i0, integer *n0, f2c_real *z__, integer *pp,
	f2c_real *dmin__, f2c_real *dmin1, f2c_real *dmin2, f2c_real *dn, f2c_real *dnm1, f2c_real *
	dnm2);

/* Subroutine */ int slasr_(char *side, char *pivot, char *direct, integer *m,
	integer *n, f2c_real *c__, f2c_real *s, f2c_real *a, integer *lda);

/* Subroutine */ int slasrt_(char *id, integer *n, f2c_real *d__, integer *info);

/* Subroutine */ int slassq_(integer *n, f2c_real *x, integer *incx, f2c_real *scale, 
	f2c_real *sumsq);

/* Subroutine */ int slasv2_(f2c_real *f, f2c_real *g, f2c_real *h__, f2c_real *ssmin, f2c_real *
	ssmax, f2c_real *snr, f2c_real *csr, f2c_real *snl, f2c_real *csl);

/* Subroutine */ int slaswp_(integer *n, f2c_real *a, integer *lda, integer *k1, 
	integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int slasy2_(logical *ltranl, logical *ltranr, integer *isgn, 
	integer *n1, integer *n2, f2c_real *tl, integer *ldtl, f2c_real *tr, integer *
	ldtr, f2c_real *b, integer *ldb, f2c_real *scale, f2c_real *x, integer *ldx, f2c_real 
	*xnorm, integer *info);

/* Subroutine */ int slasyf_(char *uplo, integer *n, integer *nb, integer *kb,
	f2c_real *a, integer *lda, integer *ipiv, f2c_real *w, integer *ldw, integer 
	*info);

/* Subroutine */ int slatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, f2c_real *ab, integer *ldab, f2c_real *x, 
	f2c_real *scale, f2c_real *cnorm, integer *info);

/* Subroutine */ int slatdf_(integer *ijob, integer *n, f2c_real *z__, integer *
	ldz, f2c_real *rhs, f2c_real *rdsum, f2c_real *rdscal, integer *ipiv, integer *
	jpiv);

/* Subroutine */ int slatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, f2c_real *ap, f2c_real *x, f2c_real *scale, f2c_real *cnorm, 
	integer *info);

/* Subroutine */ int slatrd_(char *uplo, integer *n, integer *nb, f2c_real *a, 
	integer *lda, f2c_real *e, f2c_real *tau, f2c_real *w, integer *ldw);

/* Subroutine */ int slatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, f2c_real *a, integer *lda, f2c_real *x, f2c_real *scale, f2c_real 
	*cnorm, integer *info);

/* Subroutine */ int slatrz_(integer *m, integer *n, integer *l, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work);

/* Subroutine */ int slatzm_(char *side, integer *m, integer *n, f2c_real *v, 
	integer *incv, f2c_real *tau, f2c_real *c1, f2c_real *c2, integer *ldc, f2c_real *
	work);

/* Subroutine */ int slauu2_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *info);

/* Subroutine */ int slauum_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *info);

/* Subroutine */ int sopgtr_(char *uplo, integer *n, f2c_real *ap, f2c_real *tau, 
	f2c_real *q, integer *ldq, f2c_real *work, integer *info);

/* Subroutine */ int sopmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, f2c_real *ap, f2c_real *tau, f2c_real *c__, integer *ldc, f2c_real *work, 
	integer *info);

/* Subroutine */ int sorg2l_(integer *m, integer *n, integer *k, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sorg2r_(integer *m, integer *n, integer *k, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sorgbr_(char *vect, integer *m, integer *n, integer *k, 
	f2c_real *a, integer *lda, f2c_real *tau, f2c_real *work, integer *lwork, integer 
	*info);

/* Subroutine */ int sorghr_(integer *n, integer *ilo, integer *ihi, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgl2_(integer *m, integer *n, integer *k, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sorglq_(integer *m, integer *n, integer *k, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgql_(integer *m, integer *n, integer *k, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgqr_(integer *m, integer *n, integer *k, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgr2_(integer *m, integer *n, integer *k, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *info);

/* Subroutine */ int sorgrq_(integer *m, integer *n, integer *k, f2c_real *a, 
	integer *lda, f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgtr_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sorm2l_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *info);

/* Subroutine */ int sorm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *info);

/* Subroutine */ int sormbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, 
	integer *ldc, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sormhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *
	c__, integer *ldc, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sorml2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *info);

/* Subroutine */ int sormlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sormql_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sormqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sormr2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *info);

/* Subroutine */ int sormr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, 
	integer *ldc, f2c_real *work, integer *info);

/* Subroutine */ int sormrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sormrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, 
	integer *ldc, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int sormtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, f2c_real *a, integer *lda, f2c_real *tau, f2c_real *c__, integer *ldc,
	f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int spbcon_(char *uplo, integer *n, integer *kd, f2c_real *ab, 
	integer *ldab, f2c_real *anorm, f2c_real *rcond, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int spbequ_(char *uplo, integer *n, integer *kd, f2c_real *ab, 
	integer *ldab, f2c_real *s, f2c_real *scond, f2c_real *amax, integer *info);

/* Subroutine */ int spbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, f2c_real *ab, integer *ldab, f2c_real *afb, integer *ldafb, f2c_real *b, 
	integer *ldb, f2c_real *x, integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_real *
	work, integer *iwork, integer *info);

/* Subroutine */ int spbstf_(char *uplo, integer *n, integer *kd, f2c_real *ab, 
	integer *ldab, integer *info);

/* Subroutine */ int spbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, f2c_real *ab, integer *ldab, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int spbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, f2c_real *ab, integer *ldab, f2c_real *afb, integer *ldafb, 
	char *equed, f2c_real *s, f2c_real *b, integer *ldb, f2c_real *x, integer *ldx, 
	f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int spbtf2_(char *uplo, integer *n, integer *kd, f2c_real *ab, 
	integer *ldab, integer *info);

/* Subroutine */ int spbtrf_(char *uplo, integer *n, integer *kd, f2c_real *ab, 
	integer *ldab, integer *info);

/* Subroutine */ int spbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, f2c_real *ab, integer *ldab, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int spocon_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	f2c_real *anorm, f2c_real *rcond, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int spoequ_(integer *n, f2c_real *a, integer *lda, f2c_real *s, f2c_real 
	*scond, f2c_real *amax, integer *info);

/* Subroutine */ int sporfs_(char *uplo, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *af, integer *ldaf, f2c_real *b, integer *ldb, f2c_real *x,
	integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int sposv_(char *uplo, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_real *a, integer *lda, f2c_real *af, integer *ldaf, char *equed, 
	f2c_real *s, f2c_real *b, integer *ldb, f2c_real *x, integer *ldx, f2c_real *rcond, 
	f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int spotf2_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *info);

/* Subroutine */ int spotrf_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *info);

/* Subroutine */ int spotri_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *info);

/* Subroutine */ int spotrs_(char *uplo, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sppcon_(char *uplo, integer *n, f2c_real *ap, f2c_real *anorm, 
	f2c_real *rcond, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sppequ_(char *uplo, integer *n, f2c_real *ap, f2c_real *s, f2c_real *
	scond, f2c_real *amax, integer *info);

/* Subroutine */ int spprfs_(char *uplo, integer *n, integer *nrhs, f2c_real *ap, 
	f2c_real *afp, f2c_real *b, integer *ldb, f2c_real *x, integer *ldx, f2c_real *ferr, 
	f2c_real *berr, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sppsv_(char *uplo, integer *n, integer *nrhs, f2c_real *ap, 
	f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_real *ap, f2c_real *afp, char *equed, f2c_real *s, f2c_real *b, integer *
	ldb, f2c_real *x, integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_real 
	*work, integer *iwork, integer *info);

/* Subroutine */ int spptrf_(char *uplo, integer *n, f2c_real *ap, integer *info);

/* Subroutine */ int spptri_(char *uplo, integer *n, f2c_real *ap, integer *info);

/* Subroutine */ int spptrs_(char *uplo, integer *n, integer *nrhs, f2c_real *ap, 
	f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sptcon_(integer *n, f2c_real *d__, f2c_real *e, f2c_real *anorm, 
	f2c_real *rcond, f2c_real *work, integer *info);

/* Subroutine */ int spteqr_(char *compz, integer *n, f2c_real *d__, f2c_real *e, 
	f2c_real *z__, integer *ldz, f2c_real *work, integer *info);

/* Subroutine */ int sptrfs_(integer *n, integer *nrhs, f2c_real *d__, f2c_real *e, 
	f2c_real *df, f2c_real *ef, f2c_real *b, integer *ldb, f2c_real *x, integer *ldx, 
	f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *info);

/* Subroutine */ int sptsv_(integer *n, integer *nrhs, f2c_real *d__, f2c_real *e, 
	f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sptsvx_(char *fact, integer *n, integer *nrhs, f2c_real *d__,
	f2c_real *e, f2c_real *df, f2c_real *ef, f2c_real *b, integer *ldb, f2c_real *x, integer 
	*ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *info);

/* Subroutine */ int spttrf_(integer *n, f2c_real *d__, f2c_real *e, integer *info);

/* Subroutine */ int spttrs_(integer *n, integer *nrhs, f2c_real *d__, f2c_real *e, 
	f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sptts2_(integer *n, integer *nrhs, f2c_real *d__, f2c_real *e, 
	f2c_real *b, integer *ldb);

/* Subroutine */ int srscl_(integer *n, f2c_real *sa, f2c_real *sx, integer *incx);

/* Subroutine */ int ssbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	f2c_real *ab, integer *ldab, f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work,
	integer *info);

/* Subroutine */ int ssbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	f2c_real *ab, integer *ldab, f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int ssbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, f2c_real *ab, integer *ldab, f2c_real *q, integer *ldq, f2c_real *vl,
	f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, integer *m, f2c_real *
	w, f2c_real *z__, integer *ldz, f2c_real *work, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int ssbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, f2c_real *ab, integer *ldab, f2c_real *bb, integer *ldbb, f2c_real *
	x, integer *ldx, f2c_real *work, integer *info);

/* Subroutine */ int ssbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, f2c_real *ab, integer *ldab, f2c_real *bb, integer *ldbb, f2c_real *
	w, f2c_real *z__, integer *ldz, f2c_real *work, integer *info);

/* Subroutine */ int ssbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, f2c_real *ab, integer *ldab, f2c_real *bb, integer *ldbb, f2c_real *
	w, f2c_real *z__, integer *ldz, f2c_real *work, integer *lwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int ssbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, f2c_real *ab, integer *ldab, f2c_real *bb, integer *
	ldbb, f2c_real *q, integer *ldq, f2c_real *vl, f2c_real *vu, integer *il, integer 
	*iu, f2c_real *abstol, integer *m, f2c_real *w, f2c_real *z__, integer *ldz, f2c_real 
	*work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int ssbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	f2c_real *ab, integer *ldab, f2c_real *d__, f2c_real *e, f2c_real *q, integer *ldq, 
	f2c_real *work, integer *info);

/* Subroutine */ int sspcon_(char *uplo, integer *n, f2c_real *ap, integer *ipiv, 
	f2c_real *anorm, f2c_real *rcond, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sspev_(char *jobz, char *uplo, integer *n, f2c_real *ap, 
	f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work, integer *info);

/* Subroutine */ int sspevd_(char *jobz, char *uplo, integer *n, f2c_real *ap, 
	f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work, integer *lwork, integer 
	*iwork, integer *liwork, integer *info);

/* Subroutine */ int sspevx_(char *jobz, char *range, char *uplo, integer *n, 
	f2c_real *ap, f2c_real *vl, f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, 
	integer *m, f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work, integer *
	iwork, integer *ifail, integer *info);

/* Subroutine */ int sspgst_(integer *itype, char *uplo, integer *n, f2c_real *ap,
	f2c_real *bp, integer *info);

/* Subroutine */ int sspgv_(integer *itype, char *jobz, char *uplo, integer *
	n, f2c_real *ap, f2c_real *bp, f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work, 
	integer *info);

/* Subroutine */ int sspgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, f2c_real *ap, f2c_real *bp, f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int sspgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, f2c_real *ap, f2c_real *bp, f2c_real *vl, f2c_real *vu, integer *il,
	integer *iu, f2c_real *abstol, integer *m, f2c_real *w, f2c_real *z__, integer *
	ldz, f2c_real *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int ssprfs_(char *uplo, integer *n, integer *nrhs, f2c_real *ap, 
	f2c_real *afp, integer *ipiv, f2c_real *b, integer *ldb, f2c_real *x, integer *
	ldx, f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *iwork, integer *
	info);

/* Subroutine */ int sspsv_(char *uplo, integer *n, integer *nrhs, f2c_real *ap, 
	integer *ipiv, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_real *ap, f2c_real *afp, integer *ipiv, f2c_real *b, integer *ldb, f2c_real 
	*x, integer *ldx, f2c_real *rcond, f2c_real *ferr, f2c_real *berr, f2c_real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int ssptrd_(char *uplo, integer *n, f2c_real *ap, f2c_real *d__, 
	f2c_real *e, f2c_real *tau, integer *info);

/* Subroutine */ int ssptrf_(char *uplo, integer *n, f2c_real *ap, integer *ipiv, 
	integer *info);

/* Subroutine */ int ssptri_(char *uplo, integer *n, f2c_real *ap, integer *ipiv, 
	f2c_real *work, integer *info);

/* Subroutine */ int ssptrs_(char *uplo, integer *n, integer *nrhs, f2c_real *ap, 
	integer *ipiv, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int sstebz_(char *range, char *order, integer *n, f2c_real *vl, 
	f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, f2c_real *d__, f2c_real *e, 
	integer *m, integer *nsplit, f2c_real *w, integer *iblock, integer *
	isplit, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int sstedc_(char *compz, integer *n, f2c_real *d__, f2c_real *e, 
	f2c_real *z__, integer *ldz, f2c_real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int sstegr_(char *jobz, char *range, integer *n, f2c_real *d__, 
	f2c_real *e, f2c_real *vl, f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, 
	integer *m, f2c_real *w, f2c_real *z__, integer *ldz, integer *isuppz, f2c_real *
	work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int sstein_(integer *n, f2c_real *d__, f2c_real *e, integer *m, f2c_real 
	*w, integer *iblock, integer *isplit, f2c_real *z__, integer *ldz, f2c_real *
	work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int ssteqr_(char *compz, integer *n, f2c_real *d__, f2c_real *e, 
	f2c_real *z__, integer *ldz, f2c_real *work, integer *info);

/* Subroutine */ int ssterf_(integer *n, f2c_real *d__, f2c_real *e, integer *info);

/* Subroutine */ int sstev_(char *jobz, integer *n, f2c_real *d__, f2c_real *e, f2c_real *
	z__, integer *ldz, f2c_real *work, integer *info);

/* Subroutine */ int sstevd_(char *jobz, integer *n, f2c_real *d__, f2c_real *e, f2c_real 
	*z__, integer *ldz, f2c_real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int sstevr_(char *jobz, char *range, integer *n, f2c_real *d__, 
	f2c_real *e, f2c_real *vl, f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, 
	integer *m, f2c_real *w, f2c_real *z__, integer *ldz, integer *isuppz, f2c_real *
	work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int sstevx_(char *jobz, char *range, integer *n, f2c_real *d__, 
	f2c_real *e, f2c_real *vl, f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, 
	integer *m, f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work, integer *
	iwork, integer *ifail, integer *info);

/* Subroutine */ int ssycon_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *ipiv, f2c_real *anorm, f2c_real *rcond, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int ssyev_(char *jobz, char *uplo, integer *n, f2c_real *a, 
	integer *lda, f2c_real *w, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int ssyevd_(char *jobz, char *uplo, integer *n, f2c_real *a, 
	integer *lda, f2c_real *w, f2c_real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int ssyevr_(char *jobz, char *range, char *uplo, integer *n, 
	f2c_real *a, integer *lda, f2c_real *vl, f2c_real *vu, integer *il, integer *iu, 
	f2c_real *abstol, integer *m, f2c_real *w, f2c_real *z__, integer *ldz, integer *
	isuppz, f2c_real *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int ssyevx_(char *jobz, char *range, char *uplo, integer *n, 
	f2c_real *a, integer *lda, f2c_real *vl, f2c_real *vu, integer *il, integer *iu, 
	f2c_real *abstol, integer *m, f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *
	work, integer *lwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int ssygs2_(integer *itype, char *uplo, integer *n, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int ssygst_(integer *itype, char *uplo, integer *n, f2c_real *a, 
	integer *lda, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int ssygv_(integer *itype, char *jobz, char *uplo, integer *
	n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *w, f2c_real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int ssygvd_(integer *itype, char *jobz, char *uplo, integer *
	n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *w, f2c_real *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int ssygvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *
	vl, f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, integer *m, 
	f2c_real *w, f2c_real *z__, integer *ldz, f2c_real *work, integer *lwork, integer 
	*iwork, integer *ifail, integer *info);

/* Subroutine */ int ssyrfs_(char *uplo, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, f2c_real *af, integer *ldaf, integer *ipiv, f2c_real *b, 
	integer *ldb, f2c_real *x, integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_real *
	work, integer *iwork, integer *info);

/* Subroutine */ int ssysv_(char *uplo, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, integer *ipiv, f2c_real *b, integer *ldb, f2c_real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int ssysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, f2c_real *a, integer *lda, f2c_real *af, integer *ldaf, integer *ipiv, 
	f2c_real *b, integer *ldb, f2c_real *x, integer *ldx, f2c_real *rcond, f2c_real *ferr,
	f2c_real *berr, f2c_real *work, integer *lwork, integer *iwork, integer *
	info);

/* Subroutine */ int ssytd2_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	f2c_real *d__, f2c_real *e, f2c_real *tau, integer *info);

/* Subroutine */ int ssytf2_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *ipiv, integer *info);

/* Subroutine */ int ssytrd_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	f2c_real *d__, f2c_real *e, f2c_real *tau, f2c_real *work, integer *lwork, integer *
	info);

/* Subroutine */ int ssytrf_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *ipiv, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int ssytri_(char *uplo, integer *n, f2c_real *a, integer *lda, 
	integer *ipiv, f2c_real *work, integer *info);

/* Subroutine */ int ssytrs_(char *uplo, integer *n, integer *nrhs, f2c_real *a, 
	integer *lda, integer *ipiv, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int stbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, f2c_real *ab, integer *ldab, f2c_real *rcond, f2c_real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int stbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, f2c_real *ab, integer *ldab, f2c_real *b, integer 
	*ldb, f2c_real *x, integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int stbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, f2c_real *ab, integer *ldab, f2c_real *b, integer 
	*ldb, integer *info);

/* Subroutine */ int stgevc_(char *side, char *howmny, logical *select, 
	integer *n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *vl, 
	integer *ldvl, f2c_real *vr, integer *ldvr, integer *mm, integer *m, f2c_real 
	*work, integer *info);

/* Subroutine */ int stgex2_(logical *wantq, logical *wantz, integer *n, f2c_real 
	*a, integer *lda, f2c_real *b, integer *ldb, f2c_real *q, integer *ldq, f2c_real *
	z__, integer *ldz, integer *j1, integer *n1, integer *n2, f2c_real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int stgexc_(logical *wantq, logical *wantz, integer *n, f2c_real 
	*a, integer *lda, f2c_real *b, integer *ldb, f2c_real *q, integer *ldq, f2c_real *
	z__, integer *ldz, integer *ifst, integer *ilst, f2c_real *work, integer *
	lwork, integer *info);

/* Subroutine */ int stgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, f2c_real *a, integer *lda, f2c_real *b, integer *
	ldb, f2c_real *alphar, f2c_real *alphai, f2c_real *beta, f2c_real *q, integer *ldq, 
	f2c_real *z__, integer *ldz, integer *m, f2c_real *pl, f2c_real *pr, f2c_real *dif, 
	f2c_real *work, integer *lwork, integer *iwork, integer *liwork, integer *
	info);

/* Subroutine */ int stgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, f2c_real *a, integer *lda,
	f2c_real *b, integer *ldb, f2c_real *tola, f2c_real *tolb, f2c_real *alpha, f2c_real *
	beta, f2c_real *u, integer *ldu, f2c_real *v, integer *ldv, f2c_real *q, integer *
	ldq, f2c_real *work, integer *ncycle, integer *info);

/* Subroutine */ int stgsna_(char *job, char *howmny, logical *select, 
	integer *n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *vl, 
	integer *ldvl, f2c_real *vr, integer *ldvr, f2c_real *s, f2c_real *dif, integer *
	mm, integer *m, f2c_real *work, integer *lwork, integer *iwork, integer *
	info);

/* Subroutine */ int stgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *c__, integer *
	ldc, f2c_real *d__, integer *ldd, f2c_real *e, integer *lde, f2c_real *f, integer 
	*ldf, f2c_real *scale, f2c_real *rdsum, f2c_real *rdscal, integer *iwork, integer 
	*pq, integer *info);

/* Subroutine */ int stgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *c__, integer *
	ldc, f2c_real *d__, integer *ldd, f2c_real *e, integer *lde, f2c_real *f, integer 
	*ldf, f2c_real *scale, f2c_real *dif, f2c_real *work, integer *lwork, integer *
	iwork, integer *info);

/* Subroutine */ int stpcon_(char *norm, char *uplo, char *diag, integer *n, 
	f2c_real *ap, f2c_real *rcond, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int stprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, f2c_real *ap, f2c_real *b, integer *ldb, f2c_real *x, integer *ldx,
	f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *iwork, integer *info);

/* Subroutine */ int stptri_(char *uplo, char *diag, integer *n, f2c_real *ap, 
	integer *info);

/* Subroutine */ int stptrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, f2c_real *ap, f2c_real *b, integer *ldb, integer *info);

/* Subroutine */ int strcon_(char *norm, char *uplo, char *diag, integer *n, 
	f2c_real *a, integer *lda, f2c_real *rcond, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int strevc_(char *side, char *howmny, logical *select, 
	integer *n, f2c_real *t, integer *ldt, f2c_real *vl, integer *ldvl, f2c_real *vr, 
	integer *ldvr, integer *mm, integer *m, f2c_real *work, integer *info);

/* Subroutine */ int strexc_(char *compq, integer *n, f2c_real *t, integer *ldt, 
	f2c_real *q, integer *ldq, integer *ifst, integer *ilst, f2c_real *work, 
	integer *info);

/* Subroutine */ int strrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *x, 
	integer *ldx, f2c_real *ferr, f2c_real *berr, f2c_real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int strsen_(char *job, char *compq, logical *select, integer 
	*n, f2c_real *t, integer *ldt, f2c_real *q, integer *ldq, f2c_real *wr, f2c_real *wi, 
	integer *m, f2c_real *s, f2c_real *sep, f2c_real *work, integer *lwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int strsna_(char *job, char *howmny, logical *select, 
	integer *n, f2c_real *t, integer *ldt, f2c_real *vl, integer *ldvl, f2c_real *vr, 
	integer *ldvr, f2c_real *s, f2c_real *sep, integer *mm, integer *m, f2c_real *
	work, integer *ldwork, integer *iwork, integer *info);

/* Subroutine */ int strsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, f2c_real *
	c__, integer *ldc, f2c_real *scale, integer *info);

/* Subroutine */ int strti2_(char *uplo, char *diag, integer *n, f2c_real *a, 
	integer *lda, integer *info);

/* Subroutine */ int strtri_(char *uplo, char *diag, integer *n, f2c_real *a, 
	integer *lda, integer *info);

/* Subroutine */ int strtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, f2c_real *a, integer *lda, f2c_real *b, integer *ldb, integer *
	info);

/* Subroutine */ int stzrqf_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, integer *info);

/* Subroutine */ int stzrzf_(integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *tau, f2c_real *work, integer *lwork, integer *info);

/* Subroutine */ int xerbla_(char *srname, integer *info);

/* Subroutine */ int zbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, doublereal *d__, doublereal *e, doublecomplex *vt, 
	integer *ldvt, doublecomplex *u, integer *ldu, doublecomplex *c__, 
	integer *ldc, doublereal *rwork, integer *info);

/* Subroutine */ int zdrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublereal *s);

/* Subroutine */ int zdrscl_(integer *n, doublereal *sa, doublecomplex *sx, 
	integer *incx);

/* Subroutine */ int zgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	integer *kl, integer *ku, doublecomplex *ab, integer *ldab, 
	doublereal *d__, doublereal *e, doublecomplex *q, integer *ldq, 
	doublecomplex *pt, integer *ldpt, doublecomplex *c__, integer *ldc, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	doublecomplex *ab, integer *ldab, integer *ipiv, doublereal *anorm, 
	doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info);

/* Subroutine */ int zgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *
	afb, integer *ldafb, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zgbsv_(integer *n, integer *kl, integer *ku, integer *
	nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, doublecomplex *
	b, integer *ldb, integer *info);

/* Subroutine */ int zgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	doublecomplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int zgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	doublecomplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int zgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, 
	doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *scale, integer *m, doublecomplex *v, 
	integer *ldv, integer *info);

/* Subroutine */ int zgebal_(char *job, integer *n, doublecomplex *a, integer 
	*lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);

/* Subroutine */ int zgebd2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq, 
	doublecomplex *taup, doublecomplex *work, integer *info);

/* Subroutine */ int zgebrd_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq, 
	doublecomplex *taup, doublecomplex *work, integer *lwork, integer *
	info);

/* Subroutine */ int zgecon_(char *norm, integer *n, doublecomplex *a, 
	integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zgeequ_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, integer *info);

/* Subroutine */ int zgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, 
	doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork,
	doublereal *rwork, logical *bwork, integer *info);

/* Subroutine */ int zgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, doublecomplex *a, integer *lda, integer *sdim, 
	doublecomplex *w, doublecomplex *vs, integer *ldvs, doublereal *
	rconde, doublereal *rcondv, doublecomplex *work, integer *lwork, 
	doublereal *rwork, logical *bwork, integer *info);

/* Subroutine */ int zgeev_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *w, 
	doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *scale, doublereal *abnrm, 
	doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgegs_(char *jobvsl, char *jobvsr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, 
	integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgegv_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgehd2_(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zgehrd_(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zgelq2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgelqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zgels_(char *trans, integer *m, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zgelsx_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, 
	doublereal *rwork, integer *info);

/* Subroutine */ int zgelsy_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgeql2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgeqlf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zgeqp3_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgeqpf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	doublereal *rwork, integer *info);

/* Subroutine */ int zgeqr2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgeqrf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zgerfs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	doublereal *rwork, integer *info);

/* Subroutine */ int zgerq2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgerqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zgesc2_(integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *rhs, integer *ipiv, integer *jpiv, doublereal *scale);

/* Subroutine */ int zgesv_(integer *n, integer *nrhs, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *
	info);

/* Subroutine */ int zgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zgetc2_(integer *n, doublecomplex *a, integer *lda, 
	integer *ipiv, integer *jpiv, integer *info);

/* Subroutine */ int zgetf2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int zgetrf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int zgetri_(integer *n, doublecomplex *a, integer *lda, 
	integer *ipiv, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zgetrs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int zggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublecomplex *v, integer *ldv, integer *info);

/* Subroutine */ int zggbal_(char *job, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *b, integer *ldb, integer *ilo, integer *ihi, 
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info);

/* Subroutine */ int zgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, integer *sdim, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer 
	*ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork, 
	logical *bwork, integer *info);

/* Subroutine */ int zggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, char *sense, integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha, 
	doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, 
	doublecomplex *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int zggev_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *alpha, doublecomplex *beta, 
	doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, 
	doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int zggglm_(integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *d__, doublecomplex *x, doublecomplex *y, doublecomplex 
	*work, integer *lwork, integer *info);

/* Subroutine */ int zgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *q, integer *ldq, doublecomplex *z__, 
	integer *ldz, integer *info);

/* Subroutine */ int zgglse_(integer *m, integer *n, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, doublecomplex *d__, doublecomplex *x, 
	doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zggqrf_(integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b,
	integer *ldb, doublecomplex *taub, doublecomplex *work, integer *
	lwork, integer *info);

/* Subroutine */ int zggrqf_(integer *m, integer *p, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b,
	integer *ldb, doublecomplex *taub, doublecomplex *work, integer *
	lwork, integer *info);

/* Subroutine */ int zggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v, 
	integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work, 
	doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublecomplex *a, integer *lda, doublecomplex 
	*b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k, 
	integer *l, doublecomplex *u, integer *ldu, doublecomplex *v, integer 
	*ldv, doublecomplex *q, integer *ldq, integer *iwork, doublereal *
	rwork, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgtcon_(char *norm, integer *n, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *du2, integer *
	ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work, 
	integer *info);

/* Subroutine */ int zgtrfs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zgtsv_(integer *n, integer *nrhs, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *b, integer *ldb,
	integer *info);

/* Subroutine */ int zgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zgttrf_(integer *n, doublecomplex *dl, doublecomplex *
	d__, doublecomplex *du, doublecomplex *du2, integer *ipiv, integer *
	info);

/* Subroutine */ int zgttrs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zgtts2_(integer *itrans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb);

/* Subroutine */ int zhbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zhbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zhbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublecomplex *q, 
	integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__,
	integer *ldz, doublecomplex *work, doublereal *rwork, integer *iwork,
	integer *ifail, integer *info);

/* Subroutine */ int zhbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublecomplex *x, integer *ldx, doublecomplex *work, 
	doublereal *rwork, integer *info);

/* Subroutine */ int zhbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zhbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, doublecomplex *ab, integer *ldab, 
	doublecomplex *bb, integer *ldbb, doublecomplex *q, integer *ldq, 
	doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *
	abstol, integer *m, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int zhbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *d__, doublereal *e, 
	doublecomplex *q, integer *ldq, doublecomplex *work, integer *info);

/* Subroutine */ int zhecon_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublecomplex *work, integer *info);

/* Subroutine */ int zheev_(char *jobz, char *uplo, integer *n, doublecomplex 
	*a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork, 
	doublereal *rwork, integer *info);

/* Subroutine */ int zheevd_(char *jobz, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int zheevr_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, integer *isuppz, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int zheevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *iwork, integer *ifail, integer *
	info);

/* Subroutine */ int zhegs2_(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zhegst_(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zhegv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
	integer *info);

/* Subroutine */ int zhegvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
	integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zhegvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__,
	integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork,
	integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int zherfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	doublereal *rwork, integer *info);

/* Subroutine */ int zhesv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zhesvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
	integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zhetf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int zhetrd_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, 
	doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zhetrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info);

/* Subroutine */ int zhetri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int zhetrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int zhgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
	ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	info);

/* Subroutine */ int zhpcon_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, integer *info);

/* Subroutine */ int zhpev_(char *jobz, char *uplo, integer *n, doublecomplex 
	*ap, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zhpevd_(char *jobz, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zhpevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, doublereal *
	rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int zhpgst_(integer *itype, char *uplo, integer *n, 
	doublecomplex *ap, doublecomplex *bp, integer *info);

/* Subroutine */ int zhpgv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex 
	*z__, integer *ldz, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zhpgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex 
	*z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *
	rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
	info);

/* Subroutine */ int zhpgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublecomplex *ap, doublecomplex *bp, doublereal *
	vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
	integer *m, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int zhprfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zhpsv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zhpsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zhptrd_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *d__, doublereal *e, doublecomplex *tau, integer *info);

/* Subroutine */ int zhptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, integer *info);

/* Subroutine */ int zhptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int zhptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zhsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, doublecomplex *h__, integer *ldh, doublecomplex *
	w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr,
	integer *mm, integer *m, doublecomplex *work, doublereal *rwork, 
	integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int zhseqr_(char *job, char *compz, integer *n, integer *ilo,
	integer *ihi, doublecomplex *h__, integer *ldh, doublecomplex *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zlabrd_(integer *m, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *d__, doublereal *e, 
	doublecomplex *tauq, doublecomplex *taup, doublecomplex *x, integer *
	ldx, doublecomplex *y, integer *ldy);

/* Subroutine */ int zlacgv_(integer *n, doublecomplex *x, integer *incx);

/* Subroutine */ int zlacon_(integer *n, doublecomplex *v, doublecomplex *x, 
	doublereal *est, integer *kase);

/* Subroutine */ int zlacp2_(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublecomplex *b, integer *ldb);

/* Subroutine */ int zlacpy_(char *uplo, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb);

/* Subroutine */ int zlacrm_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *b, integer *ldb, doublecomplex *c__, 
	integer *ldc, doublereal *rwork);

/* Subroutine */ int zlacrt_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublecomplex *c__, doublecomplex *
	s);

/* Subroutine */ int zlaed0_(integer *qsiz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *qstore, 
	integer *ldqs, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zlaed7_(integer *n, integer *cutpnt, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__, 
	doublecomplex *q, integer *ldq, doublereal *rho, integer *indxq, 
	doublereal *qstore, integer *qptr, integer *prmptr, integer *perm, 
	integer *givptr, integer *givcol, doublereal *givnum, doublecomplex *
	work, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zlaed8_(integer *k, integer *n, integer *qsiz, 
	doublecomplex *q, integer *ldq, doublereal *d__, doublereal *rho, 
	integer *cutpnt, doublereal *z__, doublereal *dlamda, doublecomplex *
	q2, integer *ldq2, doublereal *w, integer *indxp, integer *indx, 
	integer *indxq, integer *perm, integer *givptr, integer *givcol, 
	doublereal *givnum, integer *info);

/* Subroutine */ int zlaein_(logical *rightv, logical *noinit, integer *n, 
	doublecomplex *h__, integer *ldh, doublecomplex *w, doublecomplex *v, 
	doublecomplex *b, integer *ldb, doublereal *rwork, doublereal *eps3, 
	doublereal *smlnum, integer *info);

/* Subroutine */ int zlaesy_(doublecomplex *a, doublecomplex *b, 
	doublecomplex *c__, doublecomplex *rt1, doublecomplex *rt2, 
	doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1);

/* Subroutine */ int zlaev2_(doublecomplex *a, doublecomplex *b, 
	doublecomplex *c__, doublereal *rt1, doublereal *rt2, doublereal *cs1,
	doublecomplex *sn1);

/* Subroutine */ int zlags2_(logical *upper, doublereal *a1, doublecomplex *
	a2, doublereal *a3, doublereal *b1, doublecomplex *b2, doublereal *b3,
	doublereal *csu, doublecomplex *snu, doublereal *csv, doublecomplex *
	snv, doublereal *csq, doublecomplex *snq);

/* Subroutine */ int zlagtm_(char *trans, integer *n, integer *nrhs, 
	doublereal *alpha, doublecomplex *dl, doublecomplex *d__, 
	doublecomplex *du, doublecomplex *x, integer *ldx, doublereal *beta, 
	doublecomplex *b, integer *ldb);

/* Subroutine */ int zlahef_(char *uplo, integer *n, integer *nb, integer *kb,
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, 
	integer *ldw, integer *info);

/* Subroutine */ int zlahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, integer *info);

/* Subroutine */ int zlahrd_(integer *n, integer *k, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *t, 
	integer *ldt, doublecomplex *y, integer *ldy);

/* Subroutine */ int zlaic1_(integer *job, integer *j, doublecomplex *x, 
	doublereal *sest, doublecomplex *w, doublecomplex *gamma, doublereal *
	sestpr, doublecomplex *s, doublecomplex *c__);

/* Subroutine */ int zlals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, doublecomplex *b, integer *ldb, 
	doublecomplex *bx, integer *ldbx, integer *perm, integer *givptr, 
	integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum,
	doublereal *poles, doublereal *difl, doublereal *difr, doublereal *
	z__, integer *k, doublereal *c__, doublereal *s, doublereal *rwork, 
	integer *info);

/* Subroutine */ int zlalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, doublecomplex *b, integer *ldb, doublecomplex *bx, 
	integer *ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *
	k, doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	rwork, integer *iwork, integer *info);

/* Subroutine */ int zlapll_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *ssmin);

/* Subroutine */ int zlapmt_(logical *forwrd, integer *m, integer *n, 
	doublecomplex *x, integer *ldx, integer *k);

/* Subroutine */ int zlaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int zlaqge_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int zlaqhb_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, char *equed);

/* Subroutine */ int zlaqhe_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	char *equed);

/* Subroutine */ int zlaqhp_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int zlaqp2_(integer *m, integer *n, integer *offset, 
	doublecomplex *a, integer *lda, integer *jpvt, doublecomplex *tau, 
	doublereal *vn1, doublereal *vn2, doublecomplex *work);

/* Subroutine */ int zlaqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, doublecomplex *a, integer *lda, integer *jpvt, 
	doublecomplex *tau, doublereal *vn1, doublereal *vn2, doublecomplex *
	auxv, doublecomplex *f, integer *ldf);

/* Subroutine */ int zlaqsb_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, char *equed);

/* Subroutine */ int zlaqsp_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int zlaqsy_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	char *equed);

/* Subroutine */ int zlar1v_(integer *n, integer *b1, integer *bn, doublereal 
	*sigma, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	lld, doublereal *gersch, doublecomplex *z__, doublereal *ztz, 
	doublereal *mingma, integer *r__, integer *isuppz, doublereal *work);

/* Subroutine */ int zlar2v_(integer *n, doublecomplex *x, doublecomplex *y, 
	doublecomplex *z__, integer *incx, doublereal *c__, doublecomplex *s, 
	integer *incc);

/* Subroutine */ int zlarcm_(integer *m, integer *n, doublereal *a, integer *
	lda, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc,
	doublereal *rwork);

/* Subroutine */ int zlarf_(char *side, integer *m, integer *n, doublecomplex 
	*v, integer *incv, doublecomplex *tau, doublecomplex *c__, integer *
	ldc, doublecomplex *work);

/* Subroutine */ int zlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublecomplex *v, integer 
	*ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, integer *
	ldc, doublecomplex *work, integer *ldwork);

/* Subroutine */ int zlarfg_(integer *n, doublecomplex *alpha, doublecomplex *
	x, integer *incx, doublecomplex *tau);

/* Subroutine */ int zlarft_(char *direct, char *storev, integer *n, integer *
	k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
	t, integer *ldt);

/* Subroutine */ int zlarfx_(char *side, integer *m, integer *n, 
	doublecomplex *v, doublecomplex *tau, doublecomplex *c__, integer *
	ldc, doublecomplex *work);

/* Subroutine */ int zlargv_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c__, integer *incc);

/* Subroutine */ int zlarnv_(integer *idist, integer *iseed, integer *n, 
	doublecomplex *x);

/* Subroutine */ int zlarrv_(integer *n, doublereal *d__, doublereal *l, 
	integer *isplit, integer *m, doublereal *w, integer *iblock, 
	doublereal *gersch, doublereal *tol, doublecomplex *z__, integer *ldz,
	integer *isuppz, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int zlartg_(doublecomplex *f, doublecomplex *g, doublereal *
	cs, doublecomplex *sn, doublecomplex *r__);

/* Subroutine */ int zlartv_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c__, doublecomplex *s, 
	integer *incc);

/* Subroutine */ int zlarz_(char *side, integer *m, integer *n, integer *l, 
	doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
	c__, integer *ldc, doublecomplex *work);

/* Subroutine */ int zlarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublecomplex 
	*v, integer *ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, 
	integer *ldc, doublecomplex *work, integer *ldwork);

/* Subroutine */ int zlarzt_(char *direct, char *storev, integer *n, integer *
	k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
	t, integer *ldt);

/* Subroutine */ int zlascl_(char *type__, integer *kl, integer *ku, 
	doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int zlaset_(char *uplo, integer *m, integer *n, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, integer *
	lda);

/* Subroutine */ int zlasr_(char *side, char *pivot, char *direct, integer *m,
	integer *n, doublereal *c__, doublereal *s, doublecomplex *a, 
	integer *lda);

/* Subroutine */ int zlassq_(integer *n, doublecomplex *x, integer *incx, 
	doublereal *scale, doublereal *sumsq);

/* Subroutine */ int zlaswp_(integer *n, doublecomplex *a, integer *lda, 
	integer *k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int zlasyf_(char *uplo, integer *n, integer *nb, integer *kb,
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, 
	integer *ldw, integer *info);

/* Subroutine */ int zlatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublecomplex *ab, integer *ldab, 
	doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int zlatdf_(integer *ijob, integer *n, doublecomplex *z__, 
	integer *ldz, doublecomplex *rhs, doublereal *rdsum, doublereal *
	rdscal, integer *ipiv, integer *jpiv);

/* Subroutine */ int zlatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *ap, doublecomplex *x, doublereal *
	scale, doublereal *cnorm, integer *info);

/* Subroutine */ int zlatrd_(char *uplo, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *e, doublecomplex *tau, 
	doublecomplex *w, integer *ldw);

/* Subroutine */ int zlatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, 
	doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int zlatrz_(integer *m, integer *n, integer *l, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work);

/* Subroutine */ int zlatzm_(char *side, integer *m, integer *n, 
	doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
	c1, doublecomplex *c2, integer *ldc, doublecomplex *work);

/* Subroutine */ int zlauu2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zlauum_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zpbcon_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *anorm, doublereal *
	rcond, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zpbequ_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, integer *info);

/* Subroutine */ int zpbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *
	ldafb, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);

/* Subroutine */ int zpbstf_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int zpbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
	ldb, integer *info);

/* Subroutine */ int zpbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, 
	integer *ldafb, char *equed, doublereal *s, doublecomplex *b, integer 
	*ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *
	ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, 
	integer *info);

/* Subroutine */ int zpbtf2_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int zpbtrf_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int zpbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
	ldb, integer *info);

/* Subroutine */ int zpocon_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zpoequ_(integer *n, doublecomplex *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int zporfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);

/* Subroutine */ int zposv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, char *equed, doublereal *s, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zpotf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zpotrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zpotri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zpotrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zppcon_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info);

/* Subroutine */ int zppequ_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int zpprfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, doublecomplex *b, integer *ldb,
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zppsv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, char *equed, doublereal *
	s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zpptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *info);

/* Subroutine */ int zpptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *info);

/* Subroutine */ int zpptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zptcon_(integer *n, doublereal *d__, doublecomplex *e, 
	doublereal *anorm, doublereal *rcond, doublereal *rwork, integer *
	info);

/* Subroutine */ int zptrfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);

/* Subroutine */ int zptsv_(integer *n, integer *nrhs, doublereal *d__, 
	doublecomplex *e, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zptsvx_(char *fact, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zpttrf_(integer *n, doublereal *d__, doublecomplex *e, 
	integer *info);

/* Subroutine */ int zpttrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zptts2_(integer *iuplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb);

/* Subroutine */ int zrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s);

/* Subroutine */ int zspcon_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, integer *info);

/* Subroutine */ int zspmv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *
	beta, doublecomplex *y, integer *incy);

/* Subroutine */ int zspr_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *ap);

/* Subroutine */ int zsprfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zspsv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zsptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, integer *info);

/* Subroutine */ int zsptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int zsptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zstedc_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int zstein_(integer *n, doublereal *d__, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublecomplex *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info);

/* Subroutine */ int zsteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int zsycon_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublecomplex *work, integer *info);

/* Subroutine */ int zsymv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublecomplex *beta, doublecomplex *y, integer *incy);

/* Subroutine */ int zsyr_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *a, integer *lda);

/* Subroutine */ int zsyrfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	doublereal *rwork, integer *info);

/* Subroutine */ int zsysv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zsysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
	integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zsytf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int zsytrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info);

/* Subroutine */ int zsytri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int zsytrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int ztbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublereal *rcond, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int ztbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);

/* Subroutine */ int ztbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int ztgevc_(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer 
	*ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
	ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork,
	integer *info);

/* Subroutine */ int ztgex2_(logical *wantq, logical *wantz, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz, 
	integer *j1, integer *info);

/* Subroutine */ int ztgexc_(logical *wantq, logical *wantz, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz, 
	integer *ifst, integer *ilst, integer *info);

/* Subroutine */ int ztgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
	ldz, integer *m, doublereal *pl, doublereal *pr, doublereal *dif, 
	doublecomplex *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int ztgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *tola, 
	doublereal *tolb, doublereal *alpha, doublereal *beta, doublecomplex *
	u, integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *q, 
	integer *ldq, doublecomplex *work, integer *ncycle, integer *info);

/* Subroutine */ int ztgsna_(char *job, char *howmny, logical *select, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer 
	*ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
	ldvr, doublereal *s, doublereal *dif, integer *mm, integer *m, 
	doublecomplex *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int ztgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, 
	doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, 
	doublereal *scale, doublereal *rdsum, doublereal *rdscal, integer *
	info);

/* Subroutine */ int ztgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, 
	doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, 
	doublereal *scale, doublereal *dif, doublecomplex *work, integer *
	lwork, integer *iwork, integer *info);

/* Subroutine */ int ztpcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *ap, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info);

/* Subroutine */ int ztprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int ztptri_(char *uplo, char *diag, integer *n, 
	doublecomplex *ap, integer *info);

/* Subroutine */ int ztptrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int ztrcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int ztrevc_(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer 
	*m, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int ztrexc_(char *compq, integer *n, doublecomplex *t, 
	integer *ldt, doublecomplex *q, integer *ldq, integer *ifst, integer *
	ilst, integer *info);

/* Subroutine */ int ztrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int ztrsen_(char *job, char *compq, logical *select, integer 
	*n, doublecomplex *t, integer *ldt, doublecomplex *q, integer *ldq, 
	doublecomplex *w, integer *m, doublereal *s, doublereal *sep, 
	doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int ztrsna_(char *job, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublereal *s, 
	doublereal *sep, integer *mm, integer *m, doublecomplex *work, 
	integer *ldwork, doublereal *rwork, integer *info);

/* Subroutine */ int ztrsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *c__, integer *ldc, doublereal *scale, 
	integer *info);

/* Subroutine */ int ztrti2_(char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int ztrtri_(char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int ztrtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int ztzrqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, integer *info);

/* Subroutine */ int ztzrzf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zung2l_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zung2r_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zungbr_(char *vect, integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zunghr_(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungl2_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zunglq_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungql_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungqr_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungr2_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zungrq_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungtr_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zunm2l_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int zunm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int zunmbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	lwork, integer *info);

/* Subroutine */ int zunmhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *a, integer *lda, 
	doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zunml2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int zunmlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zunmql_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zunmqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zunmr2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int zunmr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	info);

/* Subroutine */ int zunmrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zunmrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	lwork, integer *info);

/* Subroutine */ int zunmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	integer *info);

/* Subroutine */ int zupgtr_(char *uplo, integer *n, doublecomplex *ap, 
	doublecomplex *tau, doublecomplex *q, integer *ldq, doublecomplex *
	work, integer *info);

/* Subroutine */ int zupmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *c__,
	integer *ldc, doublecomplex *work, integer *info);

/*
The following prototypes are not present in the reference clapack.h distributed via netlib.
Nevertheless, the supporting code has always been present in the library.
*/

/* Subroutine */ int cgelsd_(integer *m, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *b, integer *ldb, f2c_real *s, f2c_real *rcond, 
	integer *rank, f2c_complex *work, integer *lwork, f2c_real *rwork, integer *
	iwork, integer *info);

/* Subroutine */ int cgelss_(integer *m, integer *n, integer *nrhs, f2c_complex *
	a, integer *lda, f2c_complex *b, integer *ldb, f2c_real *s, f2c_real *rcond, 
	integer *rank, f2c_complex *work, integer *lwork, f2c_real *rwork, integer *
	info);

/* Subroutine */ int cgesdd_(char *jobz, integer *m, integer *n, f2c_complex *a, 
	integer *lda, f2c_real *s, f2c_complex *u, integer *ldu, f2c_complex *vt, integer 
	*ldvt, f2c_complex *work, integer *lwork, f2c_real *rwork, integer *iwork, 
	integer *info);

/* Subroutine */ int cgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	f2c_complex *a, integer *lda, f2c_real *s, f2c_complex *u, integer *ldu, f2c_complex *
	vt, integer *ldvt, f2c_complex *work, integer *lwork, f2c_real *rwork, 
	integer *info);

/* Subroutine */ int chbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, f2c_complex *ab, integer *ldab, f2c_complex *bb, integer *ldbb, 
	f2c_real *w, f2c_complex *z__, integer *ldz, f2c_complex *work, integer *lwork, 
	f2c_real *rwork, integer *lrwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int chetd2_(char *uplo, integer *n, f2c_complex *a, integer *lda,
	f2c_real *d__, f2c_real *e, f2c_complex *tau, integer *info);

/* Complex */ void cladiv_(f2c_complex * ret_val, f2c_complex *x, f2c_complex *y);

/* Subroutine */ int clalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, f2c_real *d__, f2c_real *e, f2c_complex *b, integer *ldb, f2c_real *rcond, 
	integer *rank, f2c_complex *work, f2c_real *rwork, integer *iwork, 
	integer *info);

doublereal clangb_(char *norm, integer *n, integer *kl, integer *ku, 
	f2c_complex *ab, integer *ldab, f2c_real *work);

doublereal clange_(char *norm, integer *m, integer *n, f2c_complex *a, 
	integer *lda, f2c_real *work);

doublereal clangt_(char *norm, integer *n, f2c_complex *dl, f2c_complex *d__, 
	f2c_complex *du);

doublereal clanhb_(char *norm, char *uplo, integer *n, integer *k, 
	f2c_complex *ab, integer *ldab, f2c_real *work);

doublereal clanhe_(char *norm, char *uplo, integer *n, f2c_complex *a, 
	integer *lda, f2c_real *work);

doublereal clanhp_(char *norm, char *uplo, integer *n, f2c_complex *ap, 
	f2c_real *work);

doublereal clanhs_(char *norm, integer *n, f2c_complex *a, integer *lda, 
	f2c_real *work);

doublereal clanht_(char *norm, integer *n, f2c_real *d__, f2c_complex *e);

doublereal clansb_(char *norm, char *uplo, integer *n, integer *k, 
	f2c_complex *ab, integer *ldab, f2c_real *work);

doublereal clansp_(char *norm, char *uplo, integer *n, f2c_complex *ap, 
	f2c_real *work);

doublereal clansy_(char *norm, char *uplo, integer *n, f2c_complex *a, 
	integer *lda, f2c_real *work);

doublereal clantb_(char *norm, char *uplo, char *diag, integer *n, 
	integer *k, f2c_complex *ab, integer *ldab, f2c_real *work);

doublereal clantp_(char *norm, char *uplo, char *diag, integer *n, 
	f2c_complex *ap, f2c_real *work);

doublereal clantr_(char *norm, char *uplo, char *diag, integer *m, 
	integer *n, f2c_complex *a, integer *lda, f2c_real *work);

/* Subroutine */ int cpteqr_(char *compz, integer *n, f2c_real *d__, f2c_real *e, 
	f2c_complex *z__, integer *ldz, f2c_real *work, integer *info);

/* Subroutine */ int cstegr_(char *jobz, char *range, integer *n, f2c_real *d__, 
	f2c_real *e, f2c_real *vl, f2c_real *vu, integer *il, integer *iu, f2c_real *abstol, 
	integer *m, f2c_real *w, f2c_complex *z__, integer *ldz, integer *isuppz, 
	f2c_real *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);

doublereal dlamc3_(doublereal *a, doublereal *b);

doublereal dlamch_(char *cmach);

doublereal dlangb_(char *norm, integer *n, integer *kl, integer *ku, 
	doublereal *ab, integer *ldab, doublereal *work);

doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, 
	integer *lda, doublereal *work);

doublereal dlangt_(char *norm, integer *n, doublereal *dl, doublereal *d__, 
	doublereal *du);

doublereal dlanhs_(char *norm, integer *n, doublereal *a, integer *lda, 
	doublereal *work);

doublereal dlansb_(char *norm, char *uplo, integer *n, integer *k, 
	doublereal *ab, integer *ldab, doublereal *work);

doublereal dlansp_(char *norm, char *uplo, integer *n, doublereal *ap, 
	doublereal *work);

doublereal dlanst_(char *norm, integer *n, doublereal *d__, doublereal *e);

doublereal dlansy_(char *norm, char *uplo, integer *n, doublereal *a, 
	integer *lda, doublereal *work);

doublereal dlantb_(char *norm, char *uplo, char *diag, integer *n, integer *k,
	doublereal *ab, integer *ldab, doublereal *work);

doublereal dlantp_(char *norm, char *uplo, char *diag, integer *n, 
	doublereal *ap, doublereal *work);

doublereal dlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n,
	doublereal *a, integer *lda, doublereal *work);

doublereal dlapy2_(doublereal *x, doublereal *y);

doublereal dlapy3_(doublereal *x, doublereal *y, doublereal *z__);

doublereal dsecnd_();

doublereal dzsum1_(integer *n, doublecomplex *cx, integer *incx);

logical lsame_(char *ca, char *cb);

logical lsamen_(integer *n, char *ca, char *cb);

doublereal scsum1_(integer *n, f2c_complex *cx, integer *incx);

doublereal second_();

doublereal slamc3_(f2c_real *a, f2c_real *b);

doublereal slamch_(char *cmach);

doublereal slangb_(char *norm, integer *n, integer *kl, integer *ku, f2c_real *ab,
	integer *ldab, f2c_real *work);

doublereal slange_(char *norm, integer *m, integer *n, f2c_real *a, integer *lda, 
	f2c_real *work);

doublereal slangt_(char *norm, integer *n, f2c_real *dl, f2c_real *d__, f2c_real *du);

doublereal slanhs_(char *norm, integer *n, f2c_real *a, integer *lda, f2c_real *work);

doublereal slansb_(char *norm, char *uplo, integer *n, integer *k, f2c_real *ab, 
	integer *ldab, f2c_real *work);

doublereal slansp_(char *norm, char *uplo, integer *n, f2c_real *ap, f2c_real *work);

doublereal slanst_(char *norm, integer *n, f2c_real *d__, f2c_real *e);

doublereal slansy_(char *norm, char *uplo, integer *n, f2c_real *a, integer *lda, 
	f2c_real *work);

doublereal slantb_(char *norm, char *uplo, char *diag, integer *n, integer *k,
	f2c_real *ab, integer *ldab, f2c_real *work);

doublereal slantp_(char *norm, char *uplo, char *diag, integer *n, f2c_real *ap, 
	f2c_real *work);

doublereal slantr_(char *norm, char *uplo, char *diag, integer *m, integer *n,
	f2c_real *a, integer *lda, f2c_real *work);

doublereal slapy2_(f2c_real *x, f2c_real *y);

doublereal slapy3_(f2c_real *x, f2c_real *y, f2c_real *z__);

/* Subroutine */ int zgelsd_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *s, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zgesdd_(char *jobz, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, 
	integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, 
	integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zhbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zhetd2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, 
	integer *info);

/* Double Complex */ void zladiv_(doublecomplex * ret_val, doublecomplex *x, 
	doublecomplex *y);

/* Subroutine */ int zlalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, doublereal *d__, doublereal *e, doublecomplex *b, integer *ldb,
	doublereal *rcond, integer *rank, doublecomplex *work, 
	doublereal *rwork, integer *iwork, integer *info);

doublereal zlangb_(char *norm, integer *n, integer *kl, integer *ku, 
	doublecomplex *ab, integer *ldab, doublereal *work);

doublereal zlange_(char *norm, integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work);

doublereal zlangt_(char *norm, integer *n, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du);

doublereal zlanhb_(char *norm, char *uplo, integer *n, integer *k, 
	doublecomplex *ab, integer *ldab, doublereal *work);

doublereal zlanhe_(char *norm, char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work);

doublereal zlanhp_(char *norm, char *uplo, integer *n, doublecomplex *ap, 
	doublereal *work);

doublereal zlanhs_(char *norm, integer *n, doublecomplex *a, integer *lda, 
	doublereal *work);

doublereal zlanht_(char *norm, integer *n, doublereal *d__, doublecomplex *e);

doublereal zlansb_(char *norm, char *uplo, integer *n, integer *k, 
	doublecomplex *ab, integer *ldab, doublereal *work);

doublereal zlansp_(char *norm, char *uplo, integer *n, doublecomplex *ap, 
	doublereal *work);

doublereal zlansy_(char *norm, char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work);

doublereal zlantb_(char *norm, char *uplo, char *diag, integer *n, integer *k,
	doublecomplex *ab, integer *ldab, doublereal *work);

doublereal zlantp_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *ap, doublereal *work);

doublereal zlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n,
	doublecomplex *a, integer *lda, doublereal *work);

/* Subroutine */ int zpteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int zstegr_(char *jobz, char *range, integer *n, 
	doublereal *d__, doublereal *e, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublecomplex *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);


#ifdef __cplusplus
}
#endif
#endif /* __CLAPACK_H */
