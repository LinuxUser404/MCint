#ifndef _SOBSEQ_H
#define _SOBSEQ_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sobseq_err(msg) {printf("sobseq_err: %s.\n",msg); exit(255);}
#define MAXDIM 6
#define MAXBIT 62

const double sobseq_boxlength = 1.;

inline
void sobseq(double *x, const int n, const int set)
/* NB: x runs from 0..n-1, NOT NR's usual 1..n */
/* n is dim of system, set => initialize if == 0 */
{
    int j,k,l;
    unsigned long long i,im,ipp;
    static double fac;
    static unsigned long long in,ix[MAXDIM],*iu[MAXBIT+1];
    static unsigned long long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
    static unsigned long long ip[MAXDIM+1]={0,0,1,1,2,1,4};
    static unsigned long long iv[MAXDIM*MAXBIT+1],ivb[MAXDIM*MAXBIT+1]={
		0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

    if (!set) {
        if (n < 1 || n > MAXDIM) 
            sobseq_err("dim parameter is out of range");
        for (k = 0; k <= MAXDIM*MAXBIT; k++)
            iv[k] = ivb[k];
	for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM)
            iu[j] = &iv[k];
	for (k=1;k<=MAXDIM;k++) {
            ix[k-1] = 0;
	    for (j=1;j<=(int)mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
	    for (j=mdeg[k]+1;j<=MAXBIT;j++) {
		ipp=ip[k];
		i=iu[j-mdeg[k]][k];
		i ^= (i >> mdeg[k]);
		for (l=mdeg[k]-1;l>=1;l--) {
		    if (ipp & 1) i ^= iu[j-l][k];
			ipp >>= 1;
		}
		iu[j][k]=i;
	    }
	}
	fac = sobseq_boxlength/(((unsigned long long)1) << MAXBIT);
	in = 0;
    }
    im=in;
    for (j = 0;j < MAXBIT;j++) {
        if (!(im & (unsigned long long)1)) break;
	im >>= 1;
    }

#if SOBSEQ_ERRCHK
    if (j == MAXBIT)
        sobseq_err("MAXBIT too small");
#endif

    im = j*MAXDIM + 1;
    for (k = 0; k < n; k++) {
        ix[k] ^= iv[im+k];
        x[k]=ix[k]*fac;
    }
    in++;
}


#undef MAXDIM

#if 1 
#define NONZEROBIT(c,v) {\
     for (c = 0; c < MAXBIT; c++) { \
         if (!((unsigned long long)v & (unsigned long long)1)) break; \
	 v >>= 1; \
     } \
     if (c == MAXBIT) sobseq_err("MAXBIT too small"); \
}
#else
#define NONZEROBIT(c,v) {\
     for (c = 0 ; ; c++) { \
         if (!(v & 1)) break; \
	 v >>= 1; \
     } \
}
#endif

#undef MAXDIM
#define MAXDIM 12

inline
void sobseqb(double *x, const int n, const int nstep, const int offset, 
	     const int set)
/* sobseq with offsets and steps along the sequence. */
/* to init, set = 0. Note: on init, the first value of x is packed. */
/* Bad Limitations: nstep must be a power of two, offset is [0..nstep-1] */
/* NB: x runs from 0..n-1, NOT NR's usual 1..n */
{
    int j,k,l,q;
    unsigned long long i,im,ipp;
    static double fac;
    static unsigned long long in,ix[MAXDIM];
    static unsigned long long *iu[MAXBIT+1],*iustep[MAXBIT+1];
    static unsigned long long mdeg[MAXDIM+1]={0,1,2,3,3,4,4,5,5,5,5,5,5};
    static unsigned long long ip[MAXDIM+1]=  {0,0,1,1,2,1,4,2,13,7,14,11,4};
    static unsigned long long iv[MAXDIM*MAXBIT+1],ivstep[MAXDIM*MAXBIT+1];
    static unsigned long long ivb[MAXDIM*MAXBIT+1]={ 0,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          3, 1, 3, 1, 3, 1, 3, 3, 1, 3, 1, 3,
          5, 7, 7, 5, 1, 3, 3, 7, 5, 1, 1, 5,
         15,11, 5, 3, 1, 7, 9,13,11, 1, 3, 7,
         17,13, 7,15, 9,31, 9, 3,27,15,29,21   
                                              };
    if (!set) {

        if (n < 1 || n > MAXDIM) 
            sobseq_err("dim parameter is out of range");
        for (k = 0; k <= MAXDIM*MAXBIT; k++)
            iv[k] = ivb[k];
	for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM)
            iu[j] = &iv[k];
	for (k=1;k<=MAXDIM;k++) {
            ix[k-1] = 0; /* not NR's scheme */
	    for (j=1;j<=(int)mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
	    for (j=mdeg[k]+1;j<=MAXBIT;j++) {
		ipp=ip[k];
		i=iu[j-mdeg[k]][k];
		i ^= (i >> mdeg[k]);
		for (l=mdeg[k]-1;l>=1;l--) {
		    if (ipp & 1) i ^= iu[j-l][k];
			ipp >>= 1;
		}
		iu[j][k]=i;
	    }
	}

	fac = sobseq_boxlength/((double)(((unsigned long long)1) << MAXBIT));
	in = 0;
        /* the NR stuff is now complete (but note new indexing on ix[]) */

        /* you might use this like the regular sobseq. just set nstep = 1... */
        if (nstep == 1) { 
            for (k = 0; k <= MAXDIM*MAXBIT; k++)
                ivstep[k] = iv[k];
            sobseqb(x,n,nstep,offset,1);
            return;
	}

        /* next we set up for skipping over seqeuence elements */
	for (j = 1, k = 0; j <= MAXBIT; j++, k += MAXDIM)
            iustep[j] = &ivstep[k];
 	for (k = 1; k <= n; k++) 
            for (j = 1; j <= MAXBIT; j++)
                iustep[j][k] = 0;


        /* j will be log2(nstep) */
        for (j = 1; (1 << j) < nstep; j++); 
	for ( ; j < MAXBIT; j++) {
            for (l = 0; l < nstep; l++) {
                im = ((unsigned long long)1 << j) - nstep + l;
                NONZEROBIT(q,im)
                q++;
 	        for (k = 1; k <= n; k++)
                    iustep[j+1][k] ^= iu[q][k];
	    }
	}

        /* now set up the offset...this needs a shitload of work */
        /* on the other hand its just init code... */
        if (offset == nstep - 1) {
            in = nstep - 1;
            sobseqb(x,n,nstep,offset,1);
	} else {
            for (j = 0, in = 0; j <= offset; j++, in++) {
                im=in;
                NONZEROBIT(k,im)
	        im = k*MAXDIM + 1;
                for (k = 0; k < n; k++)
                    ix[k] ^= iv[im+k];
	        if (j == offset)
                    for (k = 0; k < n; k++)
                        x[k] = ix[k]*fac;
	    }
            in = nstep - 1;
        }
    } else {
        im=in;
        NONZEROBIT(j,im)
	im = j*MAXDIM + 1;
	for (k = 0; k < n; k++) {
            ix[k] ^= ivstep[im+k];
            x[k]=ix[k]*fac;
	}
        in += nstep;
    }
}

static unsigned long long in_sbq,ix_sbq[MAXDIM],ivstep_sbq[MAXDIM*MAXBIT+1];
static int ndim_sbq,nstep_sbq;
static double fac;

inline
void sobseqbq_init(const int n, const int ns, const int offset)
/* sobseq with offsets and steps along the sequence. */
/* to init, set = 0. Note: on init, the first value of x is packed. */
/* Bad Limitations: nstep must be a power of two, offset is [0..nstep-1] */
/* NB: x runs from 0..n-1, NOT NR's usual 1..n */
{
    int j,k,l,q;
    unsigned long long i,im,ipp;
    static unsigned long long *iu[MAXBIT+1],*iustep[MAXBIT+1];
    static unsigned long long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
    static unsigned long long ip[MAXDIM+1]={0,0,1,1,2,1,4};
    static unsigned long long iv[MAXDIM*MAXBIT+1];
    static unsigned long long ivb[MAXDIM*MAXBIT+1]={
		0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

    ndim_sbq = n;
    nstep_sbq = ns;

    if (n < 1 || n > MAXDIM) 
        sobseq_err("dim parameter is out of range");
    for (k = 0; k <= MAXDIM*MAXBIT; k++)
        iv[k] = ivb[k];
    for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM)
        iu[j] = &iv[k];
    for (k=1;k<=MAXDIM;k++) {
        ix_sbq[k-1] = 0; /* not NR's scheme */
        for (j=1;j<=(int)mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
        for (j=mdeg[k]+1;j<=MAXBIT;j++) {
	    ipp=ip[k];
	    i=iu[j-mdeg[k]][k];
	    i ^= (i >> mdeg[k]);
	    for (l=mdeg[k]-1;l>=1;l--) {
	        if (ipp & 1) i ^= iu[j-l][k];
	            ipp >>= 1;
	    }
	    iu[j][k]=i;
	}
    }
    fac = sobseq_boxlength/((double)(((unsigned long long)1) << MAXBIT));
    in_sbq = 0;
    /* the NR stuff is now complete (but note new indexing on ix_sbq[]) */

    /* you might use this like the regular sobseq. just set nstep = 1... */
    if (ns == 1) { 
        for (k = 0; k <= MAXDIM*MAXBIT; k++)
            ivstep_sbq[k] = iv[k];
        im = 0;
        NONZEROBIT(k,im)
        im = k*MAXDIM + 1;
        for (k = 0; k < n; k++)
            ix_sbq[k] ^= iv[im+k];
        in_sbq = 1;

        return;
    }

        /* next we set up for skipping over seqeuence elements */
    for (j = 1, k = 0; j <= MAXBIT; j++, k += MAXDIM)
        iustep[j] = &ivstep_sbq[k];
    for (k = 1; k <= n; k++) 
        for (j = 1; j <= MAXBIT; j++)
            iustep[j][k] = 0;

    /* j will be log2(ns) */
    for (j = 1; (1 << j) < ns; j++); 
    for ( ; j < MAXBIT; j++) {
        for (l = 0; l < ns; l++) {
            im = (1 << j) - ns + l;
            NONZEROBIT(q,im)
            q++;
            for (k = 1; k <= n; k++)
                iustep[j+1][k] ^= iu[q][k];
        }
    }

    /* now set up the offset...this needs a shitload of work */
    /* on the other hand its just init code... */
    for (j = 0, in_sbq = 0; j <= offset; j++, in_sbq++) {
        im = in_sbq;
        NONZEROBIT(k,im)
        im = k*MAXDIM + 1;
        for (k = 0; k < n; k++)
            ix_sbq[k] ^= iv[im+k];
    }
    in_sbq = (offset == ns - 1) ? 2*ns - 1 : ns - 1;
}

inline
void sobseqbq(double *x)
{
    unsigned long long im;
    unsigned int j,k,l;

    im = in_sbq;
    NONZEROBIT(j,im)
    l = j*MAXDIM + 1;
    for (k = 0; k < (unsigned int)ndim_sbq; k++) {
        x[k]=ix_sbq[k]*fac;
        ix_sbq[k] ^= ivstep_sbq[l+k];
    }
    in_sbq += nstep_sbq;
}



#define MAXDV 10

static unsigned long long in_sbx,it_sbx;
static unsigned long long ix_sbx[MAXDIM],ivstep_sbx[MAXDV][MAXDIM*MAXBIT+1];
static int ndim_sbx,ndv_sbx,dvstep_sbx[MAXDV];
static double fac_sbx;

inline
void sobseqbx_init(int n, int ns, int offset)
/* sobseq with offsets and steps along the sequence. */
/* to init, set = 0. Note: on init, the first value of x is packed. */
/* Bad Limitations: nstep must be a power of two, offset is [0..nstep-1] */
/* NB: x runs from 0..n-1, NOT NR's usual 1..n */
{
    int j,k,l,q,dv;
    unsigned long long i,im,ipp;
    static unsigned long long *iu[MAXBIT+1],*iustep[MAXBIT+1];
    static unsigned long long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
    static unsigned long long ip[MAXDIM+1]={0,0,1,1,2,1,4};
    static unsigned long long iv[MAXDIM*MAXBIT+1];
    static unsigned long long ivb[MAXDIM*MAXBIT+1]={
		0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

    ndim_sbx = n;
    it_sbx = offset;

    if (n < 1 || n > MAXDIM) 
        sobseq_err("dim parameter is out of range");
    for (k = 0; k <= MAXDIM*MAXBIT; k++)
        iv[k] = ivb[k];
    for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM)
        iu[j] = &iv[k];
    for (k=1;k<=MAXDIM;k++) {
        ix_sbx[k-1] = 0; /* not NR's scheme */
        for (j=1;j<=(int)mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
        for (j=mdeg[k]+1;j<=MAXBIT;j++) {
	    ipp=ip[k];
	    i=iu[j-mdeg[k]][k];
	    i ^= (i >> mdeg[k]);
	    for (l=mdeg[k]-1;l>=1;l--) {
	        if (ipp & 1) i ^= iu[j-l][k];
	            ipp >>= 1;
	    }
	    iu[j][k]=i;
	}
    }
    fac_sbx = sobseq_boxlength/((double)(((unsigned long long)1) << MAXBIT));
    in_sbx = 0;

    /* the dir-vec/NR stuff is complete (but note new indexing on ix_sbx[]) */

    /* next we set up for skipping over seqeuence elements */

    for (j = ns, k = 0; j; j -= l, k++) {
        for (l = 1; l < j; l <<= 1);
        if (l > j) l >>= 1;
        if (!l) break;
        dvstep_sbx[k] = l;
    }
    ndv_sbx = k;
 
    for (dv = 0; dv < ndv_sbx; dv++) {

        if (dvstep_sbx[dv] == 1) {
            for (k = 0; k <= MAXDIM*MAXBIT; k++)
                ivstep_sbx[dv][k] = iv[k];
            continue;
        }

        for (j = 1, k = 0; j <= MAXBIT; j++, k += MAXDIM)
            iustep[j] = &ivstep_sbx[dv][k];
        for (k = 1; k <= n; k++) 
            for (j = 1; j <= MAXBIT; j++)
                iustep[j][k] = 0;

        /* j will be log2(ns) */
        for (j = 1; (1 << j) < dvstep_sbx[dv]; j++); 
        for ( ; j < MAXBIT; j++) {
            for (l = 0; l < dvstep_sbx[dv]; l++) {
                im = (1 << j) - dvstep_sbx[dv] + l;
                NONZEROBIT(q,im)
                q++;
                for (k = 1; k <= n; k++)
                    iustep[j+1][k] ^= iu[q][k];
            }
        }
    }

    for (j = 0, in_sbx = 0; j <= offset; j++, in_sbx++) {
        im = in_sbx;
        NONZEROBIT(k,im)
        im = k*MAXDIM + 1;
        for (k = 0; k < n; k++)
            ix_sbx[k] ^= iv[im+k];
    }
    in_sbx = (offset >= dvstep_sbx[0]-1) ? 2*dvstep_sbx[0]-1 : dvstep_sbx[0]-1;
}

/* THIS IS BROKEN; it only works when nstep = 2^k + 2^k' */

inline
void sobseqbx(double *x)
{
    unsigned long long im,imcum;
    unsigned int j,k,l,dv;

    for (k = 0; k < (unsigned int)ndim_sbx; k++)
        x[k]=ix_sbx[k]*fac_sbx;

    imcum = in_sbx;
    for (dv = 0; dv < (unsigned int)ndv_sbx; dv++) {

        it_sbx +=dvstep_sbx[dv];
        if (dvstep_sbx[dv] == 1) im = it_sbx;
        else {
            im = imcum;
            if (dv != (unsigned int)ndv_sbx-1) 
               while (imcum < it_sbx) 
                    imcum += dvstep_sbx[dv+1];
	}
#if 0
printf(" a: %d %d",(int)im,(int)it_sbx);
#endif
        NONZEROBIT(j,im)
        l = j*MAXDIM + 1;
        for (k = 0; k < (unsigned int)ndim_sbx; k++) 
            ix_sbx[k] ^= ivstep_sbx[dv][l+k];

#if 0
printf(" b: %g ",ix_sbx[0]*fac_sbx);
#endif
    }

    while (in_sbx <= it_sbx) 
        in_sbx += (unsigned long long)dvstep_sbx[0];
}


#undef MAXBIT
#undef MAXDIM
#undef NONZEROBIT
#undef sobseq_err
#endif
