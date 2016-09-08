#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>
#include <algorithm>

using namespace std;

inline void sdrand(unsigned long int seed);
inline double drand();
double integrand_discon(valarray<double>& x);
double integrand_smooth(valarray<double>& x);
double drand(const double x) { return drand(); }

void setsobseq(valarray<double>& x);
void setpseudo(valarray<double>& x);

int main(int argc, char *argv[])
{

  size_t dim = 12; // dimension of domain of integration
  long int n = 10000; // number of monte carlo sample points
  long int seed = 987654321; // rndgen seed. hummmm
  bool smoo = 1;    // smooth of discontinuous integrand?

  // intro stuff: 

  if (argc > 1) {
    while (argc>1) {
      string arg = *++argv;
      argc-=2;
      if (argc<1) { 
	cerr << " error reading cmd line at "<<arg<<endl;
	exit(255);
      }
      if (arg=="-dim") dim = size_t(atol(*++argv));
      else if (arg=="-n") n = atol(*++argv);
      else if (arg=="-seed") seed = atol(*++argv);
      else if (arg=="-smooth") smoo = atoi(*++argv);
      else {
	cerr << " unrecognized cmd-line option " << arg << endl;
	exit(0);
      }
    }
  } else {
    cout << "        *** monte carlo integrator *** " << endl;
    cout << " " << endl;
    cout << " this code integrates over a unit hypercube in a " << endl;
    cout << " in a DIM-dimensional space. by default the integrand" << endl;
    cout << " is a smooth function. it is sampled using both quasi-" << endl;
    cout << " random (sobol) and pseudorandom schemes. " << endl;
    cout << " " << endl;
    cout << " Parameters can be modified with command line opts: " << endl;
    cout << " " << endl;
    cout << "   -dim DIM      # sets the number of dimensions" << endl;
    cout << "   -n NPTS       # number os sample points " << endl;
    cout << "   -smooth [1|0] # [smooth|discontinuous] integrand. "<<endl;
    cout << "   -seed SEED    # seed for psuedorandom number generator."<<endl;
    cout << " " << endl;
    cout << " Defaults: " << endl;
    cout << " " << endl;
    cout << "       "<<argv[0]<<" -dim "<<dim<<" -n "<<n;
    cout << " -smooth "<<((smoo)?1:0);
    cout <<" -seed "<<seed; 
    cout << endl << endl;
  }

  //
  sdrand(seed);

  double (*integrand)(valarray<double>& x);
  integrand = (smoo) ? integrand_smooth : integrand_discon;


  double Iq   = 0.0;             // the integral, quasirandom method
  double varq = 0.0;             // variance in the answer

  double Ip   = 0.0;             // the integral, pseudo randaom method
  double varp = 0.0;             // variance in the answer

  valarray<double> x(dim);      // a sample point
   
  for (int i=0; i < n; i++) {

    // quasirandom:
    setsobseq(x);
    double Isamp = integrand(x);
    Iq += Isamp;
    varq += Isamp*Isamp;

    // pseudorandom:
    setpseudo(x);
    Isamp = integrand(x);
    Ip += Isamp;
    varp += Isamp*Isamp;

  }

  // cleanup:
  Iq /= n;
  varq = varq/n-Iq*Iq;
  double Iuncertq = sqrt(varq/n);

  Ip /= n;
  varp = varp/n-Ip*Ip;
  double Iuncertp = sqrt(varp/n);

  double Iexact;
  if (smoo) { 
    Iexact = 1.0;
    printf(" Integral is integral_0^1 dx^dim exp(-x^2/2) * norm_factor\n");
    printf(" Exact value = 1.0\n");
  } else {
    Iexact = (pow(.5,dim)*pow(acos(-1.),0.5*dim)/exp(lgamma(1.+0.5*dim))); 
    printf(" Integral is Vol of hypersphere (R=1) * norm_factor.\n");
    printf(" Exact value = %1.15g\n",Iexact);
  }

  printf(" Quasirandom method result = %1.15g +/- %1.15g\n",Iq,Iuncertq);
  printf(" Abs err = %1.15g\n",fabs(Iq-Iexact));
  printf(" Psuedorandom method result = %1.15g +/- %1.15g\n",Ip,Iuncertp);
  printf(" Abs err = %1.15g\n",fabs(Ip-Iexact));
  
  return 0;
}

double integrand_discon(valarray<double>& x)
  // WARNING! THIS DESTROYS THE ARRAY!
{
  x *= x;
  double s = sqrt(x.sum());
  return (s<=1.0) ? 1.0 : 0.0;  
}

double integrand_smooth(valarray<double>& x)
  // WARNING! THIS DESTROYS THE ARRAY!
{
  const double norm = 1.1687371345136332892;
  static size_t dimlast;
  static double normfac;
  if (x.size()!=dimlast){ dimlast=x.size(); normfac=pow(norm,x.size()); }
  x *= x;
  double x2 = x.sum();
  return normfac*exp(-0.5*x2);
}


// below is junk for generating sample points

static int highbit;

double drandmask(double x) {
  double xs = drand();
  xs *= (1<<highbit);
  xs = floor(xs);
  xs /= (1<<highbit);
  return xs;
}

void setpseudo(valarray<double>& x) 
{
  static long long ct;
  ++ct;
  if (!((1<<(highbit))&(ct/x.size()))) {
    highbit++;
  }
  x = x.apply(drandmask);

}

bool cmp_frac(const double& a, const double& b) { // increasing status 
  return (a-int(a)<b-int(b));
}


void sobseq(double *, int, int);
void sobseqb(double *, unsigned int, unsigned int, 
	     unsigned int, unsigned int);

void setsobseq(valarray<double>& x) 
{
  // this is broken!!!!
 static int init;
  size_t d = x.size();
  size_t ct = 0;
  const int dimperblock = 12;
  int blocks = d/dimperblock;
  double xb[dimperblock];
  
  sobseqb(xb,dimperblock,1,0,init++);
  for (int b=0; b<blocks; b++) {
    // for each b get a unique map of xb[j]->x[i]
    for (int j=0; j<dimperblock; j++) x[ct++]=xb[j];
  }
  int xtra=d%dimperblock;
  if (xtra) {
    for (int j=0; j<xtra; j++) x[ct++]=xb[j];
  }
  if (0) { 
    for (size_t j=0; j<x.size(); j++) cout << x[j] << " ";
    cout << endl;
  }

}

// my random # stuff

inline 
void sdrandMT(unsigned long int seed);
inline 
double drandMT();

inline 
void sdrand(unsigned long int seed) { sdrandMT(seed); }
inline 
double drand() { return drandMT(); }

typedef unsigned long uint32;
inline void seedMT(uint32 seed);
inline uint32 randomMT(void);

inline void  sdrandMT(uint32 seed)
{
  seedMT(seed);
}

#define DBL_MANTIS_BITS 52
#define DBL_EXPSGN_BITS 12


inline 
double drandMT() 
{
  unsigned long long ll;
  const double EPS = (1.0/(1ULL << DBL_MANTIS_BITS));

  ll = randomMT();
  ll <<= 32; 
  ll += randomMT();
  return (ll >> DBL_EXPSGN_BITS) * (EPS);
}

#undef DBL_MANTIS_BITS
#undef DBL_EXPSGN_BITS

// Mersenne Twister : Matsumoto & Nishimura (1996/1997)
// Code is Shawn Cokus's optimized version.
// see: http://www.math.keio.ac.jp/~matumoto/emt.html
// for more info....

#define N              (624)                 // length of state vector
#define M              (397)                 // a period parameter
#define K              (0x9908B0DFU)         // a magic constant
#define hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // move hi bit of u to hi bit of v

static uint32   MT_state[N+1];     // MT_state vector + 1 extra to not violate ANSI C
static uint32   *MT_next;          // next random value is computed from here
static int      MT_left = -1;      // can *next++ this many times before reloading


inline void seedMT(uint32 seed)
{
  register uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = MT_state;
  register int    j;

  for(MT_left=0, *s++=x, j=N; --j;
      *s++ = (x*=69069U) & 0xFFFFFFFFU);
}


inline uint32 reloadMT(void)
{
  register uint32 *p0=MT_state, *p2=MT_state+2, *pM=MT_state+M, s0, s1;
  register int    j;

  if(MT_left < -1)
    seedMT(4357U);

  MT_left=N-1, MT_next=MT_state+1;

  for(s0=MT_state[0], s1=MT_state[1], j=N-M+1; --j; s0=s1, s1=*p2++)
    *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

  for(pM=MT_state, j=M; --j; s0=s1, s1=*p2++)
    *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

  s1=MT_state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
  s1 ^= (s1 >> 11);
  s1 ^= (s1 <<  7) & 0x9D2C5680U;
  s1 ^= (s1 << 15) & 0xEFC60000U;
  return(s1 ^ (s1 >> 18));
}


inline uint32 randomMT(void)  // bcb: nixed the inline...fix this.
{
  uint32 y;

  if(--MT_left < 0)
    return(reloadMT());

  y  = *MT_next++;
  y ^= (y >> 11);
  y ^= (y <<  7) & 0x9D2C5680U;
  y ^= (y << 15) & 0xEFC60000U;
  return(y ^ (y >> 18));
}


#undef N
#undef M
#undef K
#undef hiBit
#undef loBit
#undef loBits
#undef mixBits



/*
    a sobol sequence generator, and one that can be strided for parallelism.
    set for 6-d max for now.

    bcb 1996--1999
    bcb 2000 convert to c++ compatability (whoops...)
 */

#define SOBSEQ_MAXBIT 62
#define SOBSEQ_ERRCHK 0
void sobseqb(double *, unsigned int, unsigned int, 
	     unsigned int, unsigned int);
void sobseq(double *, int, int);


#define sobseq_err(msg) {fprintf(stderr,"sobseq_err: %s.\n",msg); exit(255);}
#define MAXDIM 6
#define MAXBIT SOBSEQ_MAXBIT


double sobseq_boxlength = 1.;

void sobseq(double *x, int n, int set)
/* NB: x runs from 0..n-1, NOT NR's usual 1..n */
/* n is dim of system, set => initialize if == 0 */
{
    int k;
    unsigned long long i,j,l,im,ipp;
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
	    for (j=1;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
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

#define MAXDIM 12

void sobseqb(double *x, unsigned int n, unsigned int nstep, 
	     unsigned int offset, unsigned int set)
/* sobseq with offsets and steps along the sequence. */
/* to init, set = 0. Note: on init, the first value of x is packed. */
/* Bad Limitations: nstep must be a power of two, offset is [0..nstep-1] */
/* NB: x runs from 0..n-1, NOT NR's usual 1..n */
{
    unsigned int k,l,q;
    unsigned long long i,j,im,ipp;
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
	    for (j=1;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
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
        for (j = 1U; (1U << j) < nstep; j++); 
	for ( ; j < MAXBIT; j++) {
            for (l = 0; l < nstep; l++) {
                im = ((unsigned long long)1 << j) - (unsigned long long)(nstep + l);
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
static unsigned int ndim_sbq,nstep_sbq;
static double fac;

void sobseqbq_init(unsigned int n, unsigned int ns, unsigned int offset)
/* sobseq with offsets and steps along the sequence. */
/* to init, set = 0. Note: on init, the first value of x is packed. */
/* Bad Limitations: nstep must be a power of two, offset is [0..nstep-1] */
/* NB: x runs from 0..n-1, NOT NR's usual 1..n */
{
    unsigned int j,k,l,q;
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
        for (j=1;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
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
    for (j = 1U; (1U << j) < ns; j++); 
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


void sobseqbq(double *x)
{
    unsigned int im,j,k,l;

    im = in_sbq;
    NONZEROBIT(j,im)
    l = j*MAXDIM + 1;
    for (k = 0U; k < ndim_sbq; k++) {
        x[k]=ix_sbq[k]*fac;
        ix_sbq[k] ^= ivstep_sbq[l+k];
    }
    in_sbq += nstep_sbq;
}



#define MAXDV 10

static unsigned long long in_sbx,it_sbx;
static unsigned long long ix_sbx[MAXDIM],ivstep_sbx[MAXDV][MAXDIM*MAXBIT+1];
static unsigned int ndim_sbx,nstep_sbx,ndv_sbx,dvstep_sbx[MAXDV];
static double fac_sbx;

void sobseqbx_init(unsigned int n, unsigned int ns, unsigned int offset)
/* sobseq with offsets and steps along the sequence. */
/* to init, set = 0. Note: on init, the first value of x is packed. */
/* Bad Limitations: nstep must be a power of two, offset is [0..nstep-1] */
/* NB: x runs from 0..n-1, NOT NR's usual 1..n */
{
    unsigned int j,k,l,q,dv;
    unsigned long long i,im,ipp;
    static unsigned long long *iu[MAXBIT+1],*iustep[MAXBIT+1];
    static unsigned int mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
    static unsigned long long ip[MAXDIM+1]={0,0,1,1,2,1,4};
    static unsigned long long iv[MAXDIM*MAXBIT+1];
    static unsigned long long ivb[MAXDIM*MAXBIT+1]={
		0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

    ndim_sbx = n;
    nstep_sbx = ns;
    it_sbx = offset;

    if (n < 1 || n > MAXDIM) 
        sobseq_err("dim parameter is out of range");
    for (k = 0; k <= MAXDIM*MAXBIT; k++)
        iv[k] = ivb[k];
    for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM)
        iu[j] = &iv[k];
    for (k=1;k<=MAXDIM;k++) {
        ix_sbx[k-1] = 0; /* not NR's scheme */
        for (j=1U;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
        for (j=mdeg[k]+1;j<=MAXBIT;j++) {
	    ipp=ip[k];
	    i=iu[j-mdeg[k]][k];
	    i ^= (i >> mdeg[k]);
	    for (l=mdeg[k]-(unsigned long long)1;l>=1U;l--) {
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
        for (l = (unsigned long long)1; l < (unsigned long long)j; l <<= 1);
        if (l > j) l >>= 1ULL;
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
        for (j = 1U; (1U << j) < dvstep_sbx[dv]; j++); 
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

void sobseqbx(double *x)
{
    unsigned long long im,imcum;
    unsigned int j,k,l,dv;

    for (k = 0; k < ndim_sbx; k++)
        x[k]=ix_sbx[k]*fac_sbx;

    imcum = in_sbx;
    for (dv = 0; dv < ndv_sbx; dv++) {

        it_sbx +=dvstep_sbx[dv];
        if (dvstep_sbx[dv] == 1) im = it_sbx;
        else {
            im = imcum;
            if (dv != ndv_sbx-1) 
               while (imcum < it_sbx) 
                    imcum += dvstep_sbx[dv+1];
	}
#if 0
printf(" a: %d %d",(int)im,(int)it_sbx);
#endif
        NONZEROBIT(j,im)
        l = j*MAXDIM + 1;
        for (k = 0; k < ndim_sbx; k++) 
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

