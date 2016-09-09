#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include "mess.hpp" 
#include "sobseq.hpp"

using namespace std;

#define NOOP 0

#define SECS (1.0/(double)(CLOCKS_PER_SEC))

#define PEXIT(x) {printf0(" error: %s.\n",x);exit(255);}
#define CHK_ALL_PROCS(x) (!all_procs_OK(!(x)))
#define ERR_CHK(x,y) {if(CHK_ALL_PROCS(x)) PEXIT(y)}

/* print-n-get */
#define PNGS(a,b) { printf0(a); \
    if (gets0(b) == NULL) PEXIT("gets() failed"); }
#define PNG(a,b,c) { char tmp[BUFSIZ]; PNGS(a,tmp) \
    if (sscanf(tmp,b,c) < 1) PEXIT("scanf failed") }

#define XFER(v,n) { if (myproc==0) broadcast((void*)v,n); \
                  else recvbcast((void*)v,n,0); }

extern double integrand(double *);
extern double int_volume;

static int myproc,nproc;

void sobseqbq_init(int,int,int);
void sobseqbq(double *);

void sobseqbx_init(int,int,int);
void sobseqbx(double *);

#define MAXTEST 1000000
#define DIM 6

int main(int argc, char *argv[])
{
    int i,j,dim;
    double q[DIM];
    double p[DIM];
    double s;

    double twc = clock()*SECS;

    start_mess(&argc,&argv);
    myproc = getmyproc();
    nproc = getnproc();

    printf0(" *** test quasirand generator -- sobol\' sequence method ***\n");

    /* set up */
    dim = 3;

    printf0(" generating %d pts in %d dim\n",MAXTEST,dim);

    double tsob = clock()*SECS;
    
    for (s = 0., i = 0; i < MAXTEST/nproc; i++) {
        sobseqb(q,dim,nproc,myproc,i);
    }
    barrier();
    tsob = clock()*SECS - tsob;
    

    sobseqbq_init(dim,nproc,myproc);

    double tsobq = clock()*SECS;
    for (s = 0., i = 0; i < MAXTEST/nproc; i++) {
        sobseqbq(q);
    }
    barrier();
    tsobq = clock()*SECS - tsobq;
        
    twc = clock()*SECS - twc;
    printf0(" total wallclock time: %g\n",twc);
    printf0(" time generating sobseq (x): %g\n",tsob);
    printf0(" time generating sobseq (q): %g\n",tsobq);

    end_mess();

    sobseqbq_init(dim,1,0);
    sobseqbx_init(dim,nproc,myproc);
    for (s = 0., i = 0; i < MAXTEST; i++) {
        sobseqbq(p);
if (i < 10)
        printf(" %g ",p[0]);
        if (i % nproc == myproc) {
            sobseqbx(q);
if (i < 10)
        printf(" %g",q[0]);
            
            for (j = 0; j < dim; j++)
                s += (p[j]-q[j])*(p[j]-q[j]);
	}
if (i < 10)
        printf("\n");
    }
    printf0(" diff^2: %g\n",s);
    return 0;
}
