#ifndef _mess_h
#define _mess_h 1

#ifndef CONCURRENT
#define CONCURRENT 1
#endif

#if CONCURRENT

/*
#include <mpi++.h> 
*/
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

const int MESS_MAX_PACK = 16384;
const int MESS_MAX_PROC = 2048;

inline void begin_mess(int* argc, char ***argv)
{
    MPI_Init(argc,argv);
}

inline int start_mess(int *argc, char ***argv)
{
    MPI_Init(argc,argv);
    return 1;
}

inline int end_mess(void)
{
    MPI_Finalize();
    return 1;
}


inline void exitall(const int i) 
{
  end_mess();
  exit(i);
}


/*** pe info ***/

inline int getmyproc(void)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    return rank;
}

inline int getnproc(void)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    return size;
}

/*** execution &/| wall time ***/

inline double time_mess(void)
{
    return MPI_Wtime();
}

/*** barrier ***/


inline void barrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
}

/*** inter-pe communication ***/

inline int transmit_s(void *buf,int nbytes,int pe) 
/* send buf[nbytes] to PE # pe. ret nbyte if OK else error code <= 0. */
{
    if (MPI_Send(buf,nbytes,MPI_BYTE,pe,0,MPI_COMM_WORLD) < 0) 
        return -1;
    return nbytes;
}

inline int receive_s(void *buf, int nbytes, int pe)
/* generic blocking receive, returns nbytes recv'd, or <= 0 on error */
{
    MPI_Status status;

    if (MPI_Recv(buf,nbytes,MPI_BYTE,pe,0,MPI_COMM_WORLD,&status) < 0)
        return -1;
    return nbytes;
}

inline int broadcast_s(void *buf, int nbytes)
/* generic broadcast-to-all function */
{
    if (MPI_Bcast(buf,nbytes,MPI_BYTE,getmyproc(),MPI_COMM_WORLD) < 0) 
        return -1;
    return nbytes;
}

inline int recvbcast_s(void *buf, int nbytes, int pe)
/* blocking recv for a broadcast to all procs (less the xmitter) */
{
    if (MPI_Bcast(buf,nbytes,MPI_BYTE,pe,MPI_COMM_WORLD) < 0)
        return -1;
    return nbytes;
}

template<class T>
inline int spread(T& x, const int root)
     // bcast + bcast_recv; should use raw MPI here....
{
    T q=x;
    int err;
    int myproc = getmyproc();
    if (myproc == root) err = broadcast_s(&q,sizeof(T));
    else err = recvbcast_s(&q,sizeof(T),root);			  
    x=q;
    return err;
}

template<class T> 
inline int spread(T* x, const int n, const int root)
{
    int err;
    int myproc = getmyproc();
    if (myproc == root) err = broadcast_s((void*)x,n*sizeof(T));
    else err = recvbcast_s((void*)x,n*sizeof(T),root);
    return err;
}

inline int local_count(const int ntotal, const int proc, const int nproc)
{
    int dsi,rem;
    dsi = ntotal / nproc;
    rem = ntotal % nproc;
    return (proc < rem) ? dsi + 1 : dsi;
}

inline int global_offset(const int ntotal, const int proc, const int nproc)
{
    int dsi,rem,ret;
    dsi = ntotal / nproc;
    rem = ntotal % nproc;
    if (proc < rem)
        ret = proc * (dsi+1);
    if (proc >= rem)
        ret = rem * (dsi+1) + (proc - rem) * dsi;
    return ret;    
}

// eventually turn these into a single template:

inline
int combine_sum(double *p, const int n)
     // unbuffered sum-over-procs....
{
    int myproc = getmyproc();
    double* q = new double [n];
    int err = MPI_Reduce(p,q,n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if (myproc == 0) {
      for (int i=0; i<n; ++i) p[i]=q[i];
    }
    delete [] q;
    return err;
}

inline
int combine_sum(int *p, const int n)
{
    int myproc = getmyproc();
    int* q = new int [n];
    int err = MPI_Reduce(p,q,n,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if (myproc == 0) {
	for (int i=0; i<n; ++i) p[i]=q[i];
    }
    delete [] q;
    return err;
}

inline
int combine_sum(long *p, const int n)
{
    int myproc = getmyproc();
    long* q = new long [n];
    int err = MPI_Reduce(p,q,n,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    if (myproc == 0) {
	for (int i=0; i<n; ++i) p[i]=q[i];
    }
    delete [] q;
    return err;
}


// single element combines:

template<class T>
inline int combine_sum(T& p)
     // xxx- should fix this, need to know how to do a conditional on T
{
    T q=p;
    int err=combine_sum(&q,1);
    p=q;
    return err;
}


// highly specific gather-scatter:

inline int combine_merge(double* p, const int nloc, const int ntot)
     // nloc,ntot must be same on all procs!
{
    int myproc = getmyproc();
    int nproc = getnproc();
    double *x = new double [ntot];
    for (int m=0; m<nproc; ++m) {
	for (int i=0; i<nloc; ++i) {
	    x[i+m*nloc]=p[i+myproc*nloc];
	}
    }
    int ret(MPI_Alltoall(x,nloc,MPI_DOUBLE,p,nloc,MPI_DOUBLE,MPI_COMM_WORLD));
    delete [] x;
    return ret;
}

template <class T>
inline int combine_merge(T* p, const int nloc, const int ntot)
     // nloc,ntot must be same on all procs!
{
    int myproc = getmyproc();
    int nproc = getnproc();
    T* x = new T [ntot];
    for (int m=0; m<nproc; ++m) {
	for (int i=0; i<nloc; ++i) {
	    x[i+m*nloc]=p[i+myproc*nloc];
	}
    }
    int nbytes = sizeof(T)*nloc;
    int ret(MPI_Alltoall((void*)x,nbytes,MPI_BYTE,
		     (void*)p,nbytes,MPI_BYTE,MPI_COMM_WORLD));
    delete [] x;
    return ret;
}     


template <class T>
inline int scatter(T* p, const int ntot, const int root)
// assume p[ntot] on root, and that space is available already on other 
// nodes to handle scattered data...
// NB: data on root is NOT resized... 
{
  int myproc = getmyproc();
  int nproc = getnproc();
  int mycnt = local_count(ntot,myproc,nproc);

  if (ntot>0) {
    if (myproc == root) {
      for (int proc = 0; proc<nproc; proc++) {
	if (proc == myproc) continue;
	int thisoff = global_offset(ntot,proc,nproc);
	int thiscnt = local_count(ntot,proc,nproc);
	transmit_s(&p[thisoff],thiscnt*sizeof(T),proc);
      }
    } else {
      receive_s(&p[0],mycnt*sizeof(T),root);
    }
  }  
  return mycnt;
}

template <class T>
inline int gather(T* p, const int ntot, const int root)
// assume p[ntot] on root, has space to receive rank-ordered array data..
{
  int myproc = getmyproc();
  int nproc = getnproc();
  int mycnt = local_count(ntot,myproc,nproc);

  if (ntot>0) {
    if (myproc == root) {
      for (int proc = 0; proc<nproc; proc++) {
	if (proc == myproc) continue;
	int thisoff = global_offset(ntot,proc,nproc);
	int thiscnt = local_count(ntot,myproc,nproc);
	receive_s(&p[thisoff],thiscnt*sizeof(T),proc);
      }
    } else {
      transmit_s(&p[0],mycnt*sizeof(T),0);
    }
  }  
  return mycnt;
}

inline
double drand48all()
{
  double rnd;
  if (getmyproc() == 0) rnd = drand48(); 
  spread(rnd,0);
  return rnd;
}

inline
void srand48all(const long int seed)
{
  long int sd = seed;
  spread(sd,0);
  srand48(sd);
  return;
}     


template <class T> inline int compareall(T *p, const size_t n)
{  
  T *q = new T [n]; 
  for (size_t i=0; i<n; i++) q[i] = p[i];
  int not_same = 0;
  for (int myp = 1; myp<getnproc(); myp++) {
    if (getmyproc() == myp) {
      transmit_s((void*)q,n*sizeof(T),0);
    } else if (getmyproc()==0) {  
      receive_s((void*)q,n*sizeof(T),myp);
      for (size_t i=0; i<n; i++) {
	if (p[i] != q[i]) { 
	  not_same++; 
	  break;
	}
      }
    }
  }
  delete [] q;
  spread(not_same,0);
  return not_same;
}


template <class T> inline int compareall(T& p)
{  
  T q = p; 
  int not_same = 0;
  for (int myp = 1; myp<getnproc(); myp++) {
    if (getmyproc() == myp) {
      transmit_s((void*)&q,sizeof(T),0);
    } else if (getmyproc()==0) {  
      receive_s((void*)&q,sizeof(T),myp);
      if (p != q) not_same++; 
    }
  }
  spread<int>(not_same,0);
  return not_same;
}


template <class T> inline T maxall(const T p)
{  
  T pmax = p;
  for (int myp = 1; myp<getnproc(); myp++) {
    if (getmyproc() == myp) {
      transmit_s((void*)&pmax,sizeof(T),0);
    } else if (getmyproc()==0) {  
      T q=p; 
      receive_s((void*)&q,sizeof(T),myp);
      if (q>pmax) pmax=q;
    }
  }  
  spread(pmax,0);
  return pmax;
}

template <class T> inline T minall(const T p)
{  
  T pmin = p;
  for (int myp = 1; myp<getnproc(); myp++) {
    if (getmyproc() == myp) {
      transmit_s((void*)&pmin,sizeof(T),0);
    } else if (getmyproc()==0) {  
      T q=p; 
      receive_s((void*)&q,sizeof(T),myp);
      if (q<pmin) pmin=q; 
    }
  }
  spread(pmin,0);
  return pmin;
}


#else

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

/* these are the message passing functions when 
   no message passing can be performed... */

inline double time_mess()
// wall time, evidently; not (yet) necessarily cpu time
// this won't work on all systems; ok only if clock() is rel to startup.
// just in case, do try to init'lize but note this will fail 
// if clock is rel to something odd, and if fn is called from a module
// not containing main().
{
  static clock_t clockstart = -1;
  if (clockstart == -1) clockstart = clock();
  return double(clock()-clockstart)/double(CLOCKS_PER_SEC);
}

inline void start_mess(int* a, char*** b) { time_mess(); }
inline void begin_mess(int* a, char*** b) { time_mess(); }
inline void end_mess() {}


inline int getnproc() { return 1; }
inline int getmyproc() { return 0; }


inline int local_count(const int ntotal, const int proc, const int nproc)
{
    int dsi,rem;
    dsi = ntotal / nproc;
    rem = ntotal % nproc;
    return (proc < rem) ? dsi + 1 : dsi;
}

inline int global_offset(const int ntotal, const int proc, const int nproc)
{
    int dsi,rem,ret;
    dsi = ntotal / nproc;
    rem = ntotal % nproc;
    if (proc < rem)
        ret = proc * (dsi+1);
    if (proc >= rem)
        ret = rem * (dsi+1) + (proc - rem) * dsi;
    return ret;    
}

inline void barrier() { return; }

inline int broadcast_s(void *v, int n) { return 0;}
inline int recvbcast_s(void *v, int n, int pe) { return 0; }

template<class T> inline int combine_sum(T *p, const int n) { return 0; }
template<class T> inline int combine_sum(T& p) { return 0; }
template<class T> inline int combine_merge(T* p, const int n, const int m) {
  return 0;
}

template<class T> 
inline int spread(T& x, const int root) { return 0; }
template<class T> 
inline int spread(T* x, const int n, const int root) { return 0; }

template <class T>
inline int scatter(T* p, const int ntot, const int root) { return ntot; }
template <class T>
inline int gather(T* p, const int ntot, const int root) { return ntot; }

inline double drand48all() { return drand48(); }
inline void srand48all(const long int s) { return srand48(s); }

template <class T> inline int compareall(T *p, size_t n) { return 0; }
template <class T> inline int compareall(T &p) { return 0; }
template <class T> inline T maxall(const T p) { return p; }
template <class T> inline T min(const T p) { return p; }

inline int transmit_s(void *x, int i, int j) { return i; }
inline int receive_s(void *x, int i, int j) { return i; }
inline void exitall(const int i) { exit(i); }


#endif

#define in_sequence(a,b) \
   for (int in_seq_proc=0; in_seq_proc<(b); barrier(), in_seq_proc++) \
     if (a==in_seq_proc)


// overloads of the "*" operator to channel output through a single proc

namespace 
{
    std::ofstream zerostream;
}

inline std::ostream& operator*(std::ostream& s, const int& n)
{
    if (n == getmyproc())  return s;
    return (std::ostream&)(zerostream);
}

inline std::ostream& operator++(std::ostream& s, int ii)
{
    barrier();
    for (int n=0; n<getnproc(); n++) {
      if (n == getmyproc())  s.flush();
      barrier();
    }
    return s;
}

inline
int fprintf0(FILE *fp, const char *fmt, ...)
/* fprintf from procnum 0 */
{
    va_list args;
    int cnt = 0;

    if (fmt != NULL) {
        if (getmyproc() == 0) {
            va_start(args, fmt);
            cnt = vfprintf(fp, fmt, args);
            va_end(args);
	}
    }
    return cnt;
}

inline
int printf0(const char *fmt, ...)
/* printf from procnum 0 */
{
    va_list args;
    int cnt = 0;

    if (fmt != NULL) {
        if (getmyproc() == 0) {
            va_start(args, fmt);
            cnt = vprintf(fmt, args);
            va_end(args);
            fflush(stdout);
	}
    }
    return cnt;
}


#endif


