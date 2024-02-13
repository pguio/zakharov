
#include <blitz/array.h>
#include <blitz/array/indirect.h>
#include <list>
#include <sys/times.h>
#if defined CLK_TCK
#undef CLK_TCK
#endif
#define CLK_TCK sysconf(_SC_CLK_TCK)


int main()
{
  static const int N=4096;
  static const int M=N/3;
  blitz::Array<double,1> A(N);
  blitz::Array<double,1> B(N);
  blitz::Array<double,1> C(N);
  blitz::Array<double,1> D(N);

  struct tms start, stop;

  A=1.0;

  blitz::Range kp, kn;
  kp.setRange(0, M);    // length M+1 : k=0 k=1..M
  kn.setRange(N-M,N-1);

  B=0.0;
  times(&start);
  for (int i=1; i<1000; i++) {
    B(kp) = A(kp);
    B(kn) = A(kn);
  }
  times(&stop);
  std::cout << "time used " << (double)(stop.tms_utime-start.tms_utime)/CLK_TCK << std::endl;

  B=0.0;
  times(&start);
  for (int i=1; i<1000; i++) {
    B = A;
  }
  times(&stop);
  std::cout << "time used " << (double)(stop.tms_utime-start.tms_utime)/CLK_TCK << std::endl;

  std::list<int> I;
  for (int i=0; i<=M; i++)
    I.push_back(i);
  for (int i=N-M; i<=N-1; i++)
    I.push_back(i);

  C=0.0;
  times(&start);
  for (int i=1; i<1000; i++) {
    C[I] = A;
  }
  times(&stop);
  std::cout << "time used " << (double)(stop.tms_utime-start.tms_utime)/CLK_TCK << std::endl;

  D=B-C;
  //std::cout << "B-C=" << D << std::endl;

  return 0;
}
