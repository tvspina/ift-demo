#include "common.h"

char *AllocCharArray(int n)
{
  char *v=NULL;
  v = (char *) calloc(n,sizeof(char));
  if (v == NULL)
    Error(MSG1,"AllocCharArray");
  return(v);
}

uchar *AllocUCharArray(int n)
{
  uchar *v=NULL;
  v = (uchar *) calloc(n,sizeof(uchar));
  if (v == NULL)
    Error(MSG1,"AllocUCharArray");
  return(v);
}

ushort *AllocUShortArray(int n)
{
  ushort *v=NULL;
  v = (ushort *) calloc(n,sizeof(ushort));
  if (v == NULL)
    Error(MSG1,"AllocUShortArray");
  return(v);
}

uint *AllocUIntArray(int n)
{
  uint *v=NULL;
  v = (uint *) calloc(n,sizeof(uint));
  if (v == NULL)
    Error(MSG1,"AllocUIntArray");
  return(v);
}

int *AllocIntArray(int n)
{
  int *v=NULL;
  v = (int *) calloc(n,sizeof(int));
  if (v == NULL)
    Error(MSG1,"AllocIntArray");
  return(v);
}

float *AllocFloatArray(int n)
{
  float *v=NULL;
  v = (float *) calloc(n,sizeof(float));
  if (v == NULL)
    Error(MSG1,"AllocFloatArray");
  return(v);
}

double *AllocDoubleArray(int n)
{
  double *v=NULL;
  v = (double *) calloc(n,sizeof(double));
  if (v == NULL)
    Error(MSG1,"AllocDoubleArray");
  return(v);
}

real  *AllocRealArray(int n){
  real *v=NULL;
  v = (real *) calloc(n,sizeof(real));
  if (v == NULL)
    Error(MSG1,"AllocRealArray");
  return(v);
}

void Error(char *msg,char *func){ /* It prints error message and exits
                                    the program. */
  fprintf(stderr,"Error:%s in %s\n",msg,func);
  exit(-1);
}

void Warning(char *msg,char *func){ /* It prints warning message and
                                       leaves the routine. */
 fprintf(stdout,"Warning:%s in %s\n",msg,func);

}

void Change(int *a, int *b){ /* It changes content between a and b */
  int c;
  c  = *a;
  *a = *b;
  *b = c;
}

void FChange(float *a, float *b){ /* It changes content between floats a and b */
  float c;
  c  = *a;
  *a = *b;
  *b = c;
}

int NCFgets(char *s, int m, FILE *f) {
  while(fgets(s,m,f)!=NULL)
    if (s[0]!='#') return 1;
  return 0;
}


/**
 * Gera um número inteiro aleatório no intervalo [low,high].
http://www.ime.usp.br/~pf/algoritmos/aulas/random.html
 **/
int RandomInteger (int low, int high){
  int k;
  double d;

  d = (double) rand () / ((double) RAND_MAX);
  k = d * (high - low);
  return low + k;
}


inline real sqrtreal(real x){
  int size;

  size = sizeof(real);
  if(size == sizeof(float))
    return sqrtf(x);
  else if(size == sizeof(double))
    return sqrt(x);
  else
    return (real)sqrt((double)x);
}

int SafeMod(int a, int n)
{
	int r = a % n;

	return (r >= 0) ? r : n+r;
}

int IsPowerOf2(int x)
{
    return (x&(x-1)) == 0;
}

