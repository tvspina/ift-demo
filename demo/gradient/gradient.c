#include "ift.h"

Image* InvertValues(Image* img)
{
  const int maxval = MaximumValue(img);
  int p;
  Image* result = CreateImage(img->ncols, img->nrows);
  
  for(p = 0; p < img->ncols*img->nrows; p++)
    result->val[p] = maxval-img->val[p];

  return result;
}

Image* SigmoidalScaling(Image* img, int gain, double cutoff)
{
  const int maxval = MaximumValue(img);
  int p;
  Image* result = CreateImage(img->ncols, img->nrows);

  for(p = 0; p < img->ncols*img->nrows; p++)
    result->val[p] = 255*(1/(1 + exp(gain*(cutoff-(double)img->val[p]/(double)maxval))));
 
 return result;
}


int main(int argc, char **argv)
{
  timer    *t1=NULL,*t2=NULL;
  Image    *img=NULL, *inv_tde=NULL, *scl=NULL;

  /* The following block must the remarked when using non-linux machines */

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  
  /*----------------------------------------------------------------------*/
  
  if (argc != 2) {
    printf("Usage: %s <image.pgm>\n",argv[0]);
    exit(0);
  }

  img = ReadImage(argv[1]);
  inv_tde = InvertValues(img);
  scl = SigmoidalScaling(inv_tde, 10, 0.5);

  t1 = Tic();

  t2 = Toc();

  fprintf(stdout,"Gradient computing in %f ms\n",CTime(t1,t2));

  WriteImage(inv_tde,"inv_tde.pgm");
  WriteImage(scl,"scale.pgm");

  DestroyImage(&inv_tde);
  DestroyImage(&img);
  DestroyImage(&scl);

  /* The following block must the remarked when using non-linux machines */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   
  

  return(0);
}
