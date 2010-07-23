#include "ift.h"

// Add value to remove basins with depth less than H

Image *AddValue(Image *img, int H)
{
  Image *marker=CreateImage(img->ncols,img->nrows);
  int p,n=img->ncols*img->nrows;

  for (p=0; p < n; p++) 
    marker->val[p]=img->val[p]+H;
    
  return(marker);
}

// Watershed from grayscale marker

Image *WaterGray(Image *img, Image *marker, AdjRel *A)
{
  Image *cost=NULL,*label=NULL, *pred=NULL;
  GQueue *Q=NULL;
  int i,p,q,tmp,n,lambda=1;
  Pixel u,v;

  n     = img->ncols*img->nrows;
  cost  = CreateImage(img->ncols,img->nrows);
  label = CreateImage(img->ncols,img->nrows);
  pred  = CreateImage(img->ncols,img->nrows);
  Q     = CreateGQueue(MaximumValue(marker)+2,n,cost->val);

  // Trivial path initialization

  for (p=0; p < n; p++) {
    cost->val[p]=marker->val[p]+1; 
    pred->val[p]=NIL;
    InsertGQueue(&Q,p);
  }

  // Path propagation

  while(!EmptyGQueue(Q)) {
    p=RemoveGQueue(Q);
    if (pred->val[p]==NIL) { // on-the-fly root detection
      cost->val[p] =img->val[p];
      label->val[p]=lambda; lambda++;
    }
    u.x = p%img->ncols;
    u.y = p/img->ncols;
    for (i=1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if (ValidPixel(img,v.x,v.y)){
	q = v.x + img->tbrow[v.y];
	if (cost->val[q] > cost->val[p]){
	  tmp = MAX(cost->val[p],img->val[q]);
	  if (tmp < cost->val[q]){
	    RemoveGQueueElem(Q,q);
	    pred->val[q]  = p;
	    label->val[q] = label->val[p];
	    cost->val[q]  = tmp;
	    InsertGQueue(&Q,q);
	  }
	}
      }
    }
  }
  DestroyGQueue(&Q);
  DestroyImage(&cost);
  DestroyImage(&pred);

  return(label);
}

int main(int argc, char **argv) 
{
  timer    *t1=NULL,*t2=NULL;
  Image    *img=NULL,*grad=NULL,*marker=NULL;
  Image    *label=NULL;
  CImage   *cimg=NULL;
  AdjRel   *A=NULL;

  /*--------------------------------------------------------*/

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;

  /*--------------------------------------------------------*/

  if (argc!=3){
    fprintf(stderr,"Usage: watergray <image.pgm> <gradient.pgm>\n");
    fprintf(stderr,"image.pgm: grayscale image to overlay the watershed lines on it\n"); 
    fprintf(stderr,"gradient.pgm: gradient image to compute the watershed segmentation\n"); 
    exit(-1);
  }
 
  img    = ReadImage(argv[1]);
  grad   = ReadImage(argv[2]);

  // A grayscale marker can be created by any extensive operation:
  // A value H may be added to eliminate irrelevant basins, for instance. 

  marker = AddValue(grad,50);

  // Watershed from grayscale marker

  A = Circular(1.0); // try also higher adjacency radii: 1.5, 2.5, etc.
    
  t1 = Tic();

  label = WaterGray(grad,marker,A);

  t2 = Toc();    
  
  fprintf(stdout,"Processing time in %f ms\n",CTime(t1,t2));


  // Draw watershed lines

  cimg = DrawLabeledRegions(img,label);
  WriteCImage(cimg,"result.ppm");    
  DestroyImage(&grad);  
  DestroyImage(&img);  
  DestroyImage(&marker);  
  DestroyCImage(&cimg);
  DestroyImage(&label);
  DestroyAdjRel(&A);


  /* ---------------------------------------------------------- */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   

  return(0);
}




