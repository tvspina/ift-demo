#include "ift.h"

void ReadSeeds(char *filename, Set **Obj, Set **Bkg)
{
  FILE *fp=fopen(filename,"r");
  int i,x,y,l,nseeds, ncols, nrows;
  
  if(fscanf(fp,"%d %d %d",&nseeds, &ncols, &nrows)!=0);

  for (i=0; i < nseeds; i++){
    if(fscanf(fp,"%d %d %d",&x,&y,&l)!=0);
    if (l==0)
      InsertSet(Bkg, x + ncols*y);
    else
      InsertSet(Obj, x + ncols*y);
  }
  fclose(fp);
}

// Watershed from binary marker

Image *Watershed(Image *img, Set *Obj, Set *Bkg)
{
  AdjRel *A=NULL;
  GQueue *Q=NULL;
  Image  *cost=NULL,*label=NULL;
  Pixel   u,v;
  int     i,p,q,n,tmp,Cmax=MaximumValue(img);
  Set    *S;

  cost  = CreateImage(img->ncols,img->nrows);
  label = CreateImage(img->ncols,img->nrows);
  n     = img->ncols*img->nrows;
  Q     = CreateGQueue(Cmax+1,n,cost->val);
  A     = Circular(1.5);

  /* Trivial path initialization */

  for (p=0; p < n; p++){
    cost->val[p] =INT_MAX;
  }
  S = Obj;
  while(S != NULL){
    p=S->elem;
    label->val[p]=1;
    cost->val[p]=0;
    InsertGQueue(&Q,p);
    S = S->next;
  }
  S = Bkg;
  while(S != NULL){
    p=S->elem;
    label->val[p]=0;
    cost->val[p]=0;
    InsertGQueue(&Q,p);
    S = S->next;
  }

  /* Path propagation */

  while (!EmptyGQueue(Q)){
    p   = RemoveGQueue(Q);
    u.x = p%img->ncols;
    u.y = p/img->ncols;
    for (i=1; i < A->n; i++) {
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if (ValidPixel(img,v.x,v.y)){
	q   = v.x + img->tbrow[v.y];
	if (cost->val[p] < cost->val[q]){
	  
	  tmp = MAX(cost->val[p] , img->val[q]);
	  if (tmp < cost->val[q]){
	    if (cost->val[q]!=INT_MAX)
	      RemoveGQueueElem(Q,q);
	    cost->val[q] =tmp;
	    label->val[q]=label->val[p];
	    InsertGQueue(&Q,q);
	  }
	}

      }
    }
  }

 
  DestroyGQueue(&Q);
  DestroyImage(&cost);
  DestroyAdjRel(&A);

  return(label);
}

int main(int argc, char **argv) 
{
  timer    *t1=NULL,*t2=NULL;
  Image    *img=NULL,*grad=NULL;
  Image    *label=NULL;
  CImage   *cimg=NULL;
  Set      *Obj=NULL,*Bkg=NULL;

  /*--------------------------------------------------------*/

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;

  /*--------------------------------------------------------*/

  if (argc!=4){
    fprintf(stderr,"Usage: watershed <image.pgm> <gradient.pgm> <seeds.txt>\n");
    fprintf(stderr,"image.pgm: image to overlay the watershed lines on it\n");
    fprintf(stderr,"gradient.pgm: gradient image to compute the watershed segmentation\n");
    fprintf(stderr,"seeds.txt: seed pixels\n");
    exit(-1);
  }

  img   = ReadImage(argv[1]);
  grad  = ReadImage(argv[2]);
  ReadSeeds(argv[3],&Obj,&Bkg);

  t1 = Tic();
  
  label = Watershed(grad,Obj,Bkg);
  
  t2 = Toc();    
  

  fprintf(stdout,"Processing time in %f ms\n",CTime(t1,t2));
  
  cimg = DrawLabeledRegions(img,label);
  WriteCImage(cimg,"result.ppm");    
  DestroyImage(&grad);  
  DestroyImage(&img);  
  DestroyImage(&label);
  DestroyCImage(&cimg);  
  DestroySet(&Obj);
  DestroySet(&Bkg);

  /* ---------------------------------------------------------- */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   

  return(0);
}
