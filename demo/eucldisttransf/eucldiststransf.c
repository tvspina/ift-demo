#include "ift.h"


// Creates Empty Forest

typedef struct _forest {
  Image *P; // predecessor map
  Image *R; // root map
  Image *V; // distance (cost or connectivity) map
} Forest;

Forest *CreateForest(int ncols, int nrows)
{
  Forest *F=(Forest *)calloc(1,sizeof(Forest));

  F->P = CreateImage(ncols,nrows);
  F->R = CreateImage(ncols,nrows);
  F->V = CreateImage(ncols,nrows);

  return(F);
}

// Destroys Forest

void DestroyForest(Forest **F)
{
  Forest *tmp=*F;

  if (tmp != NULL) {
    DestroyImage(&(tmp->P));
    DestroyImage(&(tmp->R));
    DestroyImage(&(tmp->V));
    free(tmp);
    *F = NULL;
  }
}

Image  *CopyImage(Image *img)
{
    Image *imgc;

    imgc = CreateImage(img->ncols,img->nrows);
    memcpy(imgc->val,img->val,img->ncols*img->nrows*sizeof(int));

    return(imgc);
}


Image *Threshold(Image *img, int lower, int higher)
{
  Image *bin=NULL;
  int p,n;

  bin = CreateImage(img->ncols,img->nrows);
  n = img->ncols*img->nrows;
  for (p=0; p < n; p++)
    if ((img->val[p] >= lower)&&(img->val[p] <= higher))
      bin->val[p]=1;
  return(bin);
}


// Euclidean distance transform

Forest *DistTrans(Image *I) 
{
  int p,q,n=I->ncols*I->nrows,i,tmp;
  Pixel u,v,w;
  AdjRel *A=Circular(1.5),*A4=Circular(1.0);
  Forest *F=CreateForest(I->ncols,I->nrows);
  GQueue *Q=CreateGQueue(1024,n,F->V->val);
  
  // Trivial path initialization

  for (p=0; p < n; p++) {
    u.x = p % I->ncols;
    u.y = p / I->ncols;
    F->V->val[p]=INT_MAX; F->R->val[p]=p; F->P->val[p]=NIL;
    if (I->val[p]!=0){ // p belongs to an object's border
      F->V->val[p]=0;
      InsertGQueue(&Q,p);
    }
  }

  // Path propagation

  while(!EmptyGQueue(Q)){
    p = RemoveGQueue(Q);
    u.x = p % I->ncols;
    u.y = p / I->ncols;
    w.x = F->R->val[p] % I->ncols;
    w.y = F->R->val[p] / I->ncols;
    for (i=1; i < A->n; i++) {
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if (ValidPixel(I,v.x,v.y)){
	q   = v.x + I->tbrow[v.y];
	if (F->V->val[q]>F->V->val[p]){	    
	  tmp = (v.x-w.x)*(v.x-w.x)+(v.y-w.y)*(v.y-w.y);
	  if (tmp < F->V->val[q]){
	    if (F->V->val[q]!=INT_MAX) RemoveGQueueElem(Q, q);
	    F->V->val[q]=tmp; F->R->val[q]=F->R->val[p]; F->P->val[q]=p;
	    InsertGQueue(&Q,q);
	  }
	}
      }
    }
  }

  DestroyGQueue(&Q);
  DestroyAdjRel(&A);
  DestroyAdjRel(&A4);

  return(F);
}

int main(int argc, char **argv)
{
  timer    *t1=NULL,*t2=NULL;
  Image    *img,*aux;
  Forest   *tde;
  /* The following block must the remarked when using non-linux machines */

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  
  /*----------------------------------------------------------------------*/
  
  if (argc != 2) {
    printf("Usage: %s <image.pgm>\n", argv[0]);
    exit(0);
  }

  aux = ReadImage(argv[1]);

  if (MaximumValue(aux)!=1){
    fprintf(stderr,"Input image must be binary with values 0/1 \n");
    fprintf(stderr,"Assuming lower threshold 100 for this conversion\n");
    img = Threshold(aux,100,INT_MAX);
    WriteImage(img,"shape.pgm");
  }else{
    img = CopyImage(aux);
  }
  DestroyImage(&aux);
    
  t1 = Tic();

  tde = DistTrans(img);

  t2 = Toc();

  fprintf(stdout,"Euclidian Distance Transform in %f ms\n",CTime(t1,t2));

  WriteImage(tde->V,"tde.pgm");

  DestroyForest(&tde);
  DestroyImage(&img);

  /* The following block must the remarked when using non-linux machines */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   
  

  return(0);
}
