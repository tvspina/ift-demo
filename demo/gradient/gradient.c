#include "ift.h"

#define MAXGRAD 255

Image* InvertValues(Image* img)
{
  const int maxval = MaximumValue(img);
  int p;
  Image* result = CreateImage(img->ncols, img->nrows);
  
  for(p = 0; p < img->ncols*img->nrows; p++)
    result->val[p] = maxval-img->val[p];

  return result;
}

Image* SigmoidalStretch(Image* img, int thresh, int gain, double cutoff)
{
  const int maxval = MaximumValue(img);
  const int minval = MinimumValue(img);
  int p;
  Image* result = CreateImage(img->ncols, img->nrows);

  for(p = 0; p < img->ncols*img->nrows; p++)
    {
      if(img->val[p] < -thresh) 
    	result->val[p] = 0;
      else if(img->val[p] > thresh) 
	result->val[p] = maxval;
      else
	result->val[p] = (int)((double)maxval/(1 + exp(gain*((double)(img->val[p]-minval)/(double)(maxval-minval)-cutoff))));
    }
  return result;
}

Image *ObjectGradient(Image *img, float radius)
{
    real    dist,gx,gy;
    int     i,p,q,n=img->ncols*img->nrows,Imax;
    Pixel   u,v;
    AdjRel *A=Circular(radius);
    real   *md=AllocRealArray(A->n);

    Image* grad = CreateImage(img->ncols, img->nrows);

    Imax = MaximumValue(img);

    for (i=1; i < A->n; i++)
        md[i]=sqrt(A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i]);
 
    for (p=0; p < n; p++)
    {
        u.x = p%img->ncols;
        u.y = p/img->ncols;

        gx = gy = 0.0;

        for (i=1; i < A->n; i++)
        {
            v.x = u.x + A->dx[i];
            v.y = u.y + A->dy[i];
            if (ValidPixel(img,v.x,v.y))
            {
                q    = v.x + img->tbrow[v.y];
                dist = ((float)img->val[q]-(float)img->val[p])/(float)Imax;

                gx  += dist*A->dx[i]/md[i];
                gy  += dist*A->dy[i]/md[i];
            }
        }
        grad->val[p]=MAXGRAD*sqrt(gx*gx + gy*gy);
    }

    free(md);
    DestroyAdjRel(&A);

    return(grad);
}

Image *FeaturesGradient(Features *f,float radius)
{
    real    dist,gx,gy,mag;
    int     j,i,p,q,n=f->ncols*f->nrows;
    Pixel   u,v;
    AdjRel *A=Circular(radius);
    real   *md=AllocRealArray(A->n);

    Image* grad = CreateImage(f->ncols, f->nrows);  

    for (i=1; i < A->n; i++)
        md[i]=sqrt(A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i]);

    for (p=0; p < n; p++)
    {
        u.x = p%f->ncols;
        u.y = p/f->ncols;

        float max_mag = FLT_MIN;
        for (j=0; j<f->nfeats; j++)
        {
            gx = gy = 0.0;
            for (i=1; i < A->n; i++)
            {
                v.x = u.x + A->dx[i];
                v.y = u.y + A->dy[i];
                if ((v.x>=0 && v.x<f->ncols) && (v.y>=0 && v.y<f->nrows))
                {
                    q    = v.x + v.y*f->ncols;
                    dist = (f->elem[q].feat[j]-f->elem[p].feat[j]);
                    gx  += dist*A->dx[i]/md[i];
                    gy  += dist*A->dy[i]/md[i];
                }
            }
            mag = sqrt(gx*gx + gy*gy);

            if (mag > max_mag)
                max_mag = mag;
        }
        grad->val[p] = (int)MAXGRAD*max_mag;
    }

    free(md);
    DestroyAdjRel(&A);

    return(grad);
}

Image* CombineGradients(Image *objgrad, Image *imggrad, float wobj)
{
  int p;
  Image *grad = CreateImage(objgrad->ncols, objgrad->nrows);

  for(p = 0; p < objgrad->ncols*objgrad->nrows; p++)
    grad->val[p] = wobj*objgrad->val[p] + (1-wobj)*imggrad->val[p];

  return grad;
}

Image *Dilate(Image *img, AdjRel *A)
{
  Image *dil=CreateImage(img->ncols,img->nrows);
  int p,q,i;
  Pixel u,v;

  for (u.y=0; u.y < img->nrows; u.y++) 
    for (u.x=0; u.x < img->ncols; u.x++) {
      p = u.x + img->tbrow[u.y];
      dil->val[p]=img->val[p];
      for (i=1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(dil,v.x,v.y)){
	  q = v.x + img->tbrow[v.y];
	  if (img->val[q]>dil->val[p])
	    dil->val[p]=img->val[q];
	}
      }
    }
  return(dil);
}

Image *Erode(Image *img, AdjRel *A)
{
  Image *ero=CreateImage(img->ncols,img->nrows);
  int p,q,i;
  Pixel u,v;

  for (u.y=0; u.y < img->nrows; u.y++) 
    for (u.x=0; u.x < img->ncols; u.x++) {
      p = u.x + img->tbrow[u.y];
      ero->val[p]=img->val[p];
      for (i=1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(ero,v.x,v.y)){
	  q = v.x + img->tbrow[v.y];
	  if (img->val[q]<ero->val[p])
	    ero->val[p]=img->val[q];
	}
      }
    }
  return(ero);
}

Image *Open(Image *img, AdjRel *A)
{
  Image *open=NULL,*ero=NULL;

  ero  = Erode(img,A);
  open = Dilate(ero,A);
  DestroyImage(&ero);

  return(open);
}

Image *Close(Image *img, AdjRel *A)
{
  Image *close=NULL,*dil=NULL;

  dil   = Dilate(img,A);
  close = Erode(dil,A);
  DestroyImage(&dil);

  return(close);
}

Image *SupRec(Image *img, Image *marker)
{
  Image *cost=NULL;
  GQueue *Q=NULL;
  int i,p,q,tmp,n;
  Pixel u,v;
  AdjRel *A=Circular(1.0);

  n     = img->ncols*img->nrows;
  cost  = CreateImage(img->ncols,img->nrows);
  Q     = CreateGQueue(MaximumValue(marker)+1,n,cost->val);

  // Trivial path initialization

  for (p=0; p < n; p++) {
    cost->val[p]=marker->val[p]; 
    InsertGQueue(&Q,p);
  }

  // Path propagation 

  while(!EmptyGQueue(Q)) {
    p=RemoveGQueue(Q);
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
	    cost->val[q]  = tmp;
	    InsertGQueue(&Q,q);
	  }
	}
      }
    }
  }
  DestroyGQueue(&Q);
  DestroyAdjRel(&A);
  return(cost);
}

// Inferior Reconstruction img >= marker without marker imposition

Image *InfRec(Image *img, Image *marker)
{
  Image *cost=NULL;
  GQueue *Q=NULL;
  int i,p,q,tmp,n;
  Pixel u,v;
  AdjRel *A=Circular(1.0);

  n     = img->ncols*img->nrows;
  cost  = CreateImage(img->ncols,img->nrows);
  Q     = CreateGQueue(MaximumValue(img)+1,n,cost->val);
  SetRemovalPolicy(Q,MAXVALUE); 

  // Trivial path initialization

  for (p=0; p < n; p++) {
    cost->val[p]=marker->val[p]; 
    InsertGQueue(&Q,p);
  }

  // Path propagation 

  while(!EmptyGQueue(Q)) {
    p=RemoveGQueue(Q);
    u.x = p%img->ncols;
    u.y = p/img->ncols;
    for (i=1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if (ValidPixel(img,v.x,v.y)){
	q = v.x + img->tbrow[v.y];
	if (cost->val[q] < cost->val[p]){
	  tmp = MIN(cost->val[p],img->val[q]);
	  if (tmp > cost->val[q]){
	    RemoveGQueueElem(Q,q);
	    cost->val[q]  = tmp;
	    InsertGQueue(&Q,q);
	  }
	}
      }
    }
  }
  DestroyGQueue(&Q);
  DestroyAdjRel(&A);

  return(cost);
}

Image *CloseRec(Image *I, AdjRel *A)
{
  Image *J,*K;

  J = Close(I,A);
  K = SupRec(I,J);
  DestroyImage(&J);
  return(K);
}

Image *OpenRec(Image *I, AdjRel *A)
{
  Image *J,*K;

  J = Open(I,A);
  K = InfRec(I,J);
  DestroyImage(&J);
  return(K);
}


Features *MorphImageFeats(Image *img, int nscales)
{
    Features *f=CreateFeatures(img->ncols,img->nrows,nscales);
    AdjRel   *A=NULL;
    int       s,i;
    Image    *img1=NULL,*open=NULL;

    f->Imax = MaximumValue(img);

    for (s=1; s <= nscales; s=s+1)
    {
        A    = Circular(s);
	open = OpenRec(img,A);
	img1 = CloseRec(open,A);
	
        for (i=0; i < f->nelems; i++)
        {
            f->elem[i].feat[s-1] = (float)img1->val[i]/(float)f->Imax;
        }
        DestroyImage(&img1);
	DestroyImage(&open);
        DestroyAdjRel(&A);
    }

    return(f);
}

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


// Euclidean distance transform

Forest *DistTrans(Image *B, Image *I) 
{
  int p,q,n=B->ncols*B->nrows,i,tmp;
  Pixel u,v,w;
  AdjRel *A=Circular(1.5),*A4=Circular(1.0);
  Forest *F=CreateForest(B->ncols,B->nrows);
  GQueue *Q=CreateGQueue(1024,n,F->V->val);
  
  // Trivial path initialization

  for (p=0; p < n; p++) {
    u.x = p % B->ncols;
    u.y = p / B->ncols;
    F->V->val[p]=INT_MAX; F->R->val[p]=p; F->P->val[p]=NIL;
    if (B->val[p]!=0){ // p belongs to an object's border
      F->V->val[p]=0;
      InsertGQueue(&Q,p);
    }
  }

  // Path propagation

  while(!EmptyGQueue(Q)){
    p = RemoveGQueue(Q);
    u.x = p % B->ncols;
    u.y = p / B->ncols;
    w.x = F->R->val[p] % B->ncols;
    w.y = F->R->val[p] / B->ncols;
    for (i=1; i < A->n; i++) {
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if (ValidPixel(B,v.x,v.y)){
	q   = v.x + B->tbrow[v.y];
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
  for(p = 0; p < I->ncols*I->nrows; p++)
    {
      F->V->val[p] = (int)sqrt((double)F->V->val[p]);
      if(I->val[p] == 0) F->V->val[p] = -F->V->val[p];
    }
  DestroyGQueue(&Q);
  DestroyAdjRel(&A);
  DestroyAdjRel(&A4);

  return(F);
}

Image *Boundary(Image *bin)
{
  Image *bndr=NULL;
  int p=0,q,i,n;
  AdjRel *A;
  Pixel u,v;

  A     = Circular(1.5);
  n     = bin->ncols*bin->nrows;
  bndr  = CreateImage(bin->ncols,bin->nrows);
  for (p=0; p < n; p++)
    {
      if (bin->val[p]>0)
	{
	  u.x = p%bin->ncols;
	  u.y = p/bin->ncols;
	  for (i=1; i < A->n; i++)
	    {
	      v.x = u.x + A->dx[i];
	      v.y = u.y + A->dy[i];
	      if (ValidPixel(bin,v.x,v.y))
		{
		  q = v.x + bin->tbrow[v.y];
		  if (bin->val[q]==0)
		    {
		      bndr->val[p]=1;
		      break;
		    }
		}
	    }
	}
    }
  DestroyAdjRel(&A);
  
  return bndr;
}

int main(int argc, char **argv)
{
  float wobj = 0.5;
  timer    *t1=NULL,*t2=NULL;
  Image    *img=NULL, *bndr=NULL, *label=NULL, *objmap=NULL;
  Image    *objgrad=NULL, *imggrad=NULL, *grad=NULL;
  Features *feat=NULL;
  Forest   *tde=NULL;

  /* The following block must the remarked when using non-linux machines */

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  
  /*----------------------------------------------------------------------*/
  
  if (argc != 4) {
    printf("Usage: %s <image.pgm> <label.pgm> <wobj [0,1]>\n",argv[0]);
    exit(0);
  }

  img = ReadImage(argv[1]);
  label = ReadImage(argv[2]);
  wobj = atof(argv[3]);

  t1 = Tic();

  /* computing object gradient */
  bndr = Boundary(label);
  tde = DistTrans(bndr,label);
  objmap = SigmoidalStretch(tde->V, 5, -7, 0.35);
  objgrad = ObjectGradient(objmap,1.5);
  
  /* computing image gradient */
  feat = MorphImageFeats(img,2);
  imggrad = FeaturesGradient(feat, 1.5);

  /* combining gradients */
  grad = CombineGradients(objgrad, imggrad, wobj);

  t2 = Toc();

  fprintf(stdout,"Gradient computing in %f ms\n",CTime(t1,t2));
  WriteImage(bndr,"boundary.pgm");
  WriteImage(tde->V,"tde.pgm");
  WriteImage(objmap,"objmap.pgm");
  WriteImage(objgrad,"objgrad.pgm");
  WriteImage(imggrad,"imggrad.pgm");
  WriteImage(grad,"grad.pgm");
  
  DestroyForest(&tde);
  DestroyImage(&bndr);
  DestroyImage(&img);
  DestroyImage(&label);
  DestroyImage(&objmap);
  DestroyImage(&imggrad);
  DestroyImage(&objgrad);
  DestroyImage(&grad);
  DestroyFeatures(&feat);

  /* The following block must the remarked when using non-linux machines */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   
  

  return(0);
}
