#include "ift.h"

#define MAXGRAD 255

/* Object map computation */

// Extracts the boundary of an image
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

// Signed euclidean distance transform
// For pixels where the label I->val[p] == 0
// the EDT value will be negative and positive otherwise
DImage *SignedDistTrans(Image *B, Image *I) 
{
  int p,q,n=B->ncols*B->nrows,i,tmp;
  Pixel u,v,w;
  AdjRel *A=Circular(1.5),*A4=Circular(1.0);
  Image* P = CreateImage(I->ncols,I->nrows);
  Image* R = CreateImage(I->ncols,I->nrows);
  Image* V = CreateImage(I->ncols,I->nrows);
  DImage* edt = CreateDImage(I->ncols,I->nrows);

  GQueue *Q=CreateGQueue(1024,n,V->val);

  // Trivial path initialization

  for (p=0; p < n; p++) {
    u.x = p % B->ncols;
    u.y = p / B->ncols;
    V->val[p]=INT_MAX; R->val[p]=p; P->val[p]=NIL;
    if (B->val[p]!=0){ // p belongs to an object's border
      V->val[p]=0;
      InsertGQueue(&Q,p);
    }
  }

  // Path propagation

  while(!EmptyGQueue(Q)){
    p = RemoveGQueue(Q);
    u.x = p % B->ncols;
    u.y = p / B->ncols;
    w.x = R->val[p] % B->ncols;
    w.y = R->val[p] / B->ncols;
    for (i=1; i < A->n; i++) {
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if (ValidPixel(B,v.x,v.y)){
	q   = v.x + B->tbrow[v.y];
	if (V->val[q]>V->val[p]){	    
	  tmp = (v.x-w.x)*(v.x-w.x)+(v.y-w.y)*(v.y-w.y);
	  if (tmp < V->val[q]){
	    if (V->val[q]!=INT_MAX) RemoveGQueueElem(Q, q);
	    V->val[q]=tmp; R->val[q]=R->val[p]; P->val[q]=p;
	    InsertGQueue(&Q,q);
	  }
	}
      }
    }
    edt->val[p] = sqrt((double)V->val[p]);
    if(I->val[p] == 0) edt->val[p] = -edt->val[p];
  }

  DestroyGQueue(&Q);
  DestroyAdjRel(&A);
  DestroyAdjRel(&A4);
  DestroyImage(&P);
  DestroyImage(&V);
  DestroyImage(&R);

  return edt;
}

// Sigmoidal rescaling of values. Max, min, beta and alpha must be in the same range
DImage* SigmoidalStretch(DImage* img, double thresh, double min, double max, double beta, double alpha)
{
  int p;
  DImage* result = CreateDImage(img->ncols, img->nrows);

  for(p = 0; p < img->ncols*img->nrows; p++)
    {
      if(img->val[p] < -thresh) 
    	result->val[p] = min;
      else if(img->val[p] > thresh) 
	result->val[p] = max;
      else
	result->val[p] = (max - min)/(1.0 + exp(-(img->val[p] - beta)/alpha)) + min;
    }
  return result;
}

/* Gradient computation */

// Object-based gradient computation
Image *ObjectGradient(DImage *img, float radius)
{
    real    dist,gx,gy;
    int     i,p,q,n=img->ncols*img->nrows;
    Pixel   u,v;
    AdjRel *A=Circular(radius);
    real   *md=AllocRealArray(A->n);

    Image* grad = CreateImage(img->ncols, img->nrows);

    double Imax = MaximumDImageValue(img);

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
            if (ValidDImagePixel(img,v.x,v.y))
            {
                q    = v.x + img->tbrow[v.y];
                dist = ((float)img->val[q]-(float)img->val[p])/Imax;

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

// Image-based gradient computation
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

// Linear combination of gradients
Image* CombineGradients(Image *objgrad, Image *imggrad, float wobj)
{
  int p;
  Image *grad = CreateImage(objgrad->ncols, objgrad->nrows);

  for(p = 0; p < objgrad->ncols*objgrad->nrows; p++)
    grad->val[p] = wobj*objgrad->val[p] + (1-wobj)*imggrad->val[p];

  return grad;
}

int main(int argc, char **argv)
{
  float wobj = 0.5;
  char outfile[100];
  char *file_noext;
  timer    *t1=NULL,*t2=NULL;
  Image    *bndr=NULL, *label=NULL, *ch_label=NULL, *final_label=NULL;
  Image    *objgrad=NULL, *imggrad=NULL, *grad=NULL, *tmp=NULL;
  Features *feat=NULL;
  DImage   *edt=NULL, *objmap=NULL;
  Subgraph *sg=NULL, *sgtrain=NULL, *sgeval=NULL;
  Set      *Obj=NULL,*Bkg=NULL;
  AdjRel   *A=NULL;
 
  /* The following block must the remarked when using non-linux machines */

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  
  /*----------------------------------------------------------------------*/
  
  if (argc!=4){
    fprintf(stderr,"Usage: %s <image.pgm (.ppm)>  <seeds.txt> <wobj [0,1]>\n",argv[0]);
    fprintf(stderr,"image.pgm (.ppm): image to be classified\n");
    fprintf(stderr,"seeds.txt: seed pixels\n");
    fprintf(stderr,"wobj: object gradient weight\n");
    exit(-1);
  }

  char *ext = strrchr(argv[1],'.');

  if(!strcmp(ext,".pgm"))
  {
    Image   *img=NULL;
    img   = ReadImage(argv[1]);
    feat  = GaussImageFeats(img, 2);
    DestroyImage(&img);  
  }else{
    CImage   *cimg=NULL;
    cimg   = ReadCImage(argv[1]);
    feat   = GaussCImageFeats(cimg, 2);
    DestroyCImage(&cimg);
  }

  file_noext = strtok(argv[1],".");
  A = Circular(2);
  ReadSeeds(argv[2],&Obj,&Bkg);

  wobj = atof(argv[3]);

  t1 = Tic();

  sg = SubgraphFromSeeds(feat,Obj,Bkg);
  SplitSubgraph(sg, &sgtrain, &sgeval, 0.2);
  
  /* supervised binary classification using the Optimum-Path Forest classifier */
  
  OPFLearning(&sgtrain, &sgeval);
  label = OPFClassifyImage(sgtrain, feat);
  
  //reducing label noise
  ch_label = CloseHoles(label);
  final_label = OpenRec(ch_label,A);
  
  /* computing object gradient */
  bndr = Boundary(final_label);
  edt = SignedDistTrans(bndr,final_label);
  
  //thresh = 10.0, max = 1, min = 0.0, beta = 0.0 (EDT has value 0.0 on the boundary), alpha = 1.0
  objmap = SigmoidalStretch(edt, 100.0, 0.0, 1.0, 0.0, 1.0); 
  objgrad = ObjectGradient(objmap,1.5);
  
  /* computing image gradient */
  imggrad = FeaturesGradient(feat, 1.5);

  /* combining gradients */
  grad = CombineGradients(objgrad, imggrad, wobj);

  t2 = Toc();

  fprintf(stdout,"Gradient computing in %f ms\n",CTime(t1,t2));
  
  tmp = ConvertDImage2Image(edt);
  sprintf(outfile,"%s_signed_edt.pgm",file_noext);
  WriteImage(tmp,outfile);
  DestroyImage(&tmp);
  
  tmp = ConvertDImage2Image(objmap);
  sprintf(outfile,"%s_objmap.pgm",file_noext);
  WriteImage(tmp,outfile);
  DestroyImage(&tmp);
  
  sprintf(outfile,"%s_objgrad.pgm",file_noext);
  WriteImage(objgrad,outfile);
  sprintf(outfile,"%s_imggrad.pgm",file_noext);
  WriteImage(imggrad,outfile);
  sprintf(outfile,"%s_grad.pgm",file_noext);
  WriteImage(grad,outfile);
  
  DestroyDImage(&edt);
  DestroyImage(&bndr);
  DestroyImage(&label);
  DestroyImage(&ch_label);
  DestroyImage(&final_label);
  DestroyDImage(&objmap);
  DestroyImage(&imggrad);
  DestroyImage(&objgrad);
  DestroyImage(&grad);
  DestroyFeatures(&feat);
  DestroySubgraph(&sg);
  DestroySubgraph(&sgtrain);
  DestroySubgraph(&sgeval);
  DestroySet(&Obj);
  DestroySet(&Bkg);
  DestroyAdjRel(&A);

  /* The following block must the remarked when using non-linux machines */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   
  

  return(0);
}
