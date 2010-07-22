#include "ift.h"

/* Explanation for some variables

V is the cost or connectivity map
P is the predecessor map
I is the mask image 
J is the marker image

*/

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

Set *ImageBorder(Image *img)
{
  Set *S=NULL;
  Pixel u;

  for (u.x=0; u.x < img->ncols; u.x++){ 
    InsertSet(&S,u.x);
    InsertSet(&S,u.x+img->tbrow[img->nrows-1]);
  }
  for (u.y=1; u.y < img->nrows-1; u.y++){ 
    InsertSet(&S,img->tbrow[u.y]);
    InsertSet(&S,img->ncols-1+img->tbrow[u.y]);
  }
  return(S);
}

// Superior reconstruction img <= marker without marker imposition

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


// Superior reconstruction img <= marker with marker imposition

Image *SupRecMI(Image *img, Set **S)
{
  Image *cost=NULL;
  GQueue *Q=NULL;
  int i,p,q,tmp,n;
  Pixel u,v;
  AdjRel *A=Circular(1.0);

  n     = img->ncols*img->nrows;
  cost  = CreateImage(img->ncols,img->nrows);
  Q     = CreateGQueue(MaximumValue(img)+1,n,cost->val);

  // Trivial path initialization

  for (p=0; p < n; p++) {
    cost->val[p]=INT_MAX; 
  }

  while(*S!=NULL){
    p = RemoveSet(S);
    cost->val[p]=img->val[p]; 
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
	    if (cost->val[q]!=INT_MAX) RemoveGQueueElem(Q,q);
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

// Inferior Reconstruction img >= marker with marker imposition

Image *InfRecMI(Image *img, Set **S)
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
    cost->val[p]=INT_MIN; 
  }

  while(*S!=NULL){
    p = RemoveSet(S);
    cost->val[p]=img->val[p]; 
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
	    if (cost->val[q]!=INT_MIN) RemoveGQueueElem(Q,q);
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

Image *CloseBasins(Image *img)
{
  Image *cbas;
  Set *S=NULL;
  
  S    = ImageBorder(img);
  cbas = SupRecMI(img,&S);

  return(cbas);
}

Image *OpenDomes(Image *img)
{
  Image *odom;
  Set *S=NULL;
  
  S    = ImageBorder(img);
  odom = InfRecMI(img,&S);

  return(odom);
}

Image *FastDilate(Image *I,Set **S, float radius)
{
  Image  *J=CreateImage(I->ncols,I->nrows);
  AdjRel *A=Circular(sqrt(2)),*A4=Circular(1.0);
  int    nspels,*V=AllocIntArray(nspels=I->ncols*I->nrows);
  int    *R=AllocIntArray(nspels);
  GQueue *Q=CreateGQueue(1024,nspels,V);
  int    p,q,i,tmp;
  float  r2=(radius*radius);
  Pixel  u,v,w;
  char   first_time=0;

  // Initialization

  if (*S==NULL)
    first_time=1;

  for(p=0; p < nspels; p++) {
    J->val[p]=I->val[p]; R[p]=p; 
    if (I->val[p]==0) // p is background
      V[p]=INT_MAX;
    else{ // p is object
      V[p]=INT_MIN;
      if (first_time){ // detect if p is border
	u.x = p % I->ncols;
	u.y = p / I->ncols;
	for (i=1; i < A4->n; i++) {
	  v.x = u.x + A4->dx[i];
	  v.y = u.y + A4->dy[i];
	  if (ValidPixel(I,v.x,v.y)){
	    q = v.x + I->tbrow[v.y];
	    if (I->val[q]==0){ // insert p in the border set S
	      InsertSet(S,p);
	      break;
	    }
	  }
	}
      }
    }
  }

  while (*S != NULL){
    p   = RemoveSet(S); 
    V[p]= 0; 
    InsertGQueue(&Q,p);
  }


  // Euclidean IFT dilation

  while(!EmptyGQueue(Q)){
    p = RemoveGQueue(Q);
    if (V[p] <= r2) {
      J->val[p]=1; 
      u.x = p % I->ncols;
      u.y = p / I->ncols;
      w.x = R[p] % I->ncols;
      w.y = R[p] / I->ncols;
      for (i=1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(I,v.x,v.y)){
	  q   = v.x + I->tbrow[v.y];
	  if (V[q]>V[p]){	    
	    tmp = (v.x-w.x)*(v.x-w.x)+(v.y-w.y)*(v.y-w.y);
	    if (tmp < V[q]){
	      if (V[q]!=INT_MAX){
		RemoveGQueueElem(Q, q);
	      }
	      V[q]=tmp; R[q]=R[p];
	      InsertGQueue(&Q,q);
	    }
	  }
	}
      }
    }else{
      InsertSet(S,p);
    }
  }

  free(V); free(R);
  DestroyGQueue(&Q);
  DestroyAdjRel(&A);
  DestroyAdjRel(&A4);

  return(J);
}

Image *FastErode(Image *I,Set **S, float radius)
{
  Image  *J=CreateImage(I->ncols,I->nrows);
  AdjRel *A=Circular(sqrt(2)),*A4=Circular(1.0);
  int    nspels,*V=AllocIntArray(nspels=I->ncols*I->nrows);
  int    *R=AllocIntArray(nspels);
  GQueue *Q=CreateGQueue(1024,nspels,V);
  int    p,q,i,tmp;
  float  r2=(radius*radius);
  Pixel  u,v,w;
  char   first_time=0;

  // Initialization

  if (*S==NULL)
    first_time=1;

  for(p=0; p < nspels; p++) {
    J->val[p]=I->val[p]; R[p]=p; 
    if (I->val[p]==1) // p is object
      V[p]=INT_MAX;
    else{ // p is background
      V[p]=INT_MIN;
      if (first_time){ // detect if p is border
	u.x = p % I->ncols;
	u.y = p / I->ncols;
	for (i=1; i < A4->n; i++) {
	  v.x = u.x + A4->dx[i];
	  v.y = u.y + A4->dy[i];
	  if (ValidPixel(I,v.x,v.y)){
	    q = v.x + I->tbrow[v.y];
	    if (I->val[q]==1){ // insert p in the border set S
	      InsertSet(S,p);
	      break;
	    }
	  }
	}
      }
    }
  }

  while (*S != NULL){
    p   = RemoveSet(S); 
    V[p]= 0; 
    InsertGQueue(&Q,p);
  }


  // Euclidean IFT erosion

  while(!EmptyGQueue(Q)){
    p = RemoveGQueue(Q);
    if (V[p] <= r2) {
      J->val[p]=0; 
      u.x = p % I->ncols;
      u.y = p / I->ncols;
      w.x = R[p] % I->ncols;
      w.y = R[p] / I->ncols;
      for (i=1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(I,v.x,v.y)){
	  q   = v.x + I->tbrow[v.y];
	  if (V[q]>V[p]){	    
	    tmp = (v.x-w.x)*(v.x-w.x)+(v.y-w.y)*(v.y-w.y);
	    if (tmp < V[q]){
	      if (V[q]!=INT_MAX){
		RemoveGQueueElem(Q, q);
	      }
	      V[q]=tmp; J->val[q]=0; R[q]=R[p];
	      InsertGQueue(&Q,q);
	    }
	  }
	}
      }
    }else{
      InsertSet(S,p);
    }
  }

  free(V); free(R);
  DestroyGQueue(&Q);
  DestroyAdjRel(&A);
  DestroyAdjRel(&A4);
  return(J);
}

Image *FastClose(Image *I, float radius)
{
  Image *J,*K;
  Set *S=NULL;
  J = FastDilate(I,&S,radius);
  K = FastErode(J,&S,radius);
  DestroyImage(&J);
  DestroySet(&S);
  return(K);
}

Image *FastOpen(Image *I, float radius)
{
  Image *J,*K;
  Set *S=NULL;
  J = FastErode(I,&S,radius);
  K = FastDilate(J,&S,radius);
  DestroyImage(&J);
  DestroySet(&S);
  return(K);
}

Image *FastCloseRec(Image *I, float radius)
{
  Image *J,*K;

  J = FastClose(I,radius);
  K = SupRec(I,J);
  DestroyImage(&J);
  return(K);
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

Image *FastOpenRec(Image *I, float radius)
{
  Image *J,*K;
  J = FastOpen(I,radius);
  K = InfRec(I,J);
  DestroyImage(&J);
  return(K);
}


int main(int argc, char **argv)
{
  Image  *img[8];
  Set *S=NULL;
  AdjRel *A=NULL;
  int i;
  timer    *t1=NULL,*t2=NULL;

  /*------- -------------------------------------------------*/
#ifndef _WIN32
  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
#endif

  /*--------------------------------------------------------*/
  
  if (argc != 2) {
    fprintf(stderr,"usage: morphfilters <orig.pgm> \n");
    exit(-1);
  }

  img[0]= ReadImage(argv[1]);
  
  img[1]= CloseBasins(img[0]);
  WriteImage(img[1],"closebasins.pgm");

  img[2]= OpenDomes(img[0]);
  WriteImage(img[2],"opendomes.pgm");

  A = Circular(10.0);

  t1 = Tic();
  
  img[3]= CloseRec(img[0],A);

  t2 = Toc();    
  
  fprintf(stdout,"Closing by reconstruction in %f ms\n",CTime(t1,t2));
  

  WriteImage(img[3],"closerec.pgm");

  img[4]= OpenRec(img[0],A);
  WriteImage(img[4],"openrec.pgm");

  DestroyAdjRel(&A);
  
  for (i=1; i < 5; i++) 
    DestroyImage(&img[i]);
  
  if (MaximumValue(img[0])==1){ // Apply fast binary filters

    img[1] = AddFrame(img[0],10,0); // to guarantee space for dilation
    WriteImage(img[1],"shape.pgm");

    img[2]=FastDilate(img[1],&S,10.0);
    WriteImage(img[2],"fastdilation.pgm");
    DestroySet(&S);
  
    img[3]=FastErode(img[1],&S,10.0);
    WriteImage(img[3],"fasterosion.pgm");
    DestroySet(&S);
    
    img[4]=FastClose(img[1],10.0);
    WriteImage(img[4],"fastclose.pgm");
    
    img[5]= FastOpen(img[1],10.0);
    WriteImage(img[5],"fastopen.pgm");
    
    t1 = Tic();

    img[6]=FastCloseRec(img[1],10.0);

    t2 = Toc();    
  
    fprintf(stdout,"Fast closing by reconstruction in %f ms\n",CTime(t1,t2));

    WriteImage(img[6],"fastcloserec.pgm");
    
    img[7]=FastOpenRec(img[1],10.0);
    WriteImage(img[7],"fastopenrec.pgm");
  
    for (i=1; i < 8; i++) 
      DestroyImage(&img[i]);
  }

  DestroyImage(&img[0]);

  /* ---------------------------------------------------------- */

#ifndef _WIN32
  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   
#endif

  return(0);
}
