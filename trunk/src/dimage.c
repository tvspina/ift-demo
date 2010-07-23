#include "dimage.h"
#include "common.h"

DImage *CreateDImage(int ncols, int nrows)
{
  DImage *dimg=NULL;
  int i;

  dimg = (DImage *) calloc(1,sizeof(DImage));
  if (dimg == NULL){
    Error(MSG1,"CreateDImage");
  }

  dimg->val   = AllocDoubleArray(nrows*ncols);
  dimg->tbrow = AllocIntArray(nrows);

  dimg->tbrow[0]=0;
  for (i=1; i < nrows; i++)
    dimg->tbrow[i]=dimg->tbrow[i-1]+ncols;
  dimg->ncols = ncols;
  dimg->nrows = nrows;

 return(dimg);
}

void DestroyDImage(DImage **dimg)
{
  DImage *aux;

  aux = *dimg;
  if(aux != NULL){
    if (aux->val != NULL)   free(aux->val);
    if (aux->tbrow != NULL) free(aux->tbrow);
    free(aux);
    *dimg = NULL;
  }
}

DImage  *ReadDImage(char *filename){
  DImage *dimg;
  char msg[512];
  int  ncols,nrows,n;
  FILE *fp;

  fp = fopen(filename,"rb");
  if (fp == NULL){
    sprintf(msg,"Cannot open %s",filename);
    Error(msg,"ReadDImage");
  }

  fscanf(fp,"%d %d\n",&ncols,&nrows);
  dimg = CreateDImage(ncols, nrows);
  n = ncols*nrows;
  fread(dimg->val, sizeof(double), n, fp);

  return dimg;
}

void    WriteDImage(DImage *dimg, char *filename){
  char msg[512];
  int n;
  FILE *fp;

  fp = fopen(filename,"wb");
  if (fp == NULL){
    sprintf(msg,"Cannot open %s",filename);
    Error(msg,"WriteDImage");
  }

  fprintf(fp,"%d %d\n",dimg->ncols,dimg->nrows);
  n = dimg->ncols*dimg->nrows;
  fwrite(dimg->val, sizeof(double), n ,fp);
}

bool ValidDImagePixel(DImage *dimg, int x, int y)
{
  if ((x >= 0)&&(x < dimg->ncols)&&
      (y >= 0)&&(y < dimg->nrows))
    return(true);
  else
    return(false);
}

double  MaximumDImageValue(DImage *dimg)
{
  double max;
  int i,n;

  n = dimg->ncols*dimg->nrows;
  max = dimg->val[0];
  for (i=1; i < n; i++)
    if (dimg->val[i] > max)
      max = dimg->val[i];

  return(max);
}

double  MinimumDImageValue(DImage *dimg)
{
  double min;
  int i,n;

  n = dimg->ncols*dimg->nrows;
  min = dimg->val[0];
  for (i=1; i < n; i++)
    if (dimg->val[i] < min)
      min = dimg->val[i];

  return(min);
}

void SetDImage(DImage *dimg, double value)
{
  int i,n;
  n = dimg->ncols*dimg->nrows;
  for (i=0; i < n; i++){
    dimg->val[i]=value;
  }
}

DImage *CopyDImage(DImage *dimg)
{
  DImage *dimgc;

  dimgc = CreateDImage(dimg->ncols,dimg->nrows);
  memcpy(dimgc->val,dimg->val,dimg->ncols*dimg->nrows*sizeof(double));

  return(dimgc);
}


Image *ConvertDImage2Image(DImage *dimg){
  Image *img;
  double max,min;
  int p,n;

  n = dimg->ncols*dimg->nrows;
  img = CreateImage(dimg->ncols, dimg->nrows);
  max = MaximumDImageValue(dimg);
  min = MinimumDImageValue(dimg);
  if (max==min) {
    for(p=0; p<n; p++)
      img->val[p] = max;
    return img;
  }

  for(p=0; p<n; p++){
    img->val[p] = ROUND(255.0*(dimg->val[p]-min)/(max-min));
  }

  return img;
}

