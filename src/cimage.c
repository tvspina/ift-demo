#include "cimage.h"

CImage *CreateCImage(int ncols, int nrows)
{
  CImage *cimg=NULL;
  int i;

  cimg = (CImage *) calloc(1, sizeof(CImage));
  for (i=0; i < 3; i++)
    cimg->C[i] = CreateImage(ncols,nrows);
  return(cimg);
}

void    DestroyCImage(CImage **cimg)
{
  CImage *tmp;
  int i;

  tmp = *cimg;
  if (tmp != NULL) {
    for (i=0; i < 3; i++)
      DestroyImage(&(tmp->C[i]));
    free(tmp);
    *cimg = NULL;
  }
}

CImage *ReadCImage(char *filename)
{
  CImage *cimg=NULL;
  FILE *fp=NULL;
  char type[10];
  int  i,ncols,nrows,n;
  char z[256];

  fp = fopen(filename,"rb");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  if(fscanf(fp,"%s\n",type)!=0);
  if((strcmp(type,"P6")==0)){
    NCFgets(z,255,fp);
    sscanf(z,"%d %d\n",&ncols,&nrows);
    n = ncols*nrows;
    NCFgets(z,255,fp);
    sscanf(z,"%d\n",&i);
    cimg = CreateCImage(ncols,nrows);
    for (i=0; i < n; i++){
      cimg->C[0]->val[i] = fgetc(fp);
      cimg->C[1]->val[i] = fgetc(fp);
      cimg->C[2]->val[i] = fgetc(fp);
    }
    fclose(fp);
  }else{
    fprintf(stderr,"Input image must be P6\n");
    exit(-1);
  }

  return(cimg);
}

void    WriteCImage(CImage *cimg, char *filename)
{
  FILE *fp;
  int i,n;

  fp = fopen(filename,"w");
  fprintf(fp,"P6\n");
  fprintf(fp,"%d %d\n",cimg->C[0]->ncols,cimg->C[0]->nrows);
  fprintf(fp,"255\n");
  n = cimg->C[0]->ncols*cimg->C[0]->nrows;
  for (i=0; i < n; i++) {
    fputc(cimg->C[0]->val[i],fp);
    fputc(cimg->C[1]->val[i],fp);
    fputc(cimg->C[2]->val[i],fp);
  }
  fclose(fp);
}


CImage *DrawLabeledRegions(Image *img, Image *label){
  CImage *border=CreateCImage(img->ncols,img->nrows);
  int x,y,k,p,q,u,v;
  AdjRel *A;

  A = Circular(1.0);
  for(y=0;y<img->nrows;y++){
    for(x=0;x<img->ncols;x++){
      p = x + img->tbrow[y];
      border->C[0]->val[p]=
	border->C[1]->val[p]=border->C[2]->val[p]=img->val[p];

      for(k=1;k<A->n;k++){
	u = x+A->dx[k];
	v = y+A->dy[k];
	if(ValidPixel(img,u,v)){
	  q= u + img->tbrow[v];
	  if (label->val[p] != label->val[q]){
	    border->C[0]->val[p]=255;
	    border->C[1]->val[p]=0;
	    border->C[2]->val[p]=0;
	    break;
	  }
	}
      }
    }
  }
  DestroyAdjRel(&A);

  return(border);
}
