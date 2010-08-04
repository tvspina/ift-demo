/*
    Copyright (C) <2010> <Alexandre Xavier FalcÃ£o and Thiago Vallin Spina>

    This file is part of IFT-demo.

    IFT-demo is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IFT-demo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with IFT-demo.  If not, see <http://www.gnu.org/licenses/>.

    please see full copyright in COPYING file.
    -------------------------------------------------------------------------

    written by A.X. Falcão <afalcao@ic.unicamp.br> and by T.V. Spina
    <tvspina@liv.ic.unicamp.br>, 2010

*/

#include "ift.h"

/* Papers related to this program:


@inproceedings{Spina09a,
author = "T.V. Spina and J.A. Montoya-Zegarra and A.X. Falc{\~{a}}o and P.A.V. Miranda",
title = "Fast interactive segmentation of natural images using the Image Foresting Transform",
booktitle = "Proc. of the 16th Intl. Conf. on Digital Signal Processing",
address = "Santorini, Greece",
publisher = "IEEE",
doi = "10.1109/ICDSP.2009.5201044",
year = 2009
}

@article{Papa09a,
author = "J.P. Papa and A.X. Falc{\~{a}}o and C.T.N. Suzuki",
title = "Supervised Pattern Classification based on Optimum-Path
Forest", 
journal = "Intl. Journal of Imaging Systems and Technology",
publisher = "Wiley",
doi = "10.1002/ima.20188",
volume = "19",
number = "2",
pages = "120--131",
month = "Jun",
year = 2009
}

@article{Miranda10a, 
author = "P.A.V. Miranda and A.X. Falc{\~{a}}o and J.K. Udupa",
title  = "Synergistic Arc-Weight Estimation for Interactive Image Segmentation using Graphs",
journal = "Computer Vision and Image Understanding",
publisher = "Elsevier",
doi = "10.1016/j.cviu.2009.08.001", 
volume = 114, 
number = 1, 
pages = "85--99",
month = "Jan",
year = 2010
}

*/

#define MAXGRAD 255

/* Object map computation */

// Computes the optimum path forest
// on the graph *sg. For fuzzy classification
// we have two separate complete graphs, one for the
// object pixels and another one for the background pixels,
// therefore this function is called twice (see FuzzyOPFLearning)

void FuzzyOPFTraining(Subgraph *sg)
{
    int p;
    RealHeap *Q = NULL;
    float *pathval = NULL;

    // initialization

    pathval = AllocFloatArray(sg->nnodes);

    Q=CreateRealHeap(sg->nnodes, pathval);

    for (p = 0; p < sg->nnodes; p++)
    {
        if (sg->node[p].status==PROTOTYPE)
        {
            pathval[p]         = 0.0;
            sg->node[p].label  = sg->node[p].truelabel;
            sg->node[p].pred   = NIL;
            InsertRealHeap(Q, p);
        }
        else  // non-prototypes
        {
            pathval[p]  = FLT_MAX;
        }
    }

    // IFT with fmax

    int i = 0;

    while ( !IsEmptyRealHeap(Q) )
    {
        RemoveRealHeap(Q, &p);

        sg->ordered_list_of_nodes[i]=p;
        sg->node[p].pathval = pathval[p];
        i++;

        int q;
        for (q=0; q < sg->nnodes; q++)
        {
            if (p!=q)
            {
                if (pathval[p] < pathval[q])
                {
                    float weight;
		    weight = EuclDistLog(sg->node[p].feat,sg->node[q].feat,sg->nfeats);


                    float tmp  = MAX(pathval[p],weight);
                    if ( tmp < pathval[ q ] )
                    {
                        sg->node[q].pred  = p;
                        sg->node[q].label = sg->node[p].label;
                        UpdateRealHeap(Q, q, tmp);
                    }
                }
            }
        }
    }

    DestroyRealHeap( &Q );
    free( pathval );
}


// Classifies nodes of evaluation/test set
void FuzzyOPFClassify(Subgraph *sgtrainobj, Subgraph* sgtrainbkg, Subgraph *sg)
{
    register int i, j;
    int objlabel = sgtrainobj->node[0].truelabel;
    int bkglabel = sgtrainbkg->node[0].truelabel;

    for (i = 0; i < sg->nnodes; i++)
    {
        float minCostObj;
        minCostObj = FLT_MAX;
        register float weight, cost;
        for (j = 0; j < sgtrainobj->nnodes; j++)
        {

            weight = EuclDistLog(sgtrainobj->node[j].feat,sg->node[i].feat,sg->nfeats);

            cost = MAX(sgtrainobj->node[j].pathval, weight);

            if (cost < minCostObj)
            {
                minCostObj = cost;
            }
        }

        float minCostBkg = FLT_MAX;
        for (j = 0; j < sgtrainbkg->nnodes; j++)
        {
            weight = EuclDistLog(sgtrainbkg->node[j].feat,sg->node[i].feat,sg->nfeats);

            cost = MAX(sgtrainbkg->node[j].pathval, weight);

            if (cost < minCostBkg)
            {
                minCostBkg = cost;
            }
        }

        if (minCostObj == minCostBkg)
        {
            /// Forcing error
            if (sg->node[i].truelabel == objlabel)
            {
                sg->node[i].label = bkglabel;
            }
            else
            {
                sg->node[i].label = objlabel;
            }
        }
        else
        {
            if (minCostObj < minCostBkg)
            {
                sg->node[i].label = objlabel;
            }
            else
            {
                sg->node[i].label = bkglabel;
            }
        }
    }
}


// Executes the OPF learning procedure to replace 
// misclassified samples in the evaluation set by non prototypes from
// the training set
void FuzzyOPFLearning(Subgraph* sg, Subgraph** sgtrainobj, Subgraph** sgtrainbkg, float perc)
{
    int i;
    const int iterations = 5;
    float Acc,MaxAcc=FLT_MIN;
    Subgraph *sgtrainobj1 = NULL;
    Subgraph *sgtrainbkg1 = NULL;
    Subgraph *sgtrain = NULL, *sgeval = NULL;

    SplitSubgraph(sg,&sgtrain,&sgeval,perc);

    for (i = 1; i <= iterations && MaxAcc != 1.0; i++)
    {
        fprintf(stdout, "\nrunning iteration ... %d ", i);

        MSTPrototypes(sgtrain);

        *sgtrainobj = SplitSubgraphByTrueLabel(sgtrain, OPF_OBJ_LABEL);
        *sgtrainbkg = SplitSubgraphByTrueLabel(sgtrain, OPF_BKG_LABEL);

        FuzzyOPFTraining(*sgtrainobj);

        FuzzyOPFTraining(*sgtrainbkg);

        FuzzyOPFClassify(*sgtrainobj, *sgtrainbkg, sgeval);
        Acc = Accuracy(sgeval);
        if (Acc > MaxAcc)
        {
            MaxAcc = Acc;
            if (sgtrainobj1!=NULL) DestroySubgraph(&sgtrainobj1);
            if (sgtrainbkg1!=NULL) DestroySubgraph(&sgtrainbkg1);
            sgtrainobj1 = CopySubgraph(*sgtrainobj);
            sgtrainbkg1 = CopySubgraph(*sgtrainbkg);
        }
        SwapErrorsbyNonPrototypes(&sgtrain, &sgeval);

	fprintf(stdout,"Acc: %f\n", Acc);

        DestroySubgraph(sgtrainobj);
        DestroySubgraph(sgtrainbkg);
    }

    *sgtrainobj = sgtrainobj1;
    *sgtrainbkg = sgtrainbkg1;

    fprintf(stderr,"Best accuracy %f\n",MaxAcc);

    DestroySubgraph(&sgeval);
    DestroySubgraph(&sgtrain);
}

// Computes the optimum path cost for every pixel in *f
// from the classifier *sg (object or background forest)
DImage *FuzzyOPFPathCostMap(Subgraph *sg, Features *f)
{
    DImage *pvalmap=CreateDImage(f->ncols,f->nrows);

    register int p,q,n=f->nelems;

    for (q=n; q--;)
    {
        register float mincost = FLT_MAX;
        register int k;
        p = 0;
        do
        {
            k = sg->ordered_list_of_nodes[p];

	    register float weight = EuclDistLog(f->elem[q].feat,sg->node[k].feat,sg->nfeats);
	    register float tmp  = MAX(sg->node[k].pathval,weight);

	    if (tmp < mincost)
	      mincost = tmp;
	    p++;
        }
        while (p < sg->nnodes -1 &&
                mincost > sg->node[sg->ordered_list_of_nodes[k+1]].pathval);

        pvalmap->val[q] = (double)mincost;
    }
    return(pvalmap);
}

// Combines the computed path cost images *d1 and *d2
// to obtain the final result (d1 == object path cost
// d2 == background path cost, because we expect that
// the cost to the background forest of pixels that resemble
// the object be greater than the cost to the object forest)
DImage *MembershipMap(DImage *d1, DImage *d2)
{
    DImage *map;
    register int p,n;

    map = CreateDImage(d1->ncols,d1->nrows);
    n = d1->ncols*d1->nrows;
    for (p=0; p<n; p++)
    {
        register float fd1,fd2;

        fd1 = (float)d1->val[p];
        fd2 = (float)d2->val[p];
        if ((fd1+fd2)<0.00001 || fd1 == fd2)
        {
            map->val[p] = 0.5;
        }
        else
        {
            map->val[p] = fd2/(fd2+fd1);
        }
    }
    return map;
}

// Uses the object and background forests *sgrainobj and *sgtrainbkg, respectively,
// to compute an object membership map using the set of feature vectors *f
DImage* FuzzyOPFObjectMembershipMap(Subgraph *sgtrainobj, Subgraph *sgtrainbkg, Features *f)
{
  DImage *obj_pv = FuzzyOPFPathCostMap(sgtrainobj, f);
  DImage *bkg_pv = FuzzyOPFPathCostMap(sgtrainbkg, f);

  DImage* objMap = MembershipMap(obj_pv, bkg_pv);

  DestroyDImage(&obj_pv);
  DestroyDImage(&bkg_pv);

  return objMap;
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
  Image    *objgrad=NULL, *imggrad=NULL, *grad=NULL, *tmp=NULL;
  Features *feat=NULL;
  DImage   *objmap=NULL;
  Subgraph *sg=NULL, *sgtrainobj=NULL, *sgtrainbkg=NULL, *sgeval=NULL;
  Set      *Obj=NULL,*Bkg=NULL;


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

    Features *gaussfeats   = GaussCImageFeats(cimg, 2);

    int p,j;
    //converting features from [0,1] to [0,255]
    for(p = 0; p < gaussfeats->nelems; p++)
      for(j = 0; j < gaussfeats->nfeats; j++)
	gaussfeats->elem[p].feat[j] *= gaussfeats->Imax;

    feat = LabFeats(gaussfeats);
    //converting features from [0,255] to [0,1]
    for(p = 0; p < feat->nelems; p++)
      for(j = 0; j < feat->nfeats; j++)
	feat->elem[p].feat[j] /= feat->Imax;

    DestroyCImage(&cimg);
    DestroyFeatures(&gaussfeats);
  }

  file_noext = strtok(argv[1],".");

  ReadSeeds(argv[2],&Obj,&Bkg);

  wobj = atof(argv[3]);

  t1 = Tic();

  sg = SubgraphFromSeeds(feat,Obj,Bkg);

  /* supervised fuzzy classification using the Optimum-Path Forest classifier */
  FuzzyOPFLearning(sg, &sgtrainobj, &sgtrainbkg, 0.2);
  objmap = FuzzyOPFObjectMembershipMap(sgtrainobj, sgtrainbkg, feat);

  /* computing object gradient */
  objgrad = ObjectGradient(objmap,1.5);

  /* computing feature vector gradient */
  imggrad = FeaturesGradient(feat, 1.5);

  /* combining gradients */
  grad = CombineGradients(objgrad, imggrad, wobj);

  t2 = Toc();

  fprintf(stdout,"Gradient computing in %f ms\n",CTime(t1,t2));

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

  DestroyDImage(&objmap);
  DestroyImage(&imggrad);
  DestroyImage(&objgrad);
  DestroyImage(&grad);
  DestroyFeatures(&feat);
  DestroySubgraph(&sg);
  DestroySubgraph(&sgtrainobj);
  DestroySubgraph(&sgtrainbkg);
  DestroySubgraph(&sgeval);
  DestroySet(&Obj);
  DestroySet(&Bkg);

  /* The following block must the remarked when using non-linux machines */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);


  return(0);
}
