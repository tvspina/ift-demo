/*
    Copyright (C) <2010> <Alexandre Xavier Falcão and Thiago Vallin Spina>

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

int main(int argc, char **argv)
{
  timer    *t1=NULL,*t2=NULL;
  Image    *label=NULL, *final_label=NULL, *ch_label=NULL;
  Image    *img=NULL;
  CImage   *cimg=NULL;
  Features *feat=NULL;
  Subgraph *sg=NULL, *sgtrain=NULL, *sgeval=NULL;
  Set      *Obj=NULL,*Bkg=NULL;
  AdjRel   *A=NULL;
  char outfile[100];
  char *file_noext;

  /*--------------------------------------------------------*/

  void *trash = malloc(1);
  struct mallinfo info;
  int MemDinInicial, MemDinFinal;
  free(trash);
  info = mallinfo();
  MemDinInicial = info.uordblks;

  /*--------------------------------------------------------*/

  if (argc!=3){
    fprintf(stderr,"Usage: classify <image.pgm (.ppm)>  <seeds.txt>\n");
    fprintf(stderr,"image.pgm: image to be classified\n");
    fprintf(stderr,"seeds.txt: seed pixels\n");
    exit(-1);
  }

  char *ext = strrchr(argv[1],'.');


  if(!strcmp(ext,".pgm"))
  {
    img   = ReadImage(argv[1]);
    feat  = GaussImageFeats(img, 2);
  }else{
    cimg   = ReadCImage(argv[1]);
    img    = CopyImage(cimg->C[1]);
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

  A = Circular(2);
  ReadSeeds(argv[2],&Obj,&Bkg);

  sg = SubgraphFromSeeds(feat,Obj,Bkg);
  SplitSubgraph(sg, &sgtrain, &sgeval, 0.2);

  t1 = Tic();

  /* OPF-based binary classification of the image */

  OPFLearning(&sgtrain, &sgeval);
  label = OPFClassifyImage(sgtrain, feat);

  //eliminating noise
  ch_label = CloseHoles(label);
  final_label = OpenRec(ch_label,A);

  t2 = Toc();

  fprintf(stdout,"Classification and post-processing time in %f ms\n",CTime(t1,t2));

  cimg = DrawLabeledRegions(img,final_label);
  sprintf(outfile,"%s_result.ppm",strtok(argv[1],"."));
  WriteCImage(cimg,outfile);

  DestroyImage(&label);
  DestroyImage(&img);
  DestroyCImage(&cimg);
  DestroyImage(&final_label);
  DestroyImage(&ch_label);
  DestroySubgraph(&sg);
  DestroySubgraph(&sgtrain);
  DestroySubgraph(&sgeval);
  DestroySet(&Obj);
  DestroySet(&Bkg);
  DestroyFeatures(&feat);
  DestroyAdjRel(&A);

  /* ---------------------------------------------------------- */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);

  return(0);
}



