#include "ift.h"

int main(int argc, char **argv) 
{
  timer    *t1=NULL,*t2=NULL;
  Image    *label=NULL, *final_label=NULL, *ch_label=NULL;
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

  A = Circular(2);
  ReadSeeds(argv[2],&Obj,&Bkg);

  sg = SubgraphFromSeeds(feat,Obj,Bkg);
  SplitSubgraph(sg, &sgtrain, &sgeval, 0.2);
  
  t1 = Tic();

  OPFLearning(&sgtrain, &sgeval);
  label = OPFClassifyImage(sgtrain, feat);
  //eliminating noise
  ch_label = CloseHoles(label);
  final_label = OpenRec(ch_label,A);
  
  t2 = Toc();    

  fprintf(stdout,"Classification and post-processing time in %f ms\n",CTime(t1,t2));
  sprintf(outfile,"%s_label.pgm",strtok(argv[1],"."));
  WriteImage(final_label,outfile);    

  DestroyImage(&label);  
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



