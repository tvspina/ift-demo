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
    feat = GaussImageFeats(img, 2);
    DestroyImage(&img);  
  }else{
    CImage   *cimg=NULL;
    cimg   = ReadCImage(argv[1]);
    feat = GaussCImageFeats(cimg, 2);
    DestroyCImage(&cimg);
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



