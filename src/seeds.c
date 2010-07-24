#include "seeds.h"

void ReadSeeds(char *filename, Set **Obj, Set **Bkg)
{
  FILE *fp=fopen(filename,"r");
  int i,x,y,l,mk,nseeds, ncols, nrows;
  
  if(fscanf(fp,"%d %d %d",&nseeds, &ncols, &nrows)!=0);

  for (i=0; i < nseeds; i++){
    if(fscanf(fp,"%d %d %d %d",&x,&y,&mk,&l)!=0);
    if (l==0)
      InsertSet(Bkg, x + ncols*y);
    else
      InsertSet(Obj, x + ncols*y);
  }
  fclose(fp);
}
