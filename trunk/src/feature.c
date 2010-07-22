#include "feature.h"

Features* CreateFeatures(int ncols, int nrows, int nfeats)
{
    Features *f=(Features *)calloc(1,sizeof(Features));

    f->ncols  = ncols;
    f->nrows  = nrows;
    f->nelems = ncols*nrows;
    f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
    f->nfeats = nfeats;
    f->Imax = 0;

    int i;

    for (i=0; i < f->nelems; i++)
        f->elem[i].feat = AllocFloatArray(nfeats);

    return f;
}

void DestroyFeatures(Features **f)
{
    int i;

    if ((*f)!=NULL)
    {
        for (i=0; i < (*f)->nelems; i++)
            free((*f)->elem[i].feat);
        free((*f)->elem);
        free(*f);
        (*f)=NULL;
    }
}

Features* CopyFeatures(Features* feat)
{
    if(feat == NULL)
    {
        Error("Features* feat NULL ","CopyFeatures");
    }

    Features* result = CreateFeatures(feat->ncols, feat->nrows, feat->nfeats);

    result->Imax = feat->Imax;

    int i;
    int size = feat->nfeats*sizeof(float);

    for(i = 0; i < feat->nelems; i++)
        memcpy(result->elem[i].feat, feat->elem[i].feat,size);

    return result;
}

// Using Linear Convolution with Gaussian filters

Features *LMSImageFeats(Image *img, int nscales)
{
    Features *f=CreateFeatures(img->ncols, img->nrows, nscales);
    AdjRel   *A=NULL;
    int       s,i,p,q;
    Pixel     u,v;
    float    *w,d,K1,K2,sigma2,val;

    f->Imax = MaximumValue(img);

    for (s=1; s <= nscales; s=s+1)
    {
        A  = Circular(s);
        w  = AllocFloatArray(A->n);
        sigma2 = (s/3.0)*(s/3.0);
        K1     =  2.0*sigma2;
        K2     = 1.0/sqrt(2.0*PI*sigma2);
        //compute kernel coefficients

        for (i=0; i < A->n; i++)
        {
            d    = A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i];
            w[i] = K2 * exp(-d/K1); // Gaussian
        }

        // Convolution

        for (p=0; p < f->nelems; p++)
        {
            u.x = p%f->ncols;
            u.y = p/f->ncols;
            val = 0.0;
            for (i=0; i < A->n; i++)
            {
                v.x = u.x + A->dx[i];
                v.y = u.y + A->dy[i];
                if (ValidPixel(img,v.x,v.y))
                {
                    q   = v.x + img->tbrow[v.y];
                    val += (float)img->val[q]*w[i];
                }
            }
            f->elem[p].feat[s-1]=(int)val;
        }
        free(w);

        DestroyAdjRel(&A);
    }
    return(f);
}
