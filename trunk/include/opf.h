#ifndef _OPF_H_
#define _OPF_H_

#include "common.h"
#include "subgraph.h"
#include "feature.h"
#include "realheap.h"

#define MAXARCW 100000.0
#define PROTOTYPE 1

void SetSubgraphFeatures(Subgraph *sg, Features *f);

Subgraph* SubgraphFromSeeds(Features* f, Set* Si, Set* Se);

// Split subgraph into two parts such that the size of the first part
// is given by a percentual of samples.
void SplitSubgraph(Subgraph *sg, Subgraph **sg1, Subgraph **sg2, float perc1);

// Compute Euclidean distance between feature vectors
float EuclDist(float *f1, float *f2, int n);

// Discretizes original distance
float EuclDistLog(float *f1, float *f2, int n);

// Compute accuracy
float Accuracy(Subgraph *sg);

//Replace errors from evaluating set by non prototypes from training set
void SwapErrorsbyNonPrototypes(Subgraph **sgtrain, Subgraph **sgeval);

// Find prototypes by the MST approach
void MSTPrototypes(Subgraph *sg);

void OPFTraining(Subgraph *sg);

//Classification function: it simply classifies samples from sg -----
void OPFClassify(Subgraph *sgtrain, Subgraph *sg);

//Learning function: it executes the learning procedure for CompGraph replacing the
//missclassified samples in the evaluation set by non prototypes from
//training set -----
void OPFLearning(Subgraph **sgtrain, Subgraph **sgeval);

//Classification function: returns a binary image where object pixels have value 1 and background pixels value 0
Image* OPFClassifyImage(Subgraph *sgtrain, Features* feat);

#endif
