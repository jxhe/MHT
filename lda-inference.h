#ifndef lda_inference_h
#define lda_inference_h

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "utils.h"
#include "dataset.h"
#include "omp.h"


extern float VAR_CONVERGED;
extern int VAR_MAX_ITER;

double lda_inference(document *doc,vector<double> &var_gamma, vector<vector<double> > &phi,
                     vector<vector<double> > &lamda, vector<vector<double> > &sigma,lda_model *model,
                     int **topic_assign,int **topic_assign_citation,int thread_number,bool flag);

double compute_likelihood(document* doc, lda_model* model, vector<vector<double> > &phi, vector<double> &var_gamma,
                          vector<vector<double> > &sigma, vector<vector<double> > &lamda,int thread_number);

#endif /* lda_inference_h */
