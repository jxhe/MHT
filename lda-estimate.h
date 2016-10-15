#ifndef lda_estimate_h
#define lda_estimate_h

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>

#include "utils.h"
#include "dataset.h"
#include "lda-inference.h"
#include "lda-model.h"
#include "lda-data.h"

int LAG=30;

extern float EM_CONVERGED;
extern int EM_MAX_ITER;
extern int ESTIMATE_ALPHA;
extern double INITIAL_ALPHA;
extern int NTOPICS_WORD;
extern int NTOPICS_CITATION;
extern int NPERIODS;


void run_em(char* start,
            char* directory,
            dataset* corpus);

double doc_e_step(document *doc, vector<double> &gamma,
                  vector<vector<double> > &phi, vector<vector<double> > &lamda, vector<vector<double> > &sigma,
                  lda_model *model,lda_suffstats &ss,int **topic_assign,
                  int **topic_assign_citation,int thread_number,bool flag);

void read_settings(char* filename);

void write_assignment(char *directory, dataset *dataset , int** topic_assign,
                           int **topic_assign_citation,lda_model* model);

void save_gamma(char* filename,
                double** gamma,
                int num_docs,
                int num_topics);

#endif /* lda_estimate_h */
