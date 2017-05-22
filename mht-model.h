#ifndef MHT_MODEL_H
#define MHT_MODEL_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dataset.h"
#include "mht-alpha.h"
#include "cokus.h"

#define myrand() (double) (((unsigned long) randomMT()) / 4294967296.)
#define NUM_INIT 1000

void free_lda_model(lda_model*);
void save_lda_model(lda_model*, char*);
lda_model* new_lda_model(int num_terms, int num_docs,int num_topics_word,int num_topics_citatiion);
lda_suffstats new_lda_suffstats(lda_model* model);
void corpus_initialize_ss(lda_suffstats &ss, lda_model* model, dataset* c);
void random_initialize_ss(lda_suffstats &ss, lda_model* model);
void eta5_initialize_ss(lda_suffstats &ss, lda_model* model);
void fix_initialize_ss(lda_suffstats &ss, lda_model* model);
void zero_initialize_ss(lda_suffstats &ss, lda_model* model);
void lda_mle(lda_model* model, const lda_suffstats &ss, int estimate_alpha, int thread_number);
lda_model* load_lda_model(char* model_root);

#endif
