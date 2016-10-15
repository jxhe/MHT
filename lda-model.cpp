#include <stdio.h>
#include "lda-model.h"
#include <iostream>
/*
 * compute MLE lda model from sufficient statistics
 *
 */


void lda_mle(lda_model* model, const lda_suffstats &ss, int estimate_alpha, int thread_number)
{
    #pragma omp parallel for num_threads(thread_number)
    for (int k = 0; k < model->num_topics_word; k++)
    {
        for (int w = 0; w < model->num_terms; w++)
        {
            if (ss.class_word[k][w] > 0)
            {
                model->log_prob_w[k][w] =
                log(ss.class_word[k][w]) -
                log(ss.class_total[k]);
            }
            else
                model->log_prob_w[k][w] = -100;
        }
    }
    
    //newly added: update Zeta(topic-topic)
    #pragma omp parallel for num_threads(thread_number)
    for (int t1 = 0; t1 < model->num_topics_word; ++t1)
    {
        for (int t2 = 0; t2 < model->num_topics_citation; ++t2)
        {
            if (ss.class_topic[t1][t2] > 0)
            {
                model->log_eta[t1][t2] =
                log(ss.class_topic[t1][t2]) -
                log(ss.class_topic_total[t1]);
            }
            else
            {
                model->log_eta[t1][t2] = -100;
            }
        }
    }
    //
    
    //newly added: update Omega(topic-doc)
    for (int z=0; z<model->num_topics_citation; z++) {
        #pragma omp parallel for num_threads(thread_number)
        for (int d=0; d<model->num_docs; d++) {
            if (ss.class_doc[z][d]>0) {
                model->log_omega[z][d]=log(ss.class_doc[z][d])-log(ss.class_doc_total[z]);
            }
            else
                model->log_omega[z][d]=-100;
        }
    }
    //
    
    
    if (estimate_alpha == 1)
    {
        model->alpha = opt_alpha(ss.alpha_suffstats,
                                 model->num_docs,
                                 model->num_topics_word);
        
        printf("new alpha = %5.5f\n", model->alpha);
    }
}
///*
// * allocate sufficient statistics
// *
// */
//

lda_suffstats new_lda_suffstats(lda_model* model)
{
    int num_topics_words = model->num_topics_word;
    int num_topics_citation=model->num_topics_citation;
    int num_terms = model->num_terms;
    int num_docs=model->num_docs;
    int i,j;
    
    lda_suffstats ss;
    
    // ss.class_total = new double[num_topics_words];
    // ss.class_word = new double*[num_topics_words];
    ss.class_total.assign(num_topics_words,0);
    ss.class_word.assign(num_topics_words,vector<double>(num_terms,0));
    // for (i = 0; i < num_topics_words; i++)
    // {
    //     ss.class_total[i] = 0;
    //     ss.class_word[i] = new double[num_terms];
    //     for (j = 0; j < num_terms; j++)
    //     {
    //         ss.class_word[i][j] = 0;
    //     }
    // }
    
    // ss.class_doc_total=new double[num_topics_citation];
    // ss.class_doc=new double*[num_topics_citation];
    ss.class_doc_total.assign(num_topics_citation,0);
    ss.class_doc.assign(num_topics_citation,vector<double>(num_docs,0));
    // for (i=0; i<num_topics_citation; i++) {
    //     ss.class_doc_total[i]=0;
    //     ss.class_doc[i]=new double[num_docs];
    //     for (j=0; j<num_docs; j++) {
    //         ss.class_doc[i][j]=0;
    //     }
    // }
    
    // ss.class_topic_total = new double[num_topics_words];
    // ss.class_topic = new double*[num_topics_words];
    // for (i = 0; i < num_topics_words; i++)
    // {
    //     ss.class_topic_total[i] = 0;
    //     ss.class_topic[i] = new double[num_topics_citation];
    //     for (j = 0; j < num_topics_citation; j++)
    //     {
    //         ss.class_topic[i][j] = 0;
    //     }
    // }
    ss.class_topic_total.assign(num_topics_words,0);
    ss.class_topic.assign(num_topics_words,vector<double>(num_topics_citation,0));
    
    ss.alpha_suffstats = 0;
    return(ss);
}
//
///*
// * various intializations for the sufficient statistics
// *
// */
//
void zero_initialize_ss(lda_suffstats &ss, lda_model* model)
{
    int k, w;
    for (k = 0; k < model->num_topics_word; k++)
    {
        ss.class_total[k] = 0;
        for (w = 0; w < model->num_terms; w++)
        {
            ss.class_word[k][w] = 0;
        }
    }
    
    int t1, t2;
    for (t1 = 0; t1 < model->num_topics_word; ++t1)
    {
        ss.class_topic_total[t1] = 0;
        for (t2 = 0; t2 < model->num_topics_citation; ++t2)
        {
            ss.class_topic[t1][t2] = 0;
        }
    }
    
    int z,d;
    for (z=0; z<model->num_topics_citation; z++) {
        ss.class_doc_total[z]=0;
        for (d=0; d<model->num_docs; d++) {
            ss.class_doc[z][d]=0;
        }
    }
    
    ss.alpha_suffstats = 0;
}


void eta5_initialize_ss(lda_suffstats &ss, lda_model* model)
{
    int num_topics_word = model->num_topics_word;
    int num_topics_citation=model->num_topics_citation;
    int num_terms = model->num_terms;
    int num_docs=model->num_docs;
    int k, n;
    for (k = 0; k < num_topics_word; k++)
    {
        for (n = 0; n < num_terms; n++)
        {
            ss.class_word[k][n] += 1.0/num_terms + myrand();
            ss.class_total[k] += ss.class_word[k][n];
        }
    }
    
    for (k=0; k<num_topics_citation; k++) {
        for (n=0; n<num_docs; n++) {
            ss.class_doc[k][n]+=1.0/num_docs+myrand();
            ss.class_doc_total[k]+=ss.class_doc[k][n];
        }
    }
    
    int t1, t2,select;
    // for (t1 = 0; t1 < num_topics_word; ++t1)
    // {
    //     for (t2 = 0; t2 < num_topics_citation; ++t2)
    //     {
    //         ss->class_topic[t1][t2] += 1.0/num_topics_citation + myrand();
    //         ss->class_topic_total[t1] += ss->class_topic[t1][t2];
    //     }
    // }

    for (t1 = 0; t1 < num_topics_word; ++t1)
    {
        select=floor(num_topics_word*myrand());
        for (t2 = 0; t2 < num_topics_citation; ++t2)
        {
            if (t2==select)
            {
                ss.class_topic[t1][t2]=num_topics_citation;
            }
            else
                ss.class_topic[t1][t2]=1;

            ss.class_topic_total[t1] += ss.class_topic[t1][t2];
        }
        
        
    }
    
}


void fix_initialize_ss(lda_suffstats &ss, lda_model* model)
{
    int num_topics_word = model->num_topics_word;
    int num_topics_citation=model->num_topics_citation;
    int num_terms = model->num_terms;
    int num_docs=model->num_docs;
    int k, n;
    for (k = 0; k < num_topics_word; k++)
    {
        for (n = 0; n < num_terms; n++)
        {
            ss.class_word[k][n] += 1.0/num_terms + myrand();
            ss.class_total[k] += ss.class_word[k][n];
        }
    }
    
    for (k=0; k<num_topics_citation; k++) {
        for (n=0; n<num_docs; n++) {
            ss.class_doc[k][n]+=1.0/num_docs+myrand();
            ss.class_doc_total[k]+=ss.class_doc[k][n];
        }
    }
    
    int t1, t2,select;
    // for (t1 = 0; t1 < num_topics_word; ++t1)
    // {
    //     for (t2 = 0; t2 < num_topics_citation; ++t2)
    //     {
    //         ss->class_topic[t1][t2] += 1.0/num_topics_citation + myrand();
    //         ss->class_topic_total[t1] += ss->class_topic[t1][t2];
    //     }
    // }

    for (t1 = 0; t1 < num_topics_word; ++t1)
    {
        for (t2 = 0; t2 < num_topics_citation; ++t2)
        {
            if (t2==t1)
            {
                ss.class_topic[t1][t2]=100+myrand();
            }
            else
                ss.class_topic[t1][t2]=1+myrand();
            
            ss.class_topic_total[t1] += ss.class_topic[t1][t2];
        }
        
        
    }
    
}



void random_initialize_ss(lda_suffstats &ss, lda_model* model)
{
    int num_topics_word = model->num_topics_word;
    int num_topics_citation=model->num_topics_citation;
    int num_terms = model->num_terms;
    int num_docs=model->num_docs;
    int k, n;
    for (k = 0; k < num_topics_word; k++)
    {
        for (n = 0; n < num_terms; n++)
        {
            ss.class_word[k][n] += 1.0/num_terms + myrand();
            ss.class_total[k] += ss.class_word[k][n];
        }
    }
    
    for (k=0; k<num_topics_citation; k++) {
        for (n=0; n<num_docs; n++) {
            ss.class_doc[k][n]+=1.0/num_docs+myrand();
            ss.class_doc_total[k]+=ss.class_doc[k][n];
        }
    }
    
    int t1, t2;
    for (t1 = 0; t1 < num_topics_word; ++t1)
    {
        for (t2 = 0; t2 < num_topics_citation; ++t2)
        {
            ss.class_topic[t1][t2] += 1.0/num_topics_citation + myrand();
            ss.class_topic_total[t1] += ss.class_topic[t1][t2];
        }
    }

    
}


void corpus_initialize_ss(lda_suffstats &ss, lda_model* model, dataset* c)
{
   int num_topics_word = model->num_topics_word;
   int num_topics_citation=model->num_topics_citation;
   int i, k, d, n;
   document* doc;
   
   for (k = 0; k < num_topics_word; k++)
   {
       for (i = 0; i < NUM_INIT; i++)
       {
           d = floor(myrand() * c->num_docs);
           // printf("initialized beta with document %d\n", d);
           doc = c->docs[d];
           for (n = 0; n < doc->length; n++)
           {
               ss.class_word[k][doc->words[n]] += doc->counts[n];
           }
       }
       for (n = 0; n < model->num_terms; n++)
       {
           ss.class_word[k][n] += 1.0;
           ss.class_total[k] = ss.class_total[k] + ss.class_word[k][n];
       }
       
       
   }
   

   for(k = 0; k < num_topics_citation; k++)
   {
       for (i = 0; i < NUM_INIT; i++)
       {
            while (c->docs[d]->citation==0)
                d = floor(myrand() * c->num_docs);           
        
            // printf("intialized omega with document %d\n", d);
            doc = c->docs[d];
            for (n = 0; n < doc->cit_length; n++)
            {
                ss.class_doc[k][doc->citation[n]] += 1;
            }
           
       }
       for (n = 0; n < model->num_docs; n++)
       {
           ss.class_doc[k][n] += 1.0;
           ss.class_doc_total[k]+=ss.class_doc[k][n];
       }
       
   }
   
   
   //random initialize
   int t1, t2;
   for (t1 = 0; t1 < num_topics_word; ++t1)
   {
       for (t2 = 0; t2 < num_topics_citation; ++t2)
       {
           ss.class_topic[t1][t2] += 1.0/num_topics_citation + myrand();
           ss.class_topic_total[t1] += ss.class_topic[t1][t2];
       }
   }
}
/*
 * allocate new lda model
 *
 */

lda_model* new_lda_model(int num_terms,int num_docs,int num_topics_word,int num_topics_citation)
{
    int i,j,k;
    lda_model* model;
   // cout<<num_periods;
    
    model = new lda_model;
    model->num_topics_word = num_topics_word;
    model->num_topics_citation=num_topics_citation;
    model->num_terms = num_terms;
    model->num_docs=num_docs;
    model->alpha = 1.0;
    model->log_prob_w = new double*[num_topics_word];
    model->log_eta= new double*[num_topics_word];
    model->log_omega= new double*[num_topics_citation];
    
    for (i = 0; i < num_topics_word; i++)
    {
        model->log_prob_w[i] = new double[num_terms];
        model->log_eta[i]= new double[num_topics_citation];
        for (j = 0; j < num_terms; j++)
            model->log_prob_w[i][j] = 0;
        for (k=0; k<num_topics_citation; k++) {
            model->log_eta[i][k]=0;
        }
    }
    
    for (i=0; i<num_topics_citation; i++) {
        model->log_omega[i]= new double[num_docs];
        for (j=0; j<num_docs; j++) {
            model->log_omega[i][j]=0;
        }
    }
    
    
    return(model);
}

//
///*
// * deallocate new lda model
// *
// */
//
void free_lda_model(lda_model* model)
{
    int i;
    
    for (i = 0; i < model->num_topics_word; i++)
    {
        delete [] model->log_prob_w[i];
        delete [] model->log_eta[i];
    }
    
    delete model->log_prob_w;
    delete model->log_eta;
    
    for (i=0; i<model->num_topics_citation; i++) {
        delete [] model->log_omega[i];
    }
    
    delete [] model->log_omega;
    
    
}



///*
// * save an lda model
// *
// */

void save_lda_model(lda_model* model, char* model_root)
{
    char filename[100];
    FILE* fileptr;
    int i, j;
    
    sprintf(filename, "%s.beta", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics_word; i++)
    {
        fprintf(fileptr, "%5.10f", exp(model->log_prob_w[i][0]));
        for (j = 1; j < model->num_terms; j++)
        {
            fprintf(fileptr, " %5.10f", exp(model->log_prob_w[i][j]));
        }
        if (i!=model->num_topics_word-1) {
            fprintf(fileptr, "\n");
        }
        
    }
    fclose(fileptr);
    
    sprintf(filename, "%s.eta", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics_word; i++)
    {
        fprintf(fileptr, "%5.10f", exp(model->log_eta[i][0]));
        for (j = 1; j < model->num_topics_citation; j++)
        {
            fprintf(fileptr, " %5.10f", exp(model->log_eta[i][j]));
        }
        if (i!=model->num_topics_word-1) {
            fprintf(fileptr, "\n");
        }
    }
    fclose(fileptr);
    
    

    sprintf(filename, "%s.omega", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics_citation; i++)
    {
        fprintf(fileptr, "%5.10f", exp(model->log_omega[i][0]));
        for (j = 1; j < model->num_docs; j++)
        {
            fprintf(fileptr, " %5.10f", exp(model->log_omega[i][j]));
        }
        if (i!=model->num_topics_citation-1) {
            fprintf(fileptr, "\n");
        }
        
    }
    fclose(fileptr);
    
    
    
    sprintf(filename, "%s.other", model_root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "num_topics_word %d\n", model->num_topics_word);
    fprintf(fileptr, "num_topics_citation %d\n", model->num_topics_citation);
    fprintf(fileptr, "num_terms %d\n", model->num_terms);
    fprintf(fileptr, "alpha %5.10f\n", model->alpha);
    fclose(fileptr);
}

