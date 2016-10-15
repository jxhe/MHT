#ifndef dataset_h
#define dataset_h
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
using namespace std;



typedef struct{
    int *words;
    int *citation;
    int cit_length;
    int *counts;
    int ID;
    
    int length;         //number of different kinds of words
}document;


typedef struct{
    document **docs;
    
    int num_docs;
    int num_terms;
}dataset;


typedef struct
{
    double** class_doc;
    double* class_doc_total;
}omega_suffstats;


typedef struct
{
    double alpha;
    double** log_prob_w;
    int num_topics_word;
    int num_topics_citation;
    int num_terms;
    int num_docs;
    
    ///new variables are added here
    
    double **log_eta;// a matrice K*K with k=num_topic
    double **log_omega;
    
} lda_model;

class lda_suffstats
{
public:
    vector<vector<double> > class_word;
    vector<double> class_total;
    double alpha_suffstats;
    
    vector<vector<double> > class_topic;            //newly added
    vector<double> class_topic_total;       //newly added
    
//    omega_suffstats **class_omega;   //newly added
    vector<vector<double> >class_doc;
    vector<double> class_doc_total;

public:
    lda_suffstats(){}
    ~lda_suffstats(){}
    lda_suffstats operator+(const lda_suffstats &r1) const;
    lda_suffstats& operator+=(const lda_suffstats &r1);
};

dataset *read_data(char *abstract,char *citation);

#endif /* dataset_h */
