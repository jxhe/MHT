#include <stdio.h>
#include <stdio.h>
#include "dataset.h"
#include <iostream>
#include "lda-estimate.h"
#include "lda-inference.h"
#include <fstream>
#include "string.h"

float VAR_CONVERGED;
float EM_CONVERGED=0.0001;
int EM_MAX_ITER=50;
int ESTIMATE_ALPHA=0.01;
double INITIAL_ALPHA;
int NTOPICS_WORD;
int NTOPICS_CITATION;
int NPERIODS;

void save_gamma(char* filename, double** gamma, int num_docs, int num_topics)
{
    FILE* fileptr;
    int d, k;
    fileptr = fopen(filename, "w");
    
    for (d = 0; d < num_docs; d++)
    {
        fprintf(fileptr, "%5.10f", gamma[d][0]);
        for (k = 1; k < num_topics; k++)
        {
            fprintf(fileptr, " %5.10f", gamma[d][k]);
        }
        fprintf(fileptr, "\n");
    }
    fclose(fileptr);
}


void run_em(char *start, char *directory, dataset *dataset)
{
    int n,l;
    int **topic_assign;
    int **topic_assign_citation;
    topic_assign = new int*[dataset->num_docs];
    topic_assign_citation=new int*[dataset->num_docs];
    for(int d = 0; d < dataset->num_docs; d++)
    {
        topic_assign[d] = new int[dataset->docs[d]->length];
        topic_assign_citation[d]=new int[dataset->docs[d]->cit_length];
    }
    
    
    
    // double **var_gamma,**phi,**lamda,**sigma;
    // vector< vector<double> > var_gamma, phi, lamda, sigma;
    // var_gamma.assign(dataset->num_docs, vector<double>(NTOPICS_WORD, 0));


    // var_gamma = new double*[dataset->num_docs];
    //    for (d = 0; d < dataset->num_docs; d++) {
    //     var_gamma[d] = new double[NTOPICS_WORD];
    // }
    
    int max_length = max_corpus_length(dataset);
    // phi = new double*[max_length];
    
    // for (n = 0; n < max_length; n++) {
    //     phi[n] = new double[NTOPICS_WORD];
    // }
    // phi.assign(max_length, vector<double>(NTOPICS_WORD, 0));
    
    int max_citlength=max_corpus_citlength(dataset);
    // lamda=  new double*[max_citlength];
    // sigma= new double*[max_citlength];
    
    // for (l=0; l<max_citlength; l++) {
    //     lamda[l]= new double[NTOPICS_WORD];
    //     sigma[l]= new double[NTOPICS_CITATION];
    // }
    // lamda.assign(max_citlength, vector<double>(NTOPICS_WORD, 0));
    // sigma.assign(max_citlength, vector<double>(NTOPICS_CITATION, 0));

    
    lda_model *model ;
    lda_suffstats ss;
    
    // initialize model
    
    char filename[200];
    
    model = new_lda_model(dataset->num_terms, dataset->num_docs,NTOPICS_WORD,NTOPICS_CITATION);

    ss = new_lda_suffstats(model);
    // random_initialize_ss(ss, model);
    if (strcmp(start,"corpus")==0)
    {
        corpus_initialize_ss(ss,model,dataset);
        cout<<"corpus initialize ss"<<endl;
    }
    else if (strcmp(start,"random")==0)
    {
        random_initialize_ss(ss,model);
        cout<<"random initialize ss"<<endl;
    }
    else if (strcmp(start,"eta5")==0)
    {
        eta5_initialize_ss(ss,model);
        cout<<"eta 0.5 initialize ss"<<endl;
    }
    else if (strcmp(start,"fix")==0)
    {
        fix_initialize_ss(ss,model);
        cout<<"fix initialize ss"<<endl;
    }
    else
    {
        cout<<"initialize error"<<endl;
        return;
    }
    int thread_number=omp_get_num_procs();
    lda_mle(model, ss, 0, thread_number);
    model->alpha = INITIAL_ALPHA;
    cout<<1<<endl;
    
    
    int i = 0;
    double likelihood, likelihood_old = 1, converged = 1;
    sprintf(filename, "%s/likelihood.dat", directory);
    FILE* likelihood_file = fopen(filename, "w");

    
    cout<<"thread number"<<thread_number<<endl;
    while (((converged < 0) || (converged > EM_CONVERGED) || (i <= 2)) && (i <= EM_MAX_ITER))
    {
        i++; printf("**** em iteration %d ****\n", i);
        likelihood = 0;
        zero_initialize_ss(ss, model);
        // e-step,update phi and gamma
        
        //The boolean parameter is to determine assigning topics
        #pragma omp parallel num_threads(thread_number) reduction(+:likelihood)
        {
        	lda_suffstats ss_tmp;
        	ss_tmp=new_lda_suffstats(model);
        	double tmp2=0;
        	 // #pragma omp for nowait
        	// #pragma omp parallel for num_threads(thread_number) reduction(+:tmp2)
        	#pragma omp for nowait
	        for (int d=0; d<dataset->num_docs; d++) {
	        	vector< vector<double> > var_gamma, phi, lamda, sigma;
	        	phi.assign(max_length, vector<double>(NTOPICS_WORD, 0));
	        	var_gamma.assign(dataset->num_docs, vector<double>(NTOPICS_WORD, 0));
	        	lamda.assign(max_citlength, vector<double>(NTOPICS_WORD, 0));
    			sigma.assign(max_citlength, vector<double>(NTOPICS_CITATION, 0));
	            if ((d % 1000) == 0) printf("document %d\n",d);
	            likelihood += doc_e_step((dataset->docs[d]),
	                                     var_gamma[d],
	                                     phi,
	                                     lamda,
	                                     sigma,
	                                     model,
	                                     ss_tmp,
	                                     topic_assign,
	                                     topic_assign_citation,
	                                     thread_number,
	                                     false);
	        	// phi[0][0]=0;
	        	// cout<<0<<endl;
        	}

        	#pragma omp critical
        	{
        		ss+=ss_tmp;
        	}
    	}
        // m-step,update alpha and beta
        
        lda_mle(model, ss, ESTIMATE_ALPHA, thread_number);
        
        // check for convergence
        
        converged = (likelihood_old - likelihood) / (likelihood_old);
        
        likelihood_old = likelihood;
        
        // output model and likelihood
        cout<<likelihood<<' '<<converged<<endl;
        
        
        fprintf(likelihood_file, "%10.10f\t%5.5e\n", likelihood, converged);
        fflush(likelihood_file);
        if ((i % LAG) == 0)
        {
            sprintf(filename,"%s/%03d",directory, i);
            save_lda_model(model, filename);
//            sprintf(filename,"%s/%03d.gamma",directory, i);
//            save_gamma(filename, var_gamma, dataset->num_docs, model->num_topics_word);
        }
        if (converged < 0)
        {
            sprintf(filename,"%s/%s",directory, "negative");
            save_lda_model(model, filename);
//            sprintf(filename,"%s/%s.gamma",directory, "negative");
//            save_gamma(filename, var_gamma, dataset->num_docs, model->num_topics_word);
            break;
//            VAR_MAX_ITER = VAR_MAX_ITER * 2;
        }
    }
    
    //assign topics
    int tmp=0;
    zero_initialize_ss(ss, model);
    #pragma omp parallel num_threads(thread_number) reduction(+:tmp)
    {
        lda_suffstats ss_tmp;
        ss_tmp=new_lda_suffstats(model);
        #pragma omp for nowait
        for (int d=0; d<dataset->num_docs; d++) {
            vector< vector<double> > var_gamma, phi, lamda, sigma;
            phi.assign(max_length, vector<double>(NTOPICS_WORD, 0));
            var_gamma.assign(dataset->num_docs, vector<double>(NTOPICS_WORD, 0));
            lamda.assign(max_citlength, vector<double>(NTOPICS_WORD, 0));
            sigma.assign(max_citlength, vector<double>(NTOPICS_CITATION, 0));
            if ((d % 1000) == 0) printf("document %d\n",d);
            tmp += doc_e_step((dataset->docs[d]),
                                     var_gamma[d],
                                     phi,
                                     lamda,
                                     sigma,
                                     model,
                                     ss_tmp,
                                     topic_assign,
                                     topic_assign_citation,
                                     thread_number,
                                     true);
        }

        #pragma omp critical
        {
            ss+=ss_tmp;
        }
    }
    
    
    // output the final model
    
    sprintf(filename,"%s/final",directory);
    save_lda_model(model, filename);
    
    write_assignment(directory, dataset, topic_assign, topic_assign_citation,model);
    free_lda_model(model);
}

double doc_e_step(document *doc, vector<double> &gamma,
                  vector<vector<double> > &phi, vector<vector<double> > &lamda, vector<vector<double> > &sigma,
                  lda_model *model,lda_suffstats &ss,int **topic_assign,
                  int **topic_assign_citation,int thread_number,bool flag)
{
    double likelihood;
    likelihood=lda_inference(doc,gamma,phi,lamda,sigma,model,topic_assign,topic_assign_citation,thread_number,flag);
    int k=0,n=0,i=0,j=0,l=0,p=0,q=0;
    // update sufficient statistics
    double tmp1,tmp2;
    double gamma_sum = 0;
    for (k = 0; k < model->num_topics_word; k++)
    {
        gamma_sum += gamma[k];
        ss.alpha_suffstats += digamma(gamma[k]);
    }
    ss.alpha_suffstats -= model->num_topics_word * digamma(gamma_sum);
    
    for (n = 0; n < doc->length; n++)
    {
        for (k = 0; k < model->num_topics_word; k++)
        {
            tmp1=doc->counts[n]*phi[n][k];
            ss.class_word[k][doc->words[n]] += tmp1;
            ss.class_total[k] +=tmp1;
        }
    }
    
    for (i=0; i<model->num_topics_word; i++) {
        for (j=0; j<model->num_topics_citation; j++) {
            for (l=0; l<doc->cit_length; l++) {
                tmp2=sigma[l][j]*lamda[l][i];
                ss.class_topic[i][j]+=tmp2;
                ss.class_topic_total[i]+=tmp2;
            }
        }
    }
    
    for (p=0; p<doc->cit_length; p++) {
        for (q=0; q<model->num_topics_citation; q++) {
            ss.class_doc[q][doc->citation[p]]+=sigma[p][q];
            ss.class_doc_total[q]+=sigma[p][q];
        }
    }
    
    // ss.num_docs = ss.num_docs + 1;
    
    return(likelihood);
}


void write_assignment(char* directory, dataset *dataset , int** topic_assign,
                           int **topic_assign_citation,lda_model* model)
{
    int n;
    int d;
    char filename[200];
    ofstream fout1;
    ofstream fout2;
    sprintf(filename, "%s/word_assign.dat",directory);
    fout1.open(filename);
    sprintf(filename, "%s/citation_assign.dat",directory);
    fout2.open(filename);
    
    //write_assignment
    for (d = 0; d < dataset->num_docs; d++)
    {
        for (n=0; n<dataset->docs[d]->length; n++) {
            fout1<<topic_assign[d][n]<<' '<<dataset->docs[d]->counts[n]<<' ';
        }
        fout1<<'\n';
    }
    
    for (d = 0; d < dataset->num_docs; d++)
    {
        for (n=0; n<dataset->docs[d]->cit_length; n++) {
            fout2<<topic_assign_citation[d][n]<<' ';
        }
        fout2<<'\n';
    }


    fout1.close();
    fout2.close();
}




int main(int argc, const char * argv[]) {
    //cout<<"sb";
    INITIAL_ALPHA = 0.01;
    sscanf(argv[3],"%d",&NTOPICS_CITATION);
    sscanf(argv[2],"%d",&NTOPICS_WORD);
    sscanf(argv[4],"%f",&VAR_CONVERGED);

    dataset* corpus;
    long t1;
    (void) time(&t1);
    seedMT(t1);
// 
    char abstract[100];
    char citation[100];
    sprintf(abstract,"../Dataset/%s/abstract_ID.txt",argv[1]);
    sprintf(citation,"../Dataset/%s/reference.txt",argv[1]);
    cout<<"reading from "<<argv[1]<<endl;

    corpus=read_data(abstract,citation);
    
    cout<<"Read Complete!"<<endl;

    cout<<"Num_Docs:"<<corpus->num_docs<<'\n';
    cout<<"Num_Terms:"<<corpus->num_terms<<'\n';
    
    cout<<'\n';

    // for (int i=0; i<10; i++) {
    //     cout<<"This is document"<<corpus->docs[i]->ID<<'\n';
    //     cout<<"CitationLength "<<corpus->docs[i]->cit_length<<'\n';
    //     cout<<"CitationList:"<<'\n';
    //     for (int j=0; j<corpus->docs[i]->cit_length; j++) {
    //         cout << corpus->docs[i]->citation[j]<<' ';
    //     }
    //     cout<<'\n';
    //     cout<<"length:"<<corpus->docs[i]->length<<'\n';
    //     cout<<"WordList"<<'\n';
    //     for (int k=0; k<corpus->docs[i]->length; k++) {
    //         cout<<corpus->docs[i]->words[k]<<' '<<"count:"<<corpus->docs[i]->counts[k]<<'\n';
    //     }
    //     cout<<"--------------------------------------------\n";
    // }
    // cout<<'\n';
    
    cout<<argv[1]<<endl;
    cout<<"Word_Topic_Num"<<NTOPICS_WORD<<endl;
    cout<<"Citation_Topic_Num"<<NTOPICS_CITATION<<endl;
    cout<<"VAR_CONVERGED"<<VAR_CONVERGED<<endl;
    char outpath[100];
    sprintf(outpath,"../%s_%s",argv[1],argv[2]);
    mkdir(outpath, S_IRUSR|S_IWUSR|S_IXUSR);
    cout<<"start!"<<endl;
    run_em("random", outpath, corpus);
    cout<<"end"<<endl;
    return 0;
}

