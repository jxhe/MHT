#include "lda-inference.h"
#include <iostream>
#include <fstream>
using namespace std;
// float VAR_CONVERGED=0.00000001;
int VAR_MAX_ITER=100;


double lda_inference(document *doc,vector<double> &var_gamma, vector<vector<double> > &phi,
                     vector<vector<double> > &lamda, vector<vector<double> > &sigma,lda_model *model,int **topic_assign,
                     int **topic_assign_citation,int thread_number,
                     bool flag)
{
    double converged = 1;
    double phisum = 0, likelihood = 0,lamdasum=0,sigmasum=0;
    double likelihood_old = 0;
    int var_iter;
    double digamma_gam[model->num_topics_word],tmp1,tmp2,phisum_col,lamdasum_col;
    
    // compute posterior dirichlet
    
    ofstream fout;

    // #pragma omp parallel for num_threads(thread_number)
    for (int k = 0; k < model->num_topics_word; k++)
    {
        var_gamma[k] = model->alpha +doc->length/(double)model->num_topics_word+doc->cit_length/(double)model->num_topics_citation;
        digamma_gam[k] = digamma(var_gamma[k]);
        // #pragma omp parallel for num_threads(thread_number)
        for (int n = 0; n < doc->length; n++){
            phi[n][k] = 1.0/model->num_topics_word;
        }
        // #pragma omp parallel for num_threads(thread_number)
        for (int l=0; l<doc->cit_length; l++){
            lamda[l][k]=1.0/model->num_topics_word;
        }
    }
    
    // #pragma omp parallel for num_threads(thread_number)
    for (int k=0; k<model->num_topics_citation; k++) {
        // #pragma omp parallel for num_threads(thread_number)
        for (int l=0; l<doc->cit_length; l++) {
            sigma[l][k]=1.0/model->num_topics_citation;
        }
    }
    var_iter = 0;
    
    while ((converged > VAR_CONVERGED) &&
           ((var_iter < VAR_MAX_ITER) || (VAR_MAX_ITER == -1)))
    {
        var_iter++;
        
        //update phi
        // phisum=0;
        // #pragma omp parallel for num_threads(thread_number) private(phisum)
        for (int n = 0; n < doc->length; n++)
        {
            phisum = 0;
            for (int k = 0; k < model->num_topics_word; k++)
            {
                phi[n][k] =
                digamma_gam[k] +
                model->log_prob_w[k][doc->words[n]];
                
                if (k > 0)
                    phisum = log_sum(phisum, phi[n][k]);
                else
                    phisum = phi[n][k]; // note, phi is in log space
            }
            
            for (int k=0; k<model->num_topics_word; k++) {
                phi[n][k]=exp(phi[n][k]-phisum);
            }

        }
        
       
        
        // update lamda and sigma
        // #pragma omp parallel for num_threads(thread_number) private(lamdasum,sigmasum)
        for (int l=0; l<doc->cit_length; l++) {
            lamdasum=0;
            sigmasum=0;
            // #pragma omp parallel for num_threads(thread_number) private(tmp1)
            for (int i=0; i<model->num_topics_word; i++) {
                tmp1=0;
                for (int j=0; j<model->num_topics_citation; j++) {
                    tmp1+=sigma[l][j]*model->log_eta[i][j];
                }
                lamda[l][i]=digamma_gam[i]+tmp1;
                if (i>0)
                    lamdasum=log_sum(lamdasum,lamda[l][i]);     //log(exp(a)+exp(b))
                else
                    lamdasum=lamda[l][i];
            }
            
            for (int i=0;i<model->num_topics_word; i++)
                lamda[l][i]=exp(lamda[l][i]-lamdasum);
            
    
            // #pragma omp parallel for num_threads(thread_number) private(tmp2)
            for (int i=0;i<model->num_topics_citation;i++)
            {
                tmp2=0;
                for (int j=0;j<model->num_topics_word;j++)
                {
                    tmp2+=lamda[l][j]*model->log_eta[j][i];
                }
                sigma[l][i]=model->log_omega[i][doc->citation[l]]+tmp2;
                if (i>0)
                    sigmasum=log_sum(sigmasum,sigma[l][i]);
                else
                    sigmasum=sigma[l][i];
            }
            
            
            for (int i=0;i<model->num_topics_citation; i++)
                sigma[l][i]=exp(sigma[l][i]-sigmasum);
            
        }
        
        if(flag)
        {
            for (int n = 0; n < doc->length; n++)
            {

                topic_assign[doc->ID][n] =argmax(phi[n], model->num_topics_word);

            }
            
            for (int l=0; l<doc->cit_length; l++) {
                topic_assign_citation[doc->ID][l]=argmax(sigma[l],model->num_topics_citation);
            }
        }
        
        //update var_gamma
        // #pragma omp parallel for num_threads(thread_number) private(phisum_col,lamdasum_col)
        for (int i=0; i<model->num_topics_word; i++) {
            phisum_col=0;
            lamdasum_col=0;
            for (int n=0; n<doc->length; n++)
                phisum_col+=doc->counts[n]*phi[n][i];
            for (int l=0; l<doc->cit_length;l++)
                lamdasum_col+=lamda[l][i];
            var_gamma[i]=model->alpha+phisum_col+lamdasum_col;
            digamma_gam[i] = digamma(var_gamma[i]);
        }
        
        likelihood = compute_likelihood(doc, model, phi, var_gamma,sigma,lamda,thread_number);
        assert(!isnan(likelihood));
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;
        
        
    }
    
    return(likelihood);
    
}





double
compute_likelihood(document* doc, lda_model* model, vector<vector<double> > &phi, vector<double> &var_gamma,
                   vector<vector<double> > &sigma, vector<vector<double> > &lamda,int thread_number)
{
    double likelihood = 0, digsum = 0, var_gamma_sum = 0, dig[model->num_topics_word];
    for (int k = 0; k < model->num_topics_word; k++)
    {
        dig[k] = digamma(var_gamma[k]);
        var_gamma_sum += var_gamma[k];
    }
    digsum = digamma(var_gamma_sum);
    
    likelihood =
    lgamma(model->alpha * model -> num_topics_word)
    - model -> num_topics_word * lgamma(model->alpha)
    - (lgamma(var_gamma_sum));
    
    // #pragma omp parallel for num_threads(thread_number) reduction(+:likelihood)
    for (int k = 0; k < model->num_topics_word; k++)
    {
        likelihood +=
        (model->alpha - 1)*(dig[k] - digsum) + lgamma(var_gamma[k])
        - (var_gamma[k] - 1)*(dig[k] - digsum);
        
        for (int n = 0; n < doc->length; n++)
        {
            if (phi[n][k] > 0)                          //应该是不可能<=0的，记得监测一下
            {
                likelihood += doc->counts[n]*
                (phi[n][k]*((dig[k] - digsum) - log(phi[n][k])
                            + model->log_prob_w[k][doc->words[n]]));
            }
        }
    }
    
    // #pragma omp parallel for num_threads(thread_number) reduction(+:likelihood)
    for (int l=0; l<doc->cit_length; l++) {
        // #pragma omp parallel for num_threads(thread_number) reduction(+:likelihood)
        for (int i=0; i<model->num_topics_word; i++) {
            if (lamda[l][i]>0)
                likelihood+=lamda[l][i]*((dig[i]-digsum)-log(lamda[l][i]));
            
            

            
            for (int j=0; j<model->num_topics_citation; j++) {
                if (model->log_eta[i][j]>0) {
                    likelihood+=sigma[l][j]*lamda[l][i]*log(model->log_eta[i][j]);
                }
            }
        }
    }
    
    // #pragma omp parallel for num_threads(thread_number) reduction(+:likelihood)
    for (int l=0; l<doc->cit_length; l++) {
        for (int i=0; i<model->num_topics_citation; i++) {
            if  (sigma[l][i]>0)
                likelihood+=sigma[l][i]*(model->log_omega[i][doc->citation[l]]
                                          -log(sigma[l][i]));
        }
    }
    return(likelihood);
}