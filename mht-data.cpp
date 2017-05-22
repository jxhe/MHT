#include <stdio.h>
#include "mht-data.h"


int max_corpus_length(dataset* c)
{
    int n, max = 0;
    for (n = 0; n < c->num_docs; n++)
        if (c->docs[n]->length > max) max = c->docs[n]->length;
    return(max);
}

int max_corpus_citlength(dataset *c)
{
    int n,max=0;
    for (n=0; n<c->num_docs; n++) {
        if (c->docs[n]->cit_length>max) {
            max=c->docs[n]->cit_length;
        }
    }
    
    return max;
}