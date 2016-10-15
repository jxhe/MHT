
#include "dataset.h"
#include <stdio.h>
#include <stdlib.h>
#include "constant.h"
#include "strtokenizer.h"
#include <iostream>
#include <fstream>
#include <string>


using namespace std;

int transfer(string str)
{
    string tmp_str;
    tmp_str=str;
    int scale=str.length();
    if (str[0]==' ' || str[scale-1]==' ') {
        if (str[0]==' ') {
            tmp_str=str.substr(1);
        }
        if (str[scale-1]==' ') {
            tmp_str=tmp_str.substr(0,tmp_str.length()-1);
        }
    }
    int tmp=0;
    for (int i=0; i<tmp_str.length(); i++) {
            tmp=tmp*10+(str[i]-'0');
    }
    return tmp;
}

int *split(string str,char s,int &cnt)
{
    int *tmp;
    cnt=0;
    int scale;
    int start=0;
    int end=0;
    tmp=new int[2000];
    scale=str.length();
    if (str[scale-1]==' ') {
        str=str.substr(0,scale-1);
    }
    for (int i=0; i<=str.length(); i++) {
        if (str[i]==s || str[i]=='\0')
        {
            tmp[cnt]=transfer(str.substr(start,i-start));
            start=i+1;
            cnt++;
        }
    }
    
    return tmp;
}

dataset *read_data(char *abstract,char *citation)
{
    ifstream fin1;
    ifstream fin2;
    int i,j;
    
    dataset *corpus=new dataset;
    string line;
    int *storage;
    int cnt;
    fin1.open(abstract);
    fin2.open(citation);
    
    if (!fin1)
        cout<<"abstract open failed"<<endl;
    if (!fin2)
        cout<<"citation open failed"<<endl;
    
    fin1>>corpus->num_docs>>corpus->num_terms;
    fin1.get();
    cout<<corpus->num_docs<<' '<<corpus->num_terms<<endl;
    corpus->docs=new document*[corpus->num_docs];
    
    for (i=0; i<corpus->num_docs; i++) {
        corpus->docs[i]=new document;
    }
    for (i=0; i<corpus->num_docs; i++) {
        getline(fin1,line);
        cnt=0;
        storage=split(line,' ',cnt);
        corpus->docs[i]->ID=i;
        corpus->docs[i]->length=cnt/2;
        corpus->docs[i]->words=new int[cnt/2];
        corpus->docs[i]->counts=new int[cnt/2];
        for (j=0; j<cnt; j+=2) {
            corpus->docs[i]->words[(j/2)]=storage[j];
            corpus->docs[i]->counts[(j/2)]=storage[j+1];
        }
        delete [] storage;
    }
    for (i=0; i<corpus->num_docs; i++) {
        getline(fin2,line);
        cnt=0;
        storage=split(line, ' ', cnt);
        corpus->docs[i]->cit_length=cnt-1;
        corpus->docs[i]->citation=new int[cnt-1];
        for (j=0; j<cnt-1; j++) {
            corpus->docs[i]->citation[j]=storage[j+1];
        }
        delete [] storage;
    }
    
    return corpus;
}

lda_suffstats& lda_suffstats::operator+=(const lda_suffstats &r1)
{
    int word_num_topic=(int)class_word.size();
    int num_terms=(int)class_word[0].size();
    int citation_num_topic=(int)class_doc.size();
    int num_docs=(int)class_doc[0].size();

    alpha_suffstats+=r1.alpha_suffstats;
    for (int i = 0; i < word_num_topic; ++i)
    {
        class_total[i]+=r1.class_total[i];
        class_topic_total[i]+=r1.class_topic_total[i];
        for (int j = 0; j < num_terms; ++j)
        {
            class_word[i][j]+=r1.class_word[i][j];
        }

        for (int k = 0; k < citation_num_topic; ++k)
        {
            class_topic[i][k]+=r1.class_topic[i][k];
        }
    }

    for (int i = 0; i < citation_num_topic; ++i)
    {
        class_doc_total[i]+=r1.class_doc_total[i];
        for (int j = 0; j < num_docs; ++j)
        {
            class_doc[i][j]+=r1.class_doc[i][j];
        }
    }

    return *this;
}

lda_suffstats lda_suffstats::operator+(const lda_suffstats &r1) const
{
    int word_num_topic=(int)class_word.size();
    int num_terms=(int)class_word[0].size();
    int citation_num_topic=(int)class_doc.size();
    int num_docs=(int)class_doc[0].size();

    lda_suffstats tmp;
    tmp.alpha_suffstats=r1.alpha_suffstats;
    tmp.class_total.assign(word_num_topic,0);
    tmp.class_word.assign(word_num_topic,vector<double>(num_terms,0));
    tmp.class_doc_total.assign(citation_num_topic,0);
    tmp.class_doc.assign(citation_num_topic,vector<double>(num_docs,0));
    tmp.class_topic_total.assign(word_num_topic,0);
    tmp.class_topic.assign(word_num_topic,vector<double>(citation_num_topic,0));
    for (int i = 0; i < word_num_topic; ++i)
    {
        tmp.class_total[i]=class_total[i]+r1.class_total[i];
        tmp.class_topic_total[i]=class_topic_total[i]+r1.class_topic_total[i];
        for (int j = 0; j < num_terms; ++j)
        {
            tmp.class_word[i][j]=class_word[i][j]+r1.class_word[i][j];
        }

        for (int k = 0; k < citation_num_topic; ++k)
        {
            tmp.class_topic[i][k]=class_topic[i][k]+r1.class_topic[i][k];
        }
    }

    for (int i = 0; i < citation_num_topic; ++i)
    {
        tmp.class_doc_total[i]=class_doc_total[i]+r1.class_doc_total[i];
        for (int j = 0; j < num_docs; ++j)
        {
            tmp.class_doc[i][j]=class_doc[i][j]+r1.class_doc[i][j];
        }
    }

    return tmp;
}