#include <Rcpp.h>
#include "kseq.h"
#include <iostream>
#include <fcntl.h>
#include <stdio.h>
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

using namespace Rcpp;


// [[Rcpp::export]]
List import_fasta_to_vector_each_nt(std::string file) {

  //Initially run through fasta to get consensus sequence and dimensions of matrix
  int n = 0;
  int l = 0;

  gzFile fp;
  kseq_t *seq;

  const char * f = file.c_str();
  fp = gzopen(f, "r");

  seq = kseq_init(fp);

  int seq_length = -1;
  std::vector<std::vector<int> > allele_counts;

  while((l = kseq_read(seq)) >= 0) {

    if (seq_length==-1){
      // initialise
      seq_length = strlen(seq->seq.s);
      allele_counts = std::vector<std::vector<int> >(
          5,
          std::vector<int>(seq_length, 0));
    }

    for(int j=0; j<seq_length; j++){
      if((seq->seq.s[j]=='a') || (seq->seq.s[j]=='A')){
        allele_counts[0][j] += 1;
      } else if((seq->seq.s[j]=='c') || (seq->seq.s[j]=='C')){
        allele_counts[1][j] += 1;
      } else if((seq->seq.s[j]=='g') || (seq->seq.s[j]=='G')){
        allele_counts[2][j] += 1;
      } else if((seq->seq.s[j]=='t') || (seq->seq.s[j]=='T')){
        allele_counts[3][j] += 1;
      } else {
        allele_counts[4][j] += 1;
      }
    }
    n++;
  }
  kseq_destroy(seq);
  gzclose(fp);

  // Now calculate the consensus sequence
  NumericVector consensus(seq_length);
  Rcpp::StringVector seq_names(n);

  int max_allele = -1;

  for(int j=0; j<seq_length; j++){
    for(int i=0; i<5; i++){
      if(allele_counts[i][j]>max_allele){
        max_allele = allele_counts[i][j];
        consensus(j) = i;
      }
    }
    max_allele = -1;
  }

  // create matrix to store snp locations
  std::vector<int> m_i;
  std::vector<int> m_j;
  std::vector<int> m_x;
  m_i.reserve(2*10000);
  m_j.reserve(2*10000);
  m_x.reserve(2*10000);

  // open a new stream to the fasta file TODO: I think there's a cleaner way of doing this.
  fp = gzopen(f, "r");
  seq = kseq_init(fp);

  char temp_char;
  int n_snps = 0;
  n=1;
  while((l = kseq_read(seq)) >= 0) {
    // Record sequence names
    seq_names[n-1] = seq->name.s;

    if ((m_i.capacity() - 2*n_snps)<100){
      // Reserve memory
      m_i.reserve( 2*m_i.capacity() );
      m_j.reserve( 2*m_j.capacity() );
      m_x.reserve( 2*m_x.capacity() );
    }
    for(int j=0; j<seq_length; j++){
      temp_char = seq->seq.s[j];
      if((((temp_char=='A') || (temp_char=='a')) && (consensus[j]!=0)) && (consensus[j]!=4) && (allele_counts[0][j]>1)){
        m_i.push_back(n);
        m_j.push_back(j+1);
        m_x.push_back(1);
        n_snps += 1;
      } else if((((temp_char=='C') || (temp_char=='c')) && (consensus[j]!=1)) && (consensus[j]!=4) && (allele_counts[1][j]>1)){
        m_i.push_back(n);
        m_j.push_back(j+1);
        m_x.push_back(2);
        n_snps += 1;
      } else if((((temp_char=='G') || (temp_char=='g')) && (consensus[j]!=2)) && (consensus[j]!=4) && (allele_counts[2][j]>1)){
        m_i.push_back(n);
        m_j.push_back(j+1);
        m_x.push_back(3);
        n_snps += 1;
      } else if((((temp_char=='T') || (temp_char=='t')) && (consensus[j]!=3)) && (consensus[j]!=4) && (allele_counts[3][j]>1)){
        m_i.push_back(n);
        m_j.push_back(j+1);
        m_x.push_back(4);
        n_snps += 1;
      }
    }
    n += 1;
  }
  kseq_destroy(seq); // STEP 5: destroy seq
  gzclose(fp);

  return List::create(Named("i") = wrap(m_i),
                      Named("j") = wrap(m_j),
                      Named("x") = wrap(m_x),
                      Named("num.seqs") = n-1,
                      Named("consensus") = consensus,
                      Named("seq.length") = seq_length,
                      Named("seq.names") = seq_names);

}


