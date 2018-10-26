#include <Rcpp.h>
#include "kseq.h"
#include <iostream>
#include <fcntl.h>
#include <stdio.h>

using namespace Rcpp;


// [[Rcpp::export]]
List import_fasta_to_vector_each_nt(std::string file) {

  //Initially run through fasta to get consensus sequence and dimensions of matrix
  int n = 0;
  int l = 0;
  kseq seq;
  const char * f = file.c_str();

  int fp = open(f, O_RDONLY);
  FunctorRead r;
  kstream<int, FunctorRead> ks(fp, r);

  l = ks.read(seq);
  int seq_length = seq.seq.length();
  int allele_counts[5][seq_length];
  memset(allele_counts, 0, 5*seq_length*sizeof(int));


  while((l = ks.read(seq)) >= 0) {
    // std::cout << seq.name << std::endl;
    // std::cout << seq.seq << std::endl;
    for(int j=0; j<seq_length; j++){
      if((seq.seq[j]=='a') || (seq.seq[j]=='A')){
        allele_counts[0][j] += 1;
      } else if((seq.seq[j]=='c') || (seq.seq[j]=='C')){
        allele_counts[1][j] += 1;
      } else if((seq.seq[j]=='g') || (seq.seq[j]=='G')){
        allele_counts[2][j] += 1;
      } else if((seq.seq[j]=='t') || (seq.seq[j]=='T')){
        allele_counts[3][j] += 1;
      } else {
        allele_counts[4][j] += 1;
      }
    }
    n++;
  }
  close(fp);

  // Now calculate the consensus sequence
  NumericVector consensus(seq_length);
  Rcpp::StringVector seq_names(n+1);

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
  int fp2 = open(f, O_RDONLY);
  kstream<int, FunctorRead> ks2(fp2, r);
  char temp_char;
  int n_snps = 0;
  n=1;
  while((l = ks2.read(seq)) >= 0) {
    // Record sequence names
    seq_names[n-1] = seq.name;

    if ((m_i.capacity() - 2*n_snps)<100){
      // Reserve memory
      m_i.reserve( 2*m_i.capacity() );
      m_j.reserve( 2*m_j.capacity() );
      m_x.reserve( 2*m_x.capacity() );
    }
    for(int j=0; j<seq_length; j++){
      temp_char = seq.seq[j];
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
  close(fp2);

  return List::create(Named("i") = wrap(m_i),
                      Named("j") = wrap(m_j),
                      Named("x") = wrap(m_x),
                      Named("num.seqs") = n-1,
                      Named("consensus") = consensus,
                      Named("seq.length") = seq_length,
                      Named("seq.names") = seq_names);

}


