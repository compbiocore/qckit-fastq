#include <Rcpp.h>
#include <iostream>
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdint>
#include "gzstream.h"
#include "zlib.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' calculate Over Rep seqs
//'
//' @param infile  A string giving the path for the fastqfile
//' @param out_prefix A string giving the prefix to be used for outputs
//' @param buffer_size An int for the number of lines to keep in memory
//' @examples
//' infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
//' process_fastq(infile,"test",10000)
//' @return process fastq and generate sequence and quality score tables
//' @export
// [[Rcpp::export]]
void process_fastq (std::string infile, std::string out_prefix, int buffer_size)
  {
  

  std::map<std::string, int> over_rep_map;
  std::map<std::string,int>::iterator it;
  std::map<int,std::vector<int> > qual_by_column;
  std::map<int,std::vector<int> >::iterator qual_by_col_it;

  //std::string seq_out = out_prefix + ".seq.csv";
  //std::string qual_char_out = out_prefix + ".qual.char.csv";
  //std::string qual_num_out = out_prefix + ".qual.num.csv";
  //std::ofstream seq_file, qual_char_file, qual_num_file;

  std::string over_rep_out = out_prefix + ".over_rep.csv";
  std::ofstream over_rep_file;

  /*seq_file.open(seq_out.c_str());
  qual_char_file.open(qual_char_out.c_str());
  qual_num_file.open(qual_num_out.c_str());*/
  over_rep_file.open(over_rep_out.c_str());

  gz::igzstream in(infile.c_str());
  std::string line;
  int count = 1, line_count =1;
  //std::vector<int,std::vector<int> > base_counts;
  std::vector<double> gc_percent_all;
  while (std::getline(in, line)) {

    if (count == 2)
      {
        it = over_rep_map.find(line);
        if (it != over_rep_map.end())
        {
          // if found increment by 1
          over_rep_map.at(line) += 1;
        }
        else
        {
          // if not found add new key and initialize to 1
          over_rep_map.insert(std::pair<std::string, int>(line,1));
        }
        // iterating over each character in the string
        std::string base_cmp;
        //std::vector<int> counts_per_read;
        int count_A=0, count_G=0, count_T=0, count_C=0, count_N=0;
        for (std::string::iterator it = line.begin(); it != line.end(); ++it)
          {
            base_cmp.clear();
            base_cmp.push_back(*it);
            if (base_cmp.compare("A") ==0) { count_A +=1;}
            else if (base_cmp.compare("T")==0) { count_T +=1;}
            else if (base_cmp.compare("G")==0) { count_G +=1;}
            else if (base_cmp.compare("C")==0) { count_C +=1;}
            else if (base_cmp.compare("N")==0) { count_N +=1;}
          }
        double gc_percent = static_cast<double>(count_C+count_G)/static_cast<double>(count_A+count_T+count_G+count_C+count_N);
        gc_percent_all.push_back(gc_percent);
        /*counts_per_read.push_back(count_A);
        counts_per_read.push_back(count_T);
        counts_per_read.push_back(count_G);
        counts_per_read.push_back(count_C);
        counts_per_read.push_back(count_N);
        base_counts.insert(line_count,counts_per_read);*/

      }

    if (count == 4)
      {
       // iterate over each value for quality
       int pos_counter = 1;
        for (std::string::iterator it = line.begin(); it != line.end(); ++it)
          {
           qual_by_col_it = qual_by_column.find(pos_counter);
            if(qual_by_col_it != qual_by_column.end())
            {
              qual_by_column.at(pos_counter).push_back(static_cast<int>(*it));
            }
            else
            {
              std::vector<int> tmp;
              tmp.push_back(static_cast<int>(*it));
              qual_by_column.insert(std::pair<int, std::vector<int> >(pos_counter,tmp));
            }
          }
      count = 1;
      }
    else
      {
      count++;
      }
  }
  //Cleanup
  in.close();

  /*TODO
    write over_rep sequences to file
    return gc_percent vector
    return per column mean, median and quantile
   */
  over_rep_file.close();

  return ;
}


//' calculate summary of quality scores over position
//'
//' Description
//' @param inmat A matrix of score vectors per position

std::vector<std::vector<int> > qual_score_per_position (const std::map<int,std::vector<uint8_t> > &inmat)
{
  std::vector<std::vector<int> > qual_score_mat_results;
  std::vector<int> q_01,q_25, q_50,q_75, q_99;

  std::map<int,std::vector<uint8_t> >::const_iterator mat_it = inmat.begin();

  for(mat_it = inmat.begin(); mat_it !=inmat.end(); mat_it++)
  {
    std::vector<uint8_t> quantile = mat_it->second;

    int Q1 = static_cast<int> (quantile.size()*0.01);
    int Q25 = (quantile.size()+1) / 4;
    int Q50= (quantile.size()+1) / 2;
    int Q75 = Q25 + Q50;
    int Q99 = static_cast<int> (quantile.size()*0.99) ;

    std::nth_element(quantile.begin(),quantile.begin() + Q1, quantile.end());
    q_01.push_back(static_cast<int>(quantile[Q1]));
    quantile.clear();

    std::nth_element(quantile.begin(),quantile.begin() + Q25, quantile.end());
    q_25.push_back(static_cast<int>(quantile[Q25]));
    quantile.clear();

    quantile = mat_it->second;
    std::nth_element(quantile.begin(), quantile.begin() + Q50, quantile.end());
    q_50.push_back(static_cast<int>(quantile[Q50]));
    quantile.clear();

    quantile = mat_it->second;
    std::nth_element(quantile.begin(), quantile.begin() + Q75, quantile.end());
    q_75.push_back(static_cast<int>(quantile[Q75]));

    std::nth_element(quantile.begin(),quantile.begin() + Q99, quantile.end());
    q_99.push_back(static_cast<int>(quantile[Q99]));
    quantile.clear();
  }
  qual_score_mat_results.push_back(q_01);
  qual_score_mat_results.push_back(q_25);
  qual_score_mat_results.push_back(q_50);
  qual_score_mat_results.push_back(q_75);
  qual_score_mat_results.push_back(q_99);
  return qual_score_mat_results ;
}

//' calculate mean quality per read
//'
//' Calculate the mean quality score per read of the FASTQ gzipped file
//' @param infile A string giving the path for the fastqfile
//' @examples
//' infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
//' qual_score_per_read(infile)
//' @return mean quality per read
//' @export
// [[Rcpp::export]]
Rcpp::List qual_score_per_read (std::string infile)
{


  std::vector<double> quality_score_per_read;
  std::vector<uint8_t> qual_by_column;
  std::vector<uint8_t>::iterator qual_by_col_it;

  std::map<int,std::vector<uint8_t> > qual_score_matrix;

  gz::igzstream in(infile.c_str());
  std::string line;
  int count = 1;
  double quality_score_mean = 0;
  while (std::getline(in, line))
    {

    if (count == 4)
    {
      // iterate over each value for quality
      qual_by_column.clear();
      int pos_counter = 1;
      for (std::string::iterator it = line.begin(); it != line.end(); ++it)
      {
          qual_by_column.push_back(static_cast<int>(*it));
          if(pos_counter <= qual_score_matrix.size())
          {
          qual_score_matrix[pos_counter].push_back(static_cast<int>(*it));
          }
          else
          {
            std::vector<uint8_t> tmp_qual;
            tmp_qual.push_back(static_cast<int>(*it));
            std::pair<int,std::vector<uint8_t>> tmp = std::pair<int,std::vector<uint8_t>>(pos_counter,tmp_qual);
            qual_score_matrix.insert(tmp);
          }
          pos_counter++;
      }
      quality_score_mean = static_cast<double>(std::accumulate(qual_by_column.begin(), qual_by_column.end(),
                                           0.0))/static_cast<double>(qual_by_column.size());
      quality_score_per_read.push_back(quality_score_mean);
      count = 1;
    }
    else
    {
      count++;
    }
  }
  std::vector<std::vector<int> > qual_score_summary_by_position;
  qual_score_summary_by_position = qual_score_per_position(qual_score_matrix);

  std::vector<double> mu_per_position;
  std::map<int,std::vector<uint8_t> >::iterator mat_it = qual_score_matrix.begin();

  for(mat_it = qual_score_matrix.begin(); mat_it != qual_score_matrix.end(); mat_it++)
  {
    std::vector<uint8_t> quantile = mat_it->second;

    mu_per_position.push_back(static_cast<double>(std::accumulate(quantile.begin(),
                                                                  quantile.end(),0.0))/static_cast<double>(quantile.size()));
  }
  //Cleanup
  in.close();
  return Rcpp::List::create(Rcpp::Named("mu_per_read") = quality_score_per_read,
                            Rcpp::Named("mu_per_position") = mu_per_position,
                            Rcpp::Named("q01_per_position") = qual_score_summary_by_position[0],
                            Rcpp::Named("q25_per_position") = qual_score_summary_by_position[1],
                            Rcpp::Named("q50_per_position") = qual_score_summary_by_position[2],
                            Rcpp::Named("q75_per_position") = qual_score_summary_by_position[3],
                            Rcpp::Named("q99_per_position") = qual_score_summary_by_position[4]);
  }


//' calculate GC percent per read
//'
//' Calculate GC nucleotide sequence content per read of the FASTQ gzipped file
//' @param infile A string giving the path for the fastqfile
//' @examples
//' infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
//' gc_per_read(infile)
//' @return GC content perncentage per read
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gc_per_read (std::string infile)
{

  std::map<int,std::vector<int> > qual_by_column;
  std::map<int,std::vector<int> >::iterator qual_by_col_it;


  gz::igzstream in(infile.c_str());
  std::string line;
  int count = 1, line_count =1;
  //std::vector<int,std::vector<int> > base_counts;
  std::vector<double> gc_percent_per_read;
  while (std::getline(in, line))
    {

    if (count == 2)
    {
      // iterating over each character in the string
      std::string base_cmp;
      int count_A=0, count_G=0, count_T=0, count_C=0, count_N=0;
      for (std::string::iterator it = line.begin(); it != line.end(); ++it)
      {
        base_cmp.clear();
        base_cmp.push_back(*it);
        if (base_cmp.compare("A") ==0) { count_A +=1;}
        else if (base_cmp.compare("T")==0) { count_T +=1;}
        else if (base_cmp.compare("G")==0) { count_G +=1;}
        else if (base_cmp.compare("C")==0) { count_C +=1;}
        else if (base_cmp.compare("N")==0) { count_N +=1;}
      }
      double gc_percent = static_cast<double>(count_C+count_G)/static_cast<double>(count_A+count_T+count_G+count_C+count_N);
      gc_percent_per_read.push_back(gc_percent);
    }

    if (count == 4)
    {
      // Reset counter
      count = 1;
    }
    else
    {
      count++;
    }
  }
  //Cleanup
  in.close();

  return wrap(gc_percent_per_read);
}


//' calculate Over Rep seqs
//'
//' Calculate sequece counts for each unique sequence and create a table with unique sequences and corresponding counts
//' @param infile A string giving the path for the fastqfile
//' @param out_prefix A string giving the prefix to be used for outputs
//' @param min_size An int for thhresholding over representation
//' @param buffer_size An int for the number of lines to keep in memory
//' @return calculate overrepresented sequence count
//' @export
// [[Rcpp::export]]
std::map<std::string,int> calc_over_rep_seq (std::string infile, std::string out_prefix,
                                             int min_size=5, int buffer_size = 1000000)
{
  std::map<std::string, int> over_rep_map;
  std::map<std::string,int>::iterator it;
  std::map<int,std::vector<int> > qual_by_column;
  std::map<int,std::vector<int> >::iterator qual_by_col_it;

  /*
   TODO: Check to see if this needs to be written to file instead

  std::string over_rep_out = out_prefix + ".over_rep.csv";
  std::ofstream over_rep_file;
  over_rep_file.open(over_rep_out.c_str());*/

  gz::igzstream in(infile.c_str());
  std::string line;
  int count = 1, line_count =1;
  while (std::getline(in, line))
    {

    if (count == 2)
    {
      it = over_rep_map.find(line);
      if (it != over_rep_map.end())
      {
        // if found increment by 1
        over_rep_map.at(line) += 1;
      }
      else
      {
        // if not found add new key and initialize to 1
        over_rep_map.insert(std::pair<std::string, int>(line,1));
      }
    }
    if (count == 4)
    {
      // reset count
      count = 1;
    }
    else
    {
      count++;
      }
    // Reduce map after 1^e6 reads
    if ( line_count % buffer_size ==0)
    {
      it=over_rep_map.begin();
      while(it != over_rep_map.end())
      {
        if (it->second <= min_size)
        {
          std::map<std::string, int>::iterator erase_it = it;
          it++;
          over_rep_map.erase(erase_it);
        }
        else
        {
          ++it;
        }

      }
    }
   line_count++;
  }
  //Cleanup
  in.close();


  /*TODO

  over_rep_file.close();
  */

  return over_rep_map ;
}



