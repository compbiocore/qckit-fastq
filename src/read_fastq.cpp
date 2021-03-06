#include <Rcpp.h>
#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdint>
#include "gzstream.h"
#include "zlib.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' Gets quality score encoding format from the FASTQ file. Return possibilities are Sanger(/Illumina1.8),
//' Solexa(/Illumina1.0), Illumina1.3, and Illumina1.5. This encoding is heuristic based and may not be 100% accurate
//' since there is overlap in the encodings used, so it is best if you already know the format.
//'
//' @param infile  A string giving the path for the fastq file
//' @param reads_used int, the number of reads to use to determine the encoding format.
//' @examples
//' infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
//' find_format(infile,100)
//' @return A string denoting the read format. Possibilities are Sanger, Solexa, Illumina1.3, and Illumina1.5.
//' @export
// [[Rcpp::export]]
std::string find_format(std::string infile, int reads_used) {
  // TODO: test for if buffer_size & reads_used exceeds number of reads in file
    std::vector<int> qual_scores;
    std::vector<int>::iterator it;
    gz::igzstream in(infile.c_str());
    std::string line, score_format;
    int counter = 1;
    int line_count = 1;

    std::vector<double> gc_percent_all;
    while (std::getline(in, line)) {
        if (counter > reads_used) {
            break;
        } else {
            if (line_count == 4) {
                for (char &c : line) {
                    qual_scores.push_back(int(c));
                }

                /* alternative implementation basd on iterators
                 * for (std::string::iterator chit = line.begin(); chit != line.end(); ++chit)
                {
                    qual_scores.push_back(int(*chit));
                }*/

                line_count = 1;
                counter++;
            } else {
                line_count++;
            }
        }
    }
    int max_score = *max_element(qual_scores.begin(),qual_scores.end());
    int min_score = *min_element(qual_scores.begin(),qual_scores.end());

    // Based off of https://en.wikipedia.org/wiki/FASTQ_format#Encoding
    if ((min_score > 32) & (min_score < 59) & (max_score < 127)) { score_format = "Sanger"; } // Sanger is ASCII 33 to 126
    else if ((min_score > 58) & (min_score < 64) & (max_score < 127)) { score_format = "Solexa"; } // Solexa is 59 to 126
    else if ((min_score > 63) & (min_score < 66) & (max_score < 127)) { score_format = "Illumina1.3"; } // 64 to 126
    // Illumina1.5: 0/64 and 1/65 are not in use, and 2/66 is only used if read end segment is all Q15 or below
    else if ((min_score > 65) & (max_score < 127)) { score_format = "Illumina1.5"; }
    else {  throw "No plausible encoding format! Check FASTQ reads to make sure quality scores are >32 and < 127."; }

    return score_format;
}

// [[Rcpp::plugins(cpp11)]]

//' Calculate score based on Illumina format
//'
//' @param score  An ascii quality score from the fastq
//' @param score_format The illumina format
//' @examples
//' calc_format_score("A","Sanger")
//' @return a string as with the best guess as to the illumina format
//' @export
// [[Rcpp::export]]
int calc_format_score(char score, std::string score_format)
{
    // Based off of https://en.wikipedia.org/wiki/FASTQ_format#Encoding
    int calc_score = 0;
    if(score_format == "Sanger") calc_score = int(score) - 33;
    else if(score_format == "Solexa") calc_score = int(score) - 64;
    else if(score_format=="Illumina1.3") calc_score = int(score) - 64;
    else if(score_format=="Illumina1.5") calc_score = int(score) - 64;

    return calc_score;
}

// Calculate summary of quality scores over position
//
// Description
// @param inmat A matrix of score vectors per position
std::vector <std::vector<int> > qual_score_per_position(const std::map<int, std::vector<uint8_t> > &inmat) {
    std::vector <std::vector<int> > qual_score_mat_results;
    std::vector<int> q_10, q_25, q_50, q_75, q_90;

    std::map < int, std::vector < uint8_t > > ::const_iterator
    mat_it = inmat.begin();

    for (mat_it = inmat.begin(); mat_it != inmat.end(); mat_it++) {
        std::vector <uint8_t> quantile = mat_it->second;

        int Q10 = static_cast<int> (quantile.size() * 0.1);
        int Q25 = (quantile.size() + 1) / 4;
        int Q50 = (quantile.size() + 1) / 2;
        int Q75 = Q25 + Q50;
        int Q90 = static_cast<int> (quantile.size() * 0.9);

        std::nth_element(quantile.begin(), quantile.begin() + Q10, quantile.end());
        q_10.push_back(static_cast<int>(quantile[Q10]));
        quantile.clear();

        std::nth_element(quantile.begin(), quantile.begin() + Q25, quantile.end());
        q_25.push_back(static_cast<int>(quantile[Q25]));
        quantile.clear();

        quantile = mat_it->second;
        std::nth_element(quantile.begin(), quantile.begin() + Q50, quantile.end());
        q_50.push_back(static_cast<int>(quantile[Q50]));
        quantile.clear();

        quantile = mat_it->second;
        std::nth_element(quantile.begin(), quantile.begin() + Q75, quantile.end());
        q_75.push_back(static_cast<int>(quantile[Q75]));

        std::nth_element(quantile.begin(), quantile.begin() + Q90, quantile.end());
        q_90.push_back(static_cast<int>(quantile[Q90]));
        quantile.clear();
    }
    qual_score_mat_results.push_back(q_10);
    qual_score_mat_results.push_back(q_25);
    qual_score_mat_results.push_back(q_50);
    qual_score_mat_results.push_back(q_75);
    qual_score_mat_results.push_back(q_90);
    return qual_score_mat_results;
}

//' Calculate the mean quality score per read of the FASTQ gzipped file
//' 
//' @param infile A string giving the path for the fastqfile
//' @examples
//' infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
//' qual_score_per_read(infile)$q50_per_position[1:10]
//' @return mean quality per read
//' @export
// [[Rcpp::export]]
Rcpp::List qual_score_per_read(std::string infile) {

    std::vector<double> quality_score_per_read;
    std::vector <uint8_t> qual_by_column;
    std::vector<uint8_t>::iterator qual_by_col_it;

    std::map<int, std::vector<uint8_t> > qual_score_matrix;

    gz::igzstream in(infile.c_str());
    std::string line;
    int count = 1;
    double quality_score_mean = 0;
    std::string score_format = find_format(infile, 10000);
    while (std::getline(in, line))
    {

        if (count == 4)
        {
            // iterate over each value for quality
            qual_by_column.clear();
            int pos_counter = 1;
            int calc_score =0;
            for (std::string::iterator it = line.begin(); it != line.end(); ++it)
            {
                calc_score = calc_format_score(*it, score_format);

                qual_by_column.push_back(calc_score);
                if (pos_counter <= qual_score_matrix.size()) {
                    qual_score_matrix[pos_counter].push_back(calc_score);
                } else {
                    std::vector <uint8_t> tmp_qual;
                    tmp_qual.push_back(calc_score);
                    std::pair < int, std::vector < uint8_t >> tmp = std::pair < int, std::vector <
                                                                                     uint8_t >> (pos_counter, tmp_qual);
                    qual_score_matrix.insert(tmp);
                }
                pos_counter++;
            }
            quality_score_mean = static_cast<double>(std::accumulate(qual_by_column.begin(), qual_by_column.end(),
                                                                     0.0)) / static_cast<double>(qual_by_column.size());
            quality_score_per_read.push_back(quality_score_mean);
            count = 1;
        } else {
            count++;
        }
    }
    std::vector <std::vector<int>> qual_score_summary_by_position;
    qual_score_summary_by_position = qual_score_per_position(qual_score_matrix);

    std::vector<double> mu_per_position;
    std::map < int, std::vector < uint8_t > > ::iterator
    mat_it = qual_score_matrix.begin();

    for (mat_it = qual_score_matrix.begin(); mat_it != qual_score_matrix.end(); mat_it++) {
        std::vector <uint8_t> quantile = mat_it->second;

        mu_per_position.push_back(static_cast<double>(std::accumulate(quantile.begin(),
                                                                      quantile.end(), 0.0)) /
                                  static_cast<double>(quantile.size()));
    }
    //Cleanup
    in.close();
    return Rcpp::List::create(Rcpp::Named("mu_per_read") = quality_score_per_read,
                              Rcpp::Named("mu_per_position") = mu_per_position,
                              Rcpp::Named("q10_per_position") = qual_score_summary_by_position[0],
                              Rcpp::Named("q25_per_position") = qual_score_summary_by_position[1],
                              Rcpp::Named("q50_per_position") = qual_score_summary_by_position[2],
                              Rcpp::Named("q75_per_position") = qual_score_summary_by_position[3],
                              Rcpp::Named("q90_per_position") = qual_score_summary_by_position[4]);
}


//' Calculate GC nucleotide sequence content per read of the FASTQ gzipped file
//' @param infile A string giving the path for the fastqfile
//' @examples
//' infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
//' gc_per_read(infile)[1:10]
//' @return GC content perncentage per read
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gc_per_read(std::string infile) {

    std::map<int, std::vector<int> > qual_by_column;
    std::map < int, std::vector < int > > ::iterator
    qual_by_col_it;


    gz::igzstream in(infile.c_str());
    std::string line;
    int count = 1;
    //std::vector<int,std::vector<int> > base_counts;
    std::vector<double> gc_percent_per_read;
    while (std::getline(in, line)) {

        if (count == 2) {
            // iterating over each character in the string
            std::string base_cmp;
            int count_A = 0, count_G = 0, count_T = 0, count_C = 0, count_N = 0;
            for (std::string::iterator it = line.begin(); it != line.end(); ++it) {
                base_cmp.clear();
                base_cmp.push_back(*it);
                if (base_cmp.compare("A") == 0) { count_A += 1; }
                else if (base_cmp.compare("T") == 0) { count_T += 1; }
                else if (base_cmp.compare("G") == 0) { count_G += 1; }
                else if (base_cmp.compare("C") == 0) { count_C += 1; }
                else if (base_cmp.compare("N") == 0) { count_N += 1; }
            }
            double gc_percent = static_cast<double>(count_C + count_G) /
                                static_cast<double>(count_A + count_T + count_G + count_C + count_N);
            gc_percent_per_read.push_back(gc_percent);
        }

        if (count == 4) {
            // Reset counter
            count = 1;
        } else {
            count++;
        }
    }
    //Cleanup
    in.close();

    return wrap(gc_percent_per_read);
}


//' Calculate sequece counts for each unique sequence and create a table with unique sequences
//' and corresponding counts
//' @param infile A string giving the path for the fastqfile
//' @param min_size An int for thhresholding over representation
//' @param buffer_size An int for the number of lines to keep in memory
//' @return calculate overrepresented sequence count
//' @examples
//' infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
//' calc_over_rep_seq(infile)[seq_len(5)]
//' @export
// [[Rcpp::export]]
std::map<std::string, int> calc_over_rep_seq(std::string infile,
                                             int min_size = 5, int buffer_size = 1000000) {
    std::map<std::string, int> over_rep_map;
    std::map<std::string, int>::iterator it;
    std::map<int, std::vector<int> > qual_by_column;
    std::map < int, std::vector < int > > ::iterator
    qual_by_col_it;

    /*
     TODO: Check to see if this needs to be written to file instead

    std::string over_rep_out = out_prefix + ".over_rep.csv";
    std::ofstream over_rep_file;
    over_rep_file.open(over_rep_out.c_str());*/

    gz::igzstream in(infile.c_str());
    std::string line;
    int count = 1, line_count = 1;
    while (std::getline(in, line)) {

        if (count == 2) {
            it = over_rep_map.find(line);
            if (it != over_rep_map.end()) {
                // if found increment by 1
                over_rep_map.at(line) += 1;
            } else {
                // if not found add new key and initialize to 1
                over_rep_map.insert(std::pair<std::string, int>(line, 1));
            }
        }
        if (count == 4) {
            // reset count
            count = 1;
        } else {
            count++;
        }
        // Reduce map after 1^e6 reads
        if (line_count % buffer_size == 0) {
            it = over_rep_map.begin();
            while (it != over_rep_map.end()) {
                if (it->second <= min_size) {
                    std::map<std::string, int>::iterator erase_it = it;
                    it++;
                    over_rep_map.erase(erase_it);
                } else {
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

    return over_rep_map;
}