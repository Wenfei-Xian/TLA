#include <iostream>
#include <fstream>
#include <vector>
#include <htslib/sam.h>
#include "spoa/spoa.hpp"
#include <string>
#include <sstream>
#include <algorithm>
#include <unistd.h>

//g++ -o TLA.v1 TLA.v1.cpp -I/tmp/global2/wxian/software/htslib/include -L/tmp/global2/wxian/software/htslib/lib -lhts -I/tmp/global2/wxian/software/spoa/include -L/tmp/global2/wxian/software/spoa/build/lib -lspoa

int usage(){
	std::cout << "TLA: Telomere local assembly" << std::endl;
	std::cout << "Usage: TLA [options]" << std::endl;
	std::cout << "options:" << std::endl;
	std::cout << "-b string     alignment file in bam format" << std::endl;
	std::cout << "-m int        minimum length of contigs will be used (default:500kb)" << std::endl;
	std::cout << "-d int        the distance to the start or end of contigs, used to extract reads in the BAM file (default:10kb)" << std::endl;
	//std::cout << "-l int        the length to trim from the far left or far right end of the reads (default:15kb)" << std::endl;
	std::cout << "-s int        minimum length of perfect tandem repeats region(default:100bp)" << std::endl;
	std::cout << "-p int        minimum percentage of reads contained perfect tandem repeats (default:0.3)" << std::endl;
	std::cout << "-c int        coverage threshold for the generation of consensus sequences (default:2)" << std::endl;
	return 0;
}

template <typename Container>
int length(const Container& c) {
    return static_cast<int>(c.size());
}

int detect_longest_tandem_repeat( std::string& DNA, const int cutoffunit_p, const int unitlen_p) {
	unsigned int DNA_len = DNA.length();

	std::transform(DNA.begin(), DNA.end(), DNA.begin(), ::toupper);

	int max_repeat_length = 0;

	for (unsigned int start = 0; start < DNA_len; ++start) {
		for (short ssr_len = 1; ssr_len <= unitlen_p; ++ssr_len) {
			short repeat = 1;
			for (;; ++repeat) {
				short match = 0;
				for (short base = 0; base < ssr_len; ++base) {
					char left = DNA[start + base];
					char right = DNA[start + ssr_len * repeat + base];
					if (left != right || left == 'N') {
						match = 1;
						break;
					}
				}
				if (match == 1) {
					break;
				}
			}
			if( ssr_len * repeat >= cutoffunit_p ){
				start+=ssr_len*repeat-1;
				if( ssr_len >=4 ){
					max_repeat_length = std::max(max_repeat_length, ssr_len * repeat);
				}
				break;
			}
		}
	}

	return max_repeat_length;

}

struct Read {
    std::string name;
    std::string sequence;
};

std::string bam_get_seq_string(const bam1_t* b) {
	uint8_t* seq = bam_get_seq(b);
	std::string result;
	for (int i = 0; i < b->core.l_qseq; ++i) {
		result.push_back(seq_nt16_str[bam_seqi(seq, i)]);
	}
	return result;
}

int main(int argc, char** argv) {
	if(argc ==  1 ){
		return usage();
        }
	else if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
		return usage();
	}

	const char* bam_file=nullptr;
	size_t minimum_contig_length=500000;
	size_t distance_to_start_or_end=10000;
	//int length_retain=15000;
	size_t minimum_tandem_repeat_length=100;
	size_t minimum_perfect_tandem_repeat=0.3;
	size_t coverage=2;
	int a;

	while(( a=getopt( argc, argv, "b:m:d:l:s:p:c:")) >= 0 ){

		if( a == 'b' ){
			bam_file=optarg;
                }
		else if( a == 'm' ){
			minimum_contig_length=atoi(optarg);
		}
		else if( a == 'd'){
			distance_to_start_or_end=atoi(optarg);
		}
		//else if( a == 'l' ){
		//	length_retain=atoi(optarg);
		//}
		else if( a == 's' ){
			minimum_tandem_repeat_length=atoi(optarg);
		}
		else if( a == 'p' ){
			minimum_perfect_tandem_repeat=atof(optarg);
		}
		else if( a == 'c' ){
			coverage=atoi(optarg);
		}
	}

	//int minimum_contig_length = std::stoi(argv[2]);
	std::string distance_to_start_or_end_str = std::to_string(distance_to_start_or_end);
	// Open BAM for reading
	samFile* in = sam_open( bam_file, "r");
	if (!in) {
		std::cerr << "Cannot open BAM file: " << bam_file << std::endl;
		return 1;
	}

	bam_hdr_t* header = sam_hdr_read(in);
	hts_idx_t* idx = sam_index_load(in, bam_file);
	if (!idx) {
		std::cerr << "Cannot load index for BAM file: " << bam_file << std::endl;
		return 1;
	}

	bam1_t* aln = bam_init1();

	for (int i = 0; i < header->n_targets; ++i) {
		if (header->target_len[i] > minimum_contig_length) {
			std::string ref_name = header->target_name[i];

			// Collect reads for SPOA
			std::vector<std::string> start_reads, end_reads;
			// Collect reads and name for fasta output
			std::vector<Read> start_reads_name, end_reads_name;

			// Fetch reads from the start (0-5K)
			size_t long_tandem_count_start = 0;
			hts_itr_t* iter = sam_itr_queryi(idx, i, 0, distance_to_start_or_end);
			while (sam_itr_next(in, iter, aln) >= 0) {
				std::string seq = bam_get_seq_string(aln);

				Read currentRead;
				currentRead.name = bam_get_qname(aln);
				currentRead.sequence = bam_get_seq_string(aln);
				size_t len_to_extract = std::min<size_t>(minimum_contig_length, currentRead.sequence.length());
				currentRead.sequence = currentRead.sequence.substr(0, len_to_extract);
				
				//size_t len_to_extract = std::min<size_t>(minimum_contig_length, seq.length());
				std::string DNA = seq.substr(0, len_to_extract);
				
				size_t longest_tandem = detect_longest_tandem_repeat(DNA, 50, 10);
				if(longest_tandem > minimum_tandem_repeat_length) {
					start_reads.push_back(DNA);
					start_reads_name.push_back(currentRead);
					long_tandem_count_start++;
				}
				//std::cout << "Longest tandem repeat for start sequence of " << bam_get_qname(aln) << ": " << longest_tandem << "\n";
			}
			hts_itr_destroy(iter);

			if( long_tandem_count_start >= start_reads.size()*minimum_perfect_tandem_repeat ){
				std::ofstream start_out(ref_name + ".start." + distance_to_start_or_end_str + ".fa");
				for (const auto& read : start_reads_name ) {
					start_out << ">" << read.name << "\n";
					start_out << read.sequence << "\n";
				}
				std::cout << ref_name << "_start:telomere" << std::endl;
			}

			// Fetch reads from the end (-5K)
			size_t long_tandem_count_end = 0;
			size_t end_start = header->target_len[i] - distance_to_start_or_end;
			iter = sam_itr_queryi(idx, i, end_start, header->target_len[i]);
			while (sam_itr_next(in, iter, aln) >= 0) {
				std::string seq = bam_get_seq_string(aln);

				Read currentRead;
				currentRead.name = bam_get_qname(aln);
				currentRead.sequence = bam_get_seq_string(aln);
				size_t len_to_extract = std::min<size_t>(minimum_contig_length, currentRead.sequence.length());
				currentRead.sequence = currentRead.sequence.substr(0, len_to_extract);

				//size_t len_to_extract = std::min<size_t>(minimum_contig_length, seq.length());
				std::string DNA = seq.substr(seq.length() - len_to_extract);
				
				size_t longest_tandem = detect_longest_tandem_repeat(DNA, 50, 10);
				if(longest_tandem > minimum_tandem_repeat_length) {
					end_reads.push_back(DNA);
					end_reads_name.push_back(currentRead);
					long_tandem_count_end++;
				}
				//std::cout << "Longest tandem repeat for end sequence of " << bam_get_qname(aln) << ": " << longest_tandem << "\n";
			}
			hts_itr_destroy(iter);

			if( long_tandem_count_end >= end_reads.size()*minimum_perfect_tandem_repeat ){
				std::ofstream end_out(ref_name + ".end." + distance_to_start_or_end_str + ".fa");
				for (const auto& read : end_reads_name) {
					end_out << ">" << read.name << "\n";
					end_out << read.sequence << "\n";
				}
				std::cout << ref_name << "_end:telomere" << std::endl;
			}

			// Create SPOA alignment engine
			auto alignment_engine_repeat = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 7,-7,-7,-2,-9,-4 );//suitable for telomere
			//auto alignment_engine_normal = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 5,-4,-8,-6,-10,-4 );//suitable for normal
			spoa::Graph graph{};

			for (const auto& read : start_reads) {
				if( long_tandem_count_start >= start_reads.size()*minimum_perfect_tandem_repeat ){
					auto alignment = alignment_engine_repeat->Align(read, graph);
					graph.AddAlignment(alignment, read);
				}
			}
			if( long_tandem_count_start >= start_reads.size()*minimum_perfect_tandem_repeat ){
				std::string consensus = graph.GenerateConsensus(coverage);

				// Output the consensus to a file
				std::ofstream cons_out(ref_name + ".start." + distance_to_start_or_end_str + ".cons.fa");
				cons_out << ">" << ref_name << ".start." + distance_to_start_or_end_str + ".consensus\n";
				cons_out << consensus << "\n";
				cons_out.close();
			}

			// Do the same for end reads, but reset the graph first
			graph = spoa::Graph{};
			for (const auto& read : end_reads) {
				if( long_tandem_count_end >= end_reads.size()*minimum_perfect_tandem_repeat ){
					auto alignment = alignment_engine_repeat->Align(read, graph);
					graph.AddAlignment(alignment, read);
				}
			}
			if( long_tandem_count_end >= end_reads.size()*minimum_perfect_tandem_repeat ){
				std::string consensus = graph.GenerateConsensus(coverage);

				// Output the consensus to a file
				std::ofstream cons_out(ref_name + ".end." + distance_to_start_or_end_str + ".cons.fa");
				cons_out << ">" << ref_name << ".end." + distance_to_start_or_end_str + ".consensus\n";
				cons_out << consensus << "\n";
				cons_out.close();
			}
		}
	}
	bam_destroy1(aln);
	bam_hdr_destroy(header);
	sam_close(in);
	return 0;
}
