#include <iostream>
#include <vector>
#include <filesystem>
#include <algorithm>

#include "omp.h"

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "include/KMC/kmc_api/kmc_file.h"

using namespace seqan3;
using seqan3::operator""_dna5;

struct LogWriter
{
    std::ofstream log_file;

    uint16_t dir_prefix_length = 0;

    LogWriter(std::filesystem::path const & log_path, std::filesystem::path const & fasta_dir_path)
    {
        if (log_path != std::filesystem::path())
        {
            log_file.open(log_path);
            // print header
            log_file << "Path" << '\t' << "SEQ ID" << '\t' << "Bases replaced (= #kmers matched)" << '\t' << "SEQ length" << '\t' << "Note" << '\n';
            // compute common prefix length in path

            std::string const fasta_dir_path_string = fasta_dir_path.string();
            dir_prefix_length = fasta_dir_path_string.length();
            if (fasta_dir_path_string.back() != '/')
                ++dir_prefix_length;
        }
    }

    ~LogWriter()
    {
        if (log_file.is_open())
            log_file.close();
    }

    void write(std::filesystem::path const & p, auto & id, auto nbr_bases_changed, auto seq_length, auto const & note)
    {
        if (log_file.is_open())
        {
            log_file << p.string().substr(dir_prefix_length) << '\t'
                     << id << '\t'
                     << nbr_bases_changed << '\t'
                     << seq_length << '\t'
                     << note << '\n';
        }
    }
};

int main(int argc, char ** argv)
{
    argument_parser parser{"Contamination-Cleaner", argc, argv};
    parser.info.author = "Christopher Pockrandt";
    parser.info.short_description = "k-mer based approach to remove (human) contamination from fasta files. Contaminated based will be replaced with N's inplace.";
    parser.info.version = "1.0.0";

    std::filesystem::path path_fasta_files{};
    std::filesystem::path path_kmc_database{};
    std::filesystem::path path_log_file{};
    float threshold_remove_seq = 0.1;
    float threshold_remove_file = 0.1;
    unsigned threads = omp_get_max_threads();

    std::vector<std::string> const extensions = {".fna", ".fa", ".fasta", ".fas", ".fsa"};

    parser.add_positional_option(path_fasta_files, "Directory to fasta files.");
    parser.add_positional_option(path_kmc_database, "Path to KMC database.");

    parser.add_option(path_log_file, 'l', "path-log-file", "Stores a tab separated file with information which sequences and files were altered or deleted.", option_spec::DEFAULT); // TODO: make it optional in code
    parser.add_option(threshold_remove_seq , 'a', "remove-seq-threshold", "Percentage of bases replaced with Ns, to discard entire sequence.", option_spec::DEFAULT, arithmetic_range_validator{0.0001, 1.0});
    parser.add_option(threshold_remove_file, 'b', "remove-file-threshold", "Percentage of sequence discarded, to discard entire file.", option_spec::DEFAULT, arithmetic_range_validator{0.0001, 1.0});
    parser.add_option(threads , 't', "threads", "Number of threads.", option_spec::DEFAULT);

    try
    {
        parser.parse();
    }
    catch (argument_parser_error const & ext)
    {
        std::cout << "Error while parsing command line parameters: " << ext.what() << '\n';
        return -1;
    }

    // Store fasta file paths in container (for parallelization and progress output)
    std::cout << "Searching for fasta files ... " << std::flush;
    std::vector<std::filesystem::path> fasta_files;
    for (auto & p: std::filesystem::recursive_directory_iterator(path_fasta_files))
    {
        if (std::find(extensions.begin(), extensions.end(), p.path().extension()) != extensions.end())
        {
            fasta_files.push_back(p);
        }
    }
    std::cout << fasta_files.size() << " found." << std::endl;

    // Load KMC file
    CKMCFile kmc_db;
    if (!kmc_db.OpenForRA(path_kmc_database))
        std::cerr << "Could not open KMC database file.\n";
    else
        std::cout << "KMC database loaded." << std::endl;

    LogWriter logs(path_log_file, path_fasta_files);

    std::cout << "Processing file 0 of " << fasta_files.size() << " (0.00 %)" << std::flush;

    // iterate over all fasta files
    uint64_t bases_processed = 0;
    uint32_t fasta_id = 0;
    std::vector<std::vector<uint32_t>> kmer_counts_array(threads);

    #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (auto it = fasta_files.begin(); it != fasta_files.end(); ++it)
    {
        auto & p = *it;
        auto & kmer_counts = kmer_counts_array[omp_get_thread_num()];

        // print progress
        #pragma omp critical(increment_fasta_id)
        ++fasta_id;
        float const progress = 1.0f * fasta_id / fasta_files.size();

        #pragma omp critical(progress_output)
        {
            std::cout << "\rProcessing file " << fasta_id << " of " << fasta_files.size()
                      << " (" << (truncf(progress*10000)/100) << " %), ";

            if (bases_processed < 128*1024*1024) // < ~0.1 Gbases, print in Mbases
                std::cout << (truncf(1.0f * bases_processed / (1024*1024)*100)/100) << " Mbases processed in total";
            else
                std::cout << (truncf(1.0f * bases_processed / (1024*1024*1024)*100)/100) << " Gbases processed in total";

            std::cout << "\x1b[K" << std::flush; // \e[K - clr_eol (remove anything after the cursor)
        }

        uint32_t nbr_sequences = 0;
        uint32_t nbr_sequences_written = 0;

        std::filesystem::path tmp_file = p;
        tmp_file += ".tmp.fa";

        // load file and iterate over sequences
        { // make sure destructors of file classes are called before files are deleted or moved

            sequence_file_input<sequence_file_input_default_traits_dna> fin{p};
            sequence_file_output fout{tmp_file};
            for (auto & rec : fin)
            {
                ++nbr_sequences;
                uint32_t nbr_bases_changed = 0;

                auto & id = seqan3::get<seqan3::field::id>(rec);
                auto & seq = seqan3::get<seqan3::field::seq>(rec);

                auto seq_view = seq | views::to_char;
                std::string seq_string(seq_view.begin(), seq_view.end());

                // retrieve k-mer counts
                kmc_db.GetCountersForRead(seq_string, kmer_counts);
                // debug_stream << kmer_counts << '\n';

                // replace first base with an N if the k-mer occurs in the database
                // debug_stream << '\n' << seq << '\n';
                for (uint32_t i = 0; i < kmer_counts.size(); ++i)
                {
                    if (kmer_counts[i] > 0)
                    {
                        ++nbr_bases_changed;
                        seq[i] = 'N'_dna5;
                    }
                }
                // debug_stream << seq << '\n';

                float const fraction_bases_changed = 1.0f * nbr_bases_changed / seq.size();
                if (fraction_bases_changed < threshold_remove_seq)
                {
                    ++nbr_sequences_written;
                    fout.emplace_back(seq, id);
                    #pragma omp critical(writing_logs)
                    {
                        if (nbr_bases_changed > 0)
                            logs.write(p, id, nbr_bases_changed, seq.size(), ".");
                    }
                }
                else
                {
                    #pragma omp critical(writing_logs)
                    logs.write(p, id, nbr_bases_changed, seq.size(), "SEQ DELETED");
                }

                #pragma omp critical(increment_bases_processed)
                bases_processed += seq.size();
            }
        }

        std::filesystem::remove(p);

        float const fraction_sequences_discarded = 1 - (1.0f * nbr_sequences_written / nbr_sequences);
        if (nbr_sequences_written == 0 || fraction_sequences_discarded >= threshold_remove_file)
        {
            std::filesystem::remove(tmp_file);
            #pragma omp critical(writing_logs)
            logs.write(p, "", 0, 0, "FILE DELETED");
        }
        else
        {
            std::filesystem::rename(tmp_file, p);
        }
    }

    std::cout << '\n';

    return 0;
}
