#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    ArgumentParser parser("Conterminator_Dust");
    addOption(parser, ArgParseOption("F", "fasta", "Path to the input file", ArgParseArgument::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("R", "rna-file", "Are these RNA sequences? Will throw out mRna completely that come from contaminated sequences (this flag has to be set otherwise you get segfaults because the alignment positions are w.r.t. the entire seq, not the mRna."));
    //addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "Triplets (csv) with contamination (SEQ ID, from, to)"));
    addOption(parser, ArgParseOption("C", "contamination", "Contamination file (triplets with comma separated values, one triplet per line)", ArgParseArgument::INPUT_FILE, "IN"));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    bool rna_file = isSet(parser, "rna-file");;
    CharString path_in, triplet_path; //, text;
    getOptionValue(path_in, parser, "fasta");
    getOptionValue(triplet_path, parser, "contamination");
    
    std::vector<std::tuple<std::string, uint64_t, uint64_t>> triplets;
    
    std::ifstream infile(toCString(triplet_path));
    uint64_t from, to;
    std::string identifier;
    while (infile >> identifier >> from >> to)
    {
        triplets.push_back({std::move(identifier), from, to});
    }

    CharString path_out(path_in);
    path_out += ".tmp.fa";

    SeqFileIn seqFileIn(toCString(path_in));
    SeqFileOut seqFileOut(toCString(path_out));

    SequenceOutputOptions options(80, false);

    while (!atEnd(seqFileIn))
    {
        CharString _id;
        IupacString seq;

        readRecord(_id, seq, seqFileIn);

        std::string id(toCString(_id));

		// >lcl|ABGB01001532.1_mrna_569 [locus_tag=EBI_27009] [product=hypothetical protein] [location=<1..487] [gbkey=mRNA]
        auto iter = std::find_if(triplets.begin(), triplets.end(), [&id](auto const & x) {
		auto const & id_from_triplets = std::get<0>(x);
		auto hit = id.find(id_from_triplets);
		return hit != std::string::npos && (hit + length(id_from_triplets) == length(id) || !std::isdigit(id[hit + length(id_from_triplets)])); // next symbol can't be a number (NC_001.1 vs NC_001.10)
	});

        if (iter != triplets.end())
        {
            if (rna_file)
                continue;

            uint64_t from = std::get<1>(*iter);
            uint64_t to = std::get<2>(*iter);

            //seqan3::debug_stream << from << " - " << to << '\n';
            //seqan3::debug_stream << seq << '\n';

            // find flanking N's
            while (from > 0 && seq[from] != Iupac('N'))
                --from;

            while (to < length(seq) && seq[to] != Iupac('N'))
                ++to;

            float frac = static_cast<float>(to - from)/length(seq);
            if (frac > 50 || length(seq) < 50000)
                continue;

            //seqan3::debug_stream << from << " - " << to << '\n';

            for (uint64_t i = from; i < to; ++i)
                seq[i] = Iupac('N');

            //seqan3::debug_stream << seq << '\n' << '\n';

            //id = toCString(_id);
        }

        writeRecord(seqFileOut, id, seq);
    }
}
