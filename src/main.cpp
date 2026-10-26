#include "chain/fasta.hpp"
#include "chain/chainer.hpp"
#include "chain/paf.hpp"
#include <iostream>
#include <string>
#include <vector>

void print_usage() {
    std::cerr << "Usage: chain [options] -r <ref.fa> -q <query.fa>\n\n"
              << "Options:\n"
              << "  -r, --ref <file>       Target reference FASTA file (required)\n"
              << "  -q, --query <file>     Query sequence FASTA file (required)\n"
              << "  -k <int>               Minimizer k-mer size (default: 15)\n"
              << "  -w <int>               Minimizer window size (default: 10)\n"
              << "  -g, --max-gap <int>    Max gap distance between anchors (default: 5000)\n"
              << "  -c, --circular         Treat reference as circular chromosome\n"
              << "  --fenwick              Use O(N log N) Fenwick RMQ tree chaining\n"
              << "  -o, --output <file>    Output PAF file (default: stdout)\n";
}

int main(int argc, char* argv[]) {
    std::string ref_path, query_path, out_path;
    chain::ChainParams params;

    bool primary_only = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-r" || arg == "--ref") && i + 1 < argc) {
            ref_path = argv[++i];
        } else if ((arg == "-q" || arg == "--query") && i + 1 < argc) {
            query_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            out_path = argv[++i];
        } else if (arg == "-k" && i + 1 < argc) {
            params.k = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "-w" && i + 1 < argc) {
            params.w = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if ((arg == "-g" || arg == "--max-gap") && i + 1 < argc) {
            params.max_gap = std::stoll(argv[++i]);
        } else if (arg == "-c" || arg == "--circular") {
            params.circular = true;
        } else if (arg == "--fenwick") {
            params.use_fenwick = true;
        } else if (arg == "-p" || arg == "--primary-only") {
            primary_only = true;
        }
    }

    if (ref_path.empty() || query_path.empty()) {
        print_usage();
        return 1;
    }

    auto ref_records = chain::read_fasta(ref_path);
    auto query_records = chain::read_fasta(query_path);

    if (ref_records.empty() || query_records.empty()) {
        std::cerr << "Error: Reference or query FASTA has no sequences.\n";
        return 1;
    }

    std::ostream* out_stream = &std::cout;
    std::ofstream file_out;
    if (!out_path.empty()) {
        file_out.open(out_path);
        if (!file_out.is_open()) {
            std::cerr << "Error: Could not open output PAF file: " << out_path << "\n";
            return 1;
        }
        out_stream = &file_out;
    }

    chain::Chainer chainer(params);

    for (const auto& q : query_records) {
        for (const auto& r : ref_records) {
            auto chains = chainer.align(r.seq, q.seq);

            for (size_t i = 0; i < chains.size(); ++i) {
                bool is_primary = (i == 0);
                if (primary_only && !is_primary) break;
                std::string paf_line = chain::format_paf_record(q.name, q.seq.size(),
                                                                r.name, r.seq.size(),
                                                                chains[i], is_primary);
                *out_stream << paf_line << "\n";
            }
        }
    }

    return 0;
}
