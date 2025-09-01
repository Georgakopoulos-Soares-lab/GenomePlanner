// plan_one.cpp
// One-shot program: build FM-index from source FASTA (gz or plain),
// run DP planner w/ optional backtracking export, produce TSV outputs.
//
// Compile:
//   g++ -O3 -std=gnu++17 plan_one.cpp -lsdsl -ldivsufsort -ldivsufsort64 -lz -o plan_one
//
// Usage:
//   ./plan_one --source SRC.fasta[.gz] --target TGT.fasta[.gz] \
//              --W INT --pcr FLOAT --join FLOAT --synth FLOAT \
//              --summary-out out.tsv [--moves-out moves.tsv] [--progress N] [--keep-temp]
//
// Summary TSV columns (tab-separated):
//   target_file  chromosome  length  total_cost  blocks  replication_blocks  synthesis_blocks  joins
//   avg_replication_len  avg_synthesis_len  total_replication_cost  total_synthesis_cost  total_join_cost
//
// Moves TSV columns (tab-separated), if requested:
//   chromosome  block_id  start  end  length  op  acq_cost  join_applied  join_cost
//
// Notes:
//  - Progress prints to stderr as "DP progress: X% ..."
//  - Handles gz by inflating to a temp file (zlib) then cleaning up (unless --keep-temp).

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <zlib.h>

#include <sdsl/csa_wt.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/util.hpp>

namespace fs = std::filesystem;
using fm_index_t = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<256>>, 512, 1024>;

struct Args {
    std::string source;
    std::string target;
    int W = -1;
    double cost_pcr = -1.0;
    double cost_join = -1.0;
    double cost_synth = -1.0;
    std::string summary_out;      // required
    std::string moves_out;        // optional (if empty => no moves file)
    int progress = 5;             // percent step (0 to disable)
    bool keep_temp = false;
};

static void die(const std::string& msg) {
    std::cerr << "ERROR: " << msg << "\n";
    std::exit(2);
}

static void usage() {
    std::cout <<
R"(Usage:
  plan_one --source SRC.fasta[.gz] --target TGT.fasta[.gz]
           --W INT --pcr FLOAT --join FLOAT --synth FLOAT
           --summary-out PATH [--moves-out PATH]
           [--progress INT] [--keep-temp] [--help]

Required:
  --source         Source genome FASTA (plain or .gz)
  --target         Target genome FASTA (plain or .gz)
  --W              Maximum block size (int)
  --pcr            Fixed per-block cost when reusing (PCR) (double)
  --join           Fixed join cost between consecutive blocks (double)
  --synth          Per-base synthesis cost (double)
  --summary-out    TSV file path for per-chromosome stats and totals

Optional:
  --moves-out      TSV file path for full step-by-step plan (can be very large)
  --progress       Progress step in percent (default 5; 0 disables)
  --keep-temp      Keep temporary inflated FASTAs
  --help           Show this message
)";
}

static bool is_gzip_path(const std::string& p) {
    auto s = fs::path(p).extension().string();
    if (s == ".gz" || s == ".gzip") return true;
    // also handle .fa.gz / .fasta.gz double extensions
    auto fname = fs::path(p).filename().string();
    return fname.size() >= 3 && fname.rfind(".gz") == fname.size() - 3;
}

static std::string stem_until_first_dot(const std::string& path) {
    std::string base = fs::path(path).filename().string();
    auto pos = base.find('.');
    return (pos == std::string::npos) ? base : base.substr(0, pos);
}

static std::string gunzip_to_temp(const std::string& gz_path, const std::string& hint, std::vector<std::string>& cleanup) {
    fs::path tmpdir = fs::temp_directory_path();
    auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    fs::path out = tmpdir / (stem_until_first_dot(hint) + "_inflated_" + std::to_string(now) + ".fasta");

    gzFile gzf = gzopen(gz_path.c_str(), "rb");
    if (!gzf) die("cannot open gzip file: " + gz_path);
    std::ofstream ofs(out, std::ios::binary);
    if (!ofs.is_open()) die("cannot create temp file: " + out.string());

    char buf[1<<16];
    int r = 0;
    while ((r = gzread(gzf, buf, sizeof(buf))) > 0) {
        ofs.write(buf, r);
    }
    gzclose(gzf);
    ofs.close();

    cleanup.push_back(out.string());
    return out.string();
}

static std::string ensure_plain_fasta(const std::string& maybe_gz_path, std::vector<std::string>& cleanup) {
    if (is_gzip_path(maybe_gz_path)) {
        return gunzip_to_temp(maybe_gz_path, maybe_gz_path, cleanup);
    }
    return maybe_gz_path;
}

static std::map<std::string, std::string> read_fasta_and_clean(const std::string& path) {
    std::map<std::string, std::string> sequences;
    std::ifstream fasta_file(path);
    if (!fasta_file.is_open()) die("cannot open FASTA: " + path);
    std::string line, header, current_sequence;
    while (std::getline(fasta_file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!header.empty()) sequences[header] = current_sequence;
            header = line.substr(1);
            current_sequence.clear();
        } else {
            for (char c : line) {
                char uc = std::toupper(static_cast<unsigned char>(c));
                if (uc == 'A' || uc == 'T' || uc == 'C' || uc == 'G') current_sequence.push_back(uc);
            }
        }
    }
    if (!header.empty()) sequences[header] = current_sequence;
    return sequences;
}

static inline std::string sanitize_for_filename(std::string s) {
    for (char& c : s) {
        if (c==' '||c==','||c=='/'||c=='\\'||c==':'||c=='\t') c = '_';
    }
    return s;
}

static bool parse_args(int argc, char* argv[], Args& a) {
    if (argc == 1) { usage(); return false; }
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto need = [&](bool ok, const std::string& msg){
            if (!ok) die(msg);
        };
        auto grab = [&](std::string& dst){
            need(i+1 < argc, "missing value for " + k);
            dst = argv[++i];
        };
        auto grabi = [&](int& dst){
            need(i+1 < argc, "missing value for " + k);
            dst = std::stoi(argv[++i]);
        };
        auto grabd = [&](double& dst){
            need(i+1 < argc, "missing value for " + k);
            dst = std::stod(argv[++i]);
        };

        if      (k == "--source") grab(a.source);
        else if (k == "--target") grab(a.target);
        else if (k == "--W") grabi(a.W);
        else if (k == "--pcr") grabd(a.cost_pcr);
        else if (k == "--join") grabd(a.cost_join);
        else if (k == "--synth") grabd(a.cost_synth);
        else if (k == "--summary-out") grab(a.summary_out);
        else if (k == "--moves-out") grab(a.moves_out);
        else if (k == "--progress") grabi(a.progress);
        else if (k == "--keep-temp") a.keep_temp = true;
        else if (k == "--help" || k == "-h") { usage(); std::exit(0); }
        else die("unknown argument: " + k);
    }
    if (a.source.empty() || a.target.empty()) die("--source and --target are required");
    if (a.W <= 0) die("--W must be > 0");
    if (a.cost_pcr < 0.0) die("--pcr must be >= 0");
    if (a.cost_join < 0.0) die("--join must be >= 0");
    if (a.cost_synth < 0.0) die("--synth must be >= 0");
    if (a.summary_out.empty()) die("--summary-out is required");
    return true;
}

static void build_index_in_memory(const std::string& source_plain_path, fm_index_t& index) {
    std::cerr << "📚 Creating FM-index from source: " << source_plain_path << "\n";
    // Build in a temp cache directory; ensure cleanup
    fs::path tmpdir = fs::temp_directory_path();
    std::string base = "sdsl_cache_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    sdsl::cache_config config(true, tmpdir.string(), base);
    // 1 = plain text input
    sdsl::construct(index, source_plain_path, config, 1);
    // Attempt cleanup (construct uses cache; the flag true cleans on destructor as well)
    sdsl::util::delete_all_files(config.file_map);
    std::cerr << "✅ Index built (in-memory)\n";
}

struct Move {
    long long start;  // inclusive
    long long end;    // inclusive
    int length;
    char op;          // 'R' reuse (PCR) or 'S' synth
    double acq_cost;
    bool join_applied;
};

struct PlanResult {
    double total_cost = 0.0;
    long long genome_len = 0;
    std::vector<Move> moves;
};

static inline bool exists_in_source(const fm_index_t& index, const char* b, const char* e) {
    // Avoid constructing std::string for speed; sdsl::count accepts iterator pairs.
    return sdsl::count(index, b, e) > 0;
}

static PlanResult solve_dp_with_plan(const std::string& seq,
                                     int W,
                                     const fm_index_t& index,
                                     double cost_pcr,
                                     double cost_join,
                                     double cost_synth_per_base,
                                     int progress_step_percent)
{
    const long long N = static_cast<long long>(seq.size());
    PlanResult res;
    res.genome_len = N;
    if (N == 0) return res;

    const double INF = 1e300;
    std::vector<double> DP(N + 1, INF);
    std::vector<int>    best_w(N + 1, -1);
    std::vector<char>   best_op(N + 1, 'U');
    DP[0] = 0.0;

    auto t0 = std::chrono::steady_clock::now();
    int last_report = -progress_step_percent;

    for (long long i = 1; i <= N; ++i) {
        double best_cost = INF;
        int    best_w_here = -1;
        char   best_op_here = 'U';

        int wmax = std::min<long long>(W, i);
        for (int w = 1; w <= wmax; ++w) {
            long long j = i - w;
            const char* b = &seq[j];
            const char* e = &seq[i];
            bool reuse_ok = exists_in_source(index, b, e);
            double acq_cost = reuse_ok ? cost_pcr : (w * cost_synth_per_base);
            double join_cost = (j > 0) ? cost_join : 0.0;
            double candidate = DP[j] + acq_cost + join_cost;
            if (candidate < best_cost) {
                best_cost = candidate;
                best_w_here = w;
                best_op_here = reuse_ok ? 'R' : 'S';
            }
        }
        DP[i] = best_cost;
        best_w[i] = best_w_here;
        best_op[i] = best_op_here;

        if (progress_step_percent > 0) {
            int pct = static_cast<int>((100.0 * i) / N);
            if (pct >= last_report + progress_step_percent || i == N) {
                last_report = pct;
                auto now = std::chrono::steady_clock::now();
                double sec = std::chrono::duration<double>(now - t0).count();
                std::cerr << "DP progress: " << pct << "% (" << i << "/" << N << ") elapsed " << sec << "s\n";
            }
        }
    }

    res.total_cost = DP[N];

    // Backtrack
    long long i = N;
    while (i > 0) {
        int w = best_w[i];
        char op = best_op[i];
        if (w <= 0 || (op != 'R' && op != 'S')) {
            std::cerr << "Backtrack error at i=" << i << " (w="<<w<<", op="<<op<<")\n";
            break;
        }
        long long j = i - w;
        double acq_cost = (op == 'R') ? cost_pcr : (w * cost_synth_per_base);
        bool join_applied = (j > 0);
        res.moves.push_back(Move{j, i - 1, w, op, acq_cost, join_applied});
        i = j;
    }
    std::reverse(res.moves.begin(), res.moves.end());
    return res;
}

static void write_moves_tsv(std::ofstream& out,
                            const std::string& chrom_name,
                            double cost_join,
                            const PlanResult& plan)
{
    for (size_t k = 0; k < plan.moves.size(); ++k) {
        const auto& m = plan.moves[k];
        out << chrom_name << '\t'
            << (k+1) << '\t'
            << m.start << '\t'
            << m.end << '\t'
            << m.length << '\t'
            << (m.op == 'R' ? "REUSE" : "SYNTH") << '\t'
            << m.acq_cost << '\t'
            << (m.join_applied ? 1 : 0) << '\t'
            << (m.join_applied ? cost_join : 0.0) << '\n';
    }
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    Args a;
    if (!parse_args(argc, argv, a)) return 2;

    std::vector<std::string> temp_files; // for cleanup

    // Inflate gz if needed
    std::cerr << "🧩 Phase 1/3: Preparing FASTAs...\n";
    std::string source_plain = ensure_plain_fasta(a.source, temp_files);
    std::string target_plain = ensure_plain_fasta(a.target, temp_files);
    std::cerr << "✅ FASTAs ready\n\n";

    // Build index in memory
    std::cerr << "📚 Phase 2/3: FM-index build...\n";
    fm_index_t index;
    build_index_in_memory(source_plain, index);
    std::cerr << "\n";

    // Read target and run DP per chromosome
    std::cerr << "🧮 Phase 3/3: DP planning...\n";
    auto chroms = read_fasta_and_clean(target_plain);
    if (chroms.empty()) die("no sequences read from target FASTA: " + target_plain);

    std::ofstream summary(a.summary_out);
    if (!summary.is_open()) die("cannot open summary-out: " + a.summary_out);
    summary << "target_file\tchromosome\tlength\ttotal_cost\tblocks\t"
               "replication_blocks\tsynthesis_blocks\tjoins\t"
               "avg_replication_len\tavg_synthesis_len\t"
               "total_replication_cost\ttotal_synthesis_cost\ttotal_join_cost\n";

    std::ofstream moves;
    bool emit_moves = !a.moves_out.empty();
    if (emit_moves) {
        moves.open(a.moves_out);
        if (!moves.is_open()) die("cannot open moves-out: " + a.moves_out);
        moves << "chromosome\tblock_id\tstart\tend\tlength\top\tacq_cost\tjoin_applied\tjoin_cost\n";
    }

    // Totals across chromosomes
    long long total_len = 0, total_blocks = 0, total_rep_blocks = 0, total_syn_blocks = 0, total_joins = 0;
    long long total_rep_bases = 0, total_syn_bases = 0;
    double total_cost_overall = 0.0, total_rep_cost = 0.0, total_syn_cost = 0.0, total_join_cost = 0.0;

    const std::string target_file_name = fs::path(a.target).filename().string();

    for (const auto& kv : chroms) {
        std::string chrom = sanitize_for_filename(kv.first);
        const std::string& seq = kv.second;
        if (seq.empty()) continue;

        std::cerr << "---- Planning chromosome: " << chrom << " (len=" << seq.size() << ") ----\n";
        auto plan = solve_dp_with_plan(seq, a.W, index, a.cost_pcr, a.cost_join, a.cost_synth, a.progress);

        long long blocks = static_cast<long long>(plan.moves.size());
        long long rep_blocks = 0, syn_blocks = 0, joins = 0;
        long long rep_bases = 0, syn_bases = 0;
        double rep_cost = 0.0, syn_cost = 0.0;

        for (const auto& m : plan.moves) {
            if (m.op == 'R') { rep_blocks++; rep_bases += m.length; rep_cost += m.acq_cost; }
            else             { syn_blocks++; syn_bases += m.length; syn_cost += m.acq_cost; }
            if (m.join_applied) joins++;
        }
        double avg_rep_len = rep_blocks ? (double)rep_bases / rep_blocks : 0.0;
        double avg_syn_len = syn_blocks ? (double)syn_bases / syn_blocks : 0.0;
        double join_cost_total = joins * a.cost_join;

        // Per-chromosome summary
        summary << target_file_name << '\t'
                << chrom << '\t'
                << plan.genome_len << '\t'
                << plan.total_cost << '\t'
                << blocks << '\t'
                << rep_blocks << '\t'
                << syn_blocks << '\t'
                << joins << '\t'
                << avg_rep_len << '\t'
                << avg_syn_len << '\t'
                << rep_cost << '\t'
                << syn_cost << '\t'
                << join_cost_total << '\n';

        // Moves (optional)
        if (emit_moves) write_moves_tsv(moves, chrom, a.cost_join, plan);

        // Accumulate totals
        total_len += plan.genome_len;
        total_blocks += blocks;
        total_rep_blocks += rep_blocks;
        total_syn_blocks += syn_blocks;
        total_joins += joins;
        total_rep_bases += rep_bases;
        total_syn_bases += syn_bases;
        total_rep_cost += rep_cost;
        total_syn_cost += syn_cost;
        total_join_cost += join_cost_total;
        total_cost_overall += plan.total_cost;
    }

    // TOTAL row
    double total_avg_rep_len = total_rep_blocks ? (double)total_rep_bases / total_rep_blocks : 0.0;
    double total_avg_syn_len = total_syn_blocks ? (double)total_syn_bases / total_syn_blocks : 0.0;
    summary << "TOTAL\tALL\t" << total_len << '\t' << total_cost_overall << '\t'
            << total_blocks << '\t'
            << total_rep_blocks << '\t'
            << total_syn_blocks << '\t'
            << total_joins << '\t'
            << total_avg_rep_len << '\t'
            << total_avg_syn_len << '\t'
            << total_rep_cost << '\t'
            << total_syn_cost << '\t'
            << total_join_cost << '\n';

    summary.close();
    if (moves.is_open()) moves.close();

    std::cerr << "✅ DP planning done.\n";
    std::cerr << "📊 Summary TSV: " << a.summary_out << "\n";
    if (emit_moves) std::cerr << "📄 Moves TSV:   " << a.moves_out << "\n";

    // Cleanup temp files
    if (!a.keep_temp) {
        for (auto& p : temp_files) {
            std::error_code ec;
            fs::remove(p, ec);
        }
    } else {
        std::cerr << "ℹ️  Keeping temp files (--keep-temp)\n";
        for (auto& p : temp_files) std::cerr << "  " << p << "\n";
    }

    return 0;
}

