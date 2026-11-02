// Standalone serial binary-QUBO simulated annealing. No third-party code.
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

struct Qubo {
    int n = 0;
    double offset = 0;
    std::vector<double> diagonal;
    std::vector<std::vector<std::pair<int, double>>> neighbors;
};

static bool read_line(std::istream& in, std::string& line) {
    if (!std::getline(in, line)) return false;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return true;
}

static std::string json_string(const std::string& value) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c == '\n') out << "\\n";
        else if (c == '\r') out << "\\r";
        else if (c == '\t') out << "\\t";
        else if (c < 32) {
            const char* hex = "0123456789abcdef";
            out << "\\u00" << hex[c >> 4] << hex[c & 15];
        } else out << c;
    }
    return out.str() + '"';
}

static Qubo load_qubo(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open QUBO input: " + path);
    std::string line, word;
    if (!read_line(in, line) || line != "QUBO1") throw std::runtime_error("expected QUBO1 header");
    Qubo q;
    if (!read_line(in, line)) throw std::runtime_error("missing variable count");
    {
        std::istringstream row(line);
        if (!(row >> word >> q.n) || word != "n" || q.n < 1 || (row >> word))
            throw std::runtime_error("invalid variable count");
    }
    if (!read_line(in, line)) throw std::runtime_error("missing offset");
    {
        std::istringstream row(line);
        if (!(row >> word >> q.offset) || word != "offset" || !std::isfinite(q.offset) || (row >> word))
            throw std::runtime_error("invalid offset");
    }
    if (!read_line(in, line) || line.rfind("meta ", 0) != 0)
        throw std::runtime_error("missing metadata record");
    q.diagonal.assign(q.n, 0);
    q.neighbors.resize(q.n);
    std::unordered_set<std::uint64_t> seen;
    while (read_line(in, line)) {
        std::istringstream row(line);
        int i, j;
        double coefficient;
        if (!(row >> word >> i >> j >> coefficient) || word != "q" || (row >> word)
            || i < 0 || i > j || j >= q.n || !std::isfinite(coefficient))
            throw std::runtime_error("invalid QUBO coefficient record");
        std::uint64_t key = std::uint64_t(i) * std::uint64_t(q.n) + std::uint64_t(j);
        if (!seen.insert(key).second) throw std::runtime_error("duplicate QUBO coefficient");
        if (i == j) q.diagonal[i] = coefficient;
        else {
            q.neighbors[i].push_back({j, coefficient});
            q.neighbors[j].push_back({i, coefficient});
        }
    }
    return q;
}

static std::vector<unsigned char> parse_bits(const std::string& bits, int n) {
    if (bits.size() != static_cast<std::size_t>(n)) throw std::runtime_error("bitstring has wrong length");
    std::vector<unsigned char> state(n);
    for (int i = 0; i < n; ++i) {
        if (bits[i] != '0' && bits[i] != '1') throw std::runtime_error("bitstring must contain only 0 and 1");
        state[i] = static_cast<unsigned char>(bits[i] - '0');
    }
    return state;
}

static std::string bit_string(const std::vector<unsigned char>& state) {
    std::string bits;
    bits.reserve(state.size());
    for (auto bit : state) bits += char('0' + bit);
    return bits;
}

static double energy(const Qubo& q, const std::vector<unsigned char>& x) {
    double result = q.offset;
    for (int i = 0; i < q.n; ++i) {
        if (!x[i]) continue;
        result += q.diagonal[i];
        for (const auto& [j, weight] : q.neighbors[i])
            if (j > i && x[j]) result += weight;
    }
    return result;
}

static std::vector<double> fields(const Qubo& q, const std::vector<unsigned char>& x) {
    std::vector<double> result = q.diagonal;
    for (int i = 0; i < q.n; ++i)
        for (const auto& [j, weight] : q.neighbors[i])
            if (x[j]) result[i] += weight;
    return result;
}

struct Options {
    std::string input, output, schedule = "geometric", initial, energy_bits, delta_bits;
    int delta_index = -1;
    std::uint64_t sweeps = 1000, trials = 1, seed = 1, checkpoint = 100;
    double start = 10, end = 0.1, target = 0;
    bool has_target = false;
};

static void usage() {
    std::cout << "Usage: qubo_sa INSTANCE.qubo [options]\n"
              << "  --sweeps N --trials N --seed N --initial BITS\n"
              << "  --schedule linear|geometric|fixed --temp-start T --temp-end T\n"
              << "  --checkpoint N (0 disables intermediate checkpoints)\n"
              << "  --target ENERGY --output results.jsonl\n"
              << "  --energy BITS | --delta BITS INDEX (diagnostic modes)\n";
}

static std::uint64_t integer(const std::string& value) {
    if (value.empty() || value[0] == '-') throw std::runtime_error("expected nonnegative integer");
    std::size_t used = 0;
    auto result = std::stoull(value, &used);
    if (used != value.size()) throw std::runtime_error("invalid integer: " + value);
    return result;
}

static double real(const std::string& value) {
    std::size_t used = 0;
    double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result)) throw std::runtime_error("invalid number: " + value);
    return result;
}

static Options options(int argc, char** argv) {
    Options o;
    for (int a = 1; a < argc; ++a) {
        std::string key = argv[a];
        if (key == "--help") { usage(); std::exit(0); }
        if (key.empty() || key[0] != '-') {
            if (!o.input.empty()) throw std::runtime_error("only one input file is allowed");
            o.input = key;
            continue;
        }
        if (++a >= argc) throw std::runtime_error("missing value for " + key);
        std::string value = argv[a];
        if (key == "--sweeps") o.sweeps = integer(value);
        else if (key == "--trials") o.trials = integer(value);
        else if (key == "--seed") o.seed = integer(value);
        else if (key == "--checkpoint") o.checkpoint = integer(value);
        else if (key == "--temp-start") o.start = real(value);
        else if (key == "--temp-end") o.end = real(value);
        else if (key == "--schedule") o.schedule = value;
        else if (key == "--initial") o.initial = value;
        else if (key == "--output") o.output = value;
        else if (key == "--target") { o.target = real(value); o.has_target = true; }
        else if (key == "--energy") o.energy_bits = value;
        else if (key == "--delta") {
            o.delta_bits = value;
            if (++a >= argc) throw std::runtime_error("--delta requires a bit index");
            auto index = integer(argv[a]);
            if (index > std::uint64_t(std::numeric_limits<int>::max())) throw std::runtime_error("index too large");
            o.delta_index = static_cast<int>(index);
        } else throw std::runtime_error("unknown option: " + key);
    }
    if (o.input.empty()) throw std::runtime_error("input file required");
    if (o.trials == 0 || o.start <= 0 || o.end <= 0 ||
        (o.schedule != "linear" && o.schedule != "geometric" && o.schedule != "fixed"))
        throw std::runtime_error("invalid trials, temperature, or schedule");
    if (!o.energy_bits.empty() && !o.delta_bits.empty()) throw std::runtime_error("choose one diagnostic mode");
    return o;
}

static double temperature(const Options& o, std::uint64_t sweep) {
    if (o.schedule == "fixed") return o.start;
    double fraction = o.sweeps <= 1 ? 1.0 : double(sweep) / double(o.sweeps - 1);
    if (o.schedule == "linear") return o.start + (o.end - o.start) * fraction;
    return o.start * std::pow(o.end / o.start, fraction);
}

int main(int argc, char** argv) {
    try {
        const Options o = options(argc, argv);
        const Qubo q = load_qubo(o.input);
        std::ofstream file;
        if (!o.output.empty()) {
            file.open(o.output);
            if (!file) throw std::runtime_error("cannot open output file");
        }
        std::ostream& out = o.output.empty() ? std::cout : file;
        out << std::setprecision(17);
        if (!o.energy_bits.empty()) {
            out << "{\"energy\":" << energy(q, parse_bits(o.energy_bits, q.n)) << "}\n";
            return 0;
        }
        if (!o.delta_bits.empty()) {
            auto x = parse_bits(o.delta_bits, q.n);
            if (o.delta_index >= q.n) throw std::runtime_error("delta index out of range");
            auto field = fields(q, x);
            out << "{\"delta\":" << (1 - 2 * int(x[o.delta_index])) * field[o.delta_index] << "}\n";
            return 0;
        }
        if (!o.initial.empty()) parse_bits(o.initial, q.n);
        for (std::uint64_t trial = 0; trial < o.trials; ++trial) {
            const auto clock_start = std::chrono::steady_clock::now();
            const auto elapsed = [&]() {
                return std::chrono::duration<double>(std::chrono::steady_clock::now() - clock_start).count();
            };
            std::uint64_t trial_seed = o.seed + trial;
            std::mt19937_64 rng(trial_seed);
            std::vector<unsigned char> x(q.n);
            if (!o.initial.empty()) x = parse_bits(o.initial, q.n);
            else for (auto& bit : x) bit = static_cast<unsigned char>(rng() & 1);
            auto field = fields(q, x);
            double current = energy(q, x), best = current;
            auto best_state = x;
            std::uint64_t attempts = 0, accepted = 0;
            bool hit = o.has_target && best <= o.target + 1e-9;
            std::uint64_t hit_sweep = 0, hit_attempt = 0;
            double hit_time = hit ? elapsed() : 0;
            std::ostringstream progress;
            progress << std::setprecision(17);
            progress << "[{\"sweep\":0,\"attempted_updates\":0,\"best_energy\":" << best
                     << ",\"elapsed_seconds\":" << elapsed() << "}";
            for (std::uint64_t sweep = 0; sweep < o.sweeps; ++sweep) {
                double temp = temperature(o, sweep);
                for (int i = 0; i < q.n; ++i) {
                    double delta = (1 - 2 * int(x[i])) * field[i];
                    ++attempts;
                    if (delta > 0 && std::generate_canonical<double, 53>(rng) >= std::exp(-delta / temp))
                        continue;
                    ++accepted;
                    int change = 1 - 2 * int(x[i]);
                    for (const auto& [j, weight] : q.neighbors[i]) field[j] += weight * change;
                    x[i] ^= 1;
                    current += delta;
                    if (current < best - 1e-9) { best = current; best_state = x; }
                    if (!hit && o.has_target && current <= o.target + 1e-9) {
                        hit = true; hit_sweep = sweep + 1; hit_attempt = attempts; hit_time = elapsed();
                    }
                }
                if ((o.checkpoint && (sweep + 1) % o.checkpoint == 0) || sweep + 1 == o.sweeps)
                    progress << ",{\"sweep\":" << sweep + 1 << ",\"attempted_updates\":" << attempts
                             << ",\"best_energy\":" << best << ",\"elapsed_seconds\":" << elapsed() << "}";
            }
            progress << ']';
            // Recompute reported energies, avoiding accumulated floating-point drift.
            double final_energy = energy(q, x);
            double best_energy = energy(q, best_state);
            out << "{\"trial\":" << trial << ",\"seed\":" << trial_seed
                << ",\"input\":" << json_string(o.input) << ",\"variables\":" << q.n
                << ",\"schedule\":" << json_string(o.schedule)
                << ",\"temperature_start\":" << o.start << ",\"temperature_end\":" << o.end
                << ",\"sweeps\":" << o.sweeps << ",\"attempted_updates\":" << attempts
                << ",\"accepted_updates\":" << accepted << ",\"best_energy\":" << best_energy
                << ",\"best_bits\":" << json_string(bit_string(best_state))
                << ",\"final_energy\":" << final_energy
                << ",\"final_bits\":" << json_string(bit_string(x))
                << ",\"progress\":" << progress.str()
                << ",\"first_hit\":";
            if (hit) out << "{\"sweep\":" << hit_sweep << ",\"attempted_updates\":" << hit_attempt
                         << ",\"elapsed_seconds\":" << hit_time << "}";
            else out << "null";
            out << ",\"runtime_seconds\":" << elapsed() << "}\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "qubo_sa: " << error.what() << '\n';
        return 2;
    }
}
