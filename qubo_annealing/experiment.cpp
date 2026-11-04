#include "qubo.h"

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "annealing.h"
#include "experiment.h"
// Runtime settings and optional diagnostic requests parsed from the command line.
struct Options {
    std::string input, output, schedule = "geometric", initial, energy_bits, delta_bits;
    std::string selection = "sequential";
    std::uint64_t sweeps_per_temperature = 1000;
    int delta_index = -1;
    std::uint64_t sweeps = 1000, trials = 1, seed = 1, checkpoint = 100;
    double start = 10, end = 0.1, target = 0;
    bool has_target = false;
};

// Validate a CLI bitstring and convert it to the solver's byte-per-bit state.
static std::vector<unsigned char> parse_bits(const std::string& bits, int n) {
    if (bits.size() != static_cast<std::size_t>(n)) throw std::runtime_error("bitstring has wrong length");
    std::vector<unsigned char> state(n);
    for (int i = 0; i < n; ++i) {
        if (bits[i] != '0' && bits[i] != '1') throw std::runtime_error("bitstring must contain only 0 and 1");
        state[i] = static_cast<unsigned char>(bits[i] - '0');
    }
    return state;
}

// Escape a value before writing it as a JSON string in a trial record.
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

// Convert the solver's byte-per-bit state to the saved 0/1 bitstring.
static std::string bit_string(const std::vector<unsigned char>& state) {
    std::string bits;
    bits.reserve(state.size());
    for (auto bit : state) bits += char('0' + bit);
    return bits;
}

// Use the input filename without its extension as the benchmark identifier.
static std::string benchmark_id(const std::string& path) {
    const auto slash = path.find_last_of("/\\");
    const std::string name = path.substr(slash == std::string::npos ? 0 : slash + 1);
    const auto dot = name.find_last_of('.');
    return name.substr(0, dot);
}

// Print the supported command-line modes and experiment settings.
static void usage() {
    std::cout << "Usage: qubo_anneal INSTANCE.qubo [options]\n"
              << "  --sweeps N --trials N --seed N --initial BITS\n"
              << "  --schedule linear|geometric|fixed|paper25 --temp-start T --temp-end T\n"
              << "  --selection sequential|random --sweeps-per-temperature N (paper25)\n"
              << "  --checkpoint N (0 disables intermediate checkpoints)\n"
              << "  --target ENERGY --output results.jsonl\n"
              << "  --energy BITS | --delta BITS INDEX (diagnostic modes)\n";
}

// Parse a nonnegative integer option, rejecting partial or negative input.
static std::uint64_t integer(const std::string& value) {
    if (value.empty() || value[0] == '-') throw std::runtime_error("expected nonnegative integer");
    std::size_t used = 0;
    auto result = std::stoull(value, &used);
    if (used != value.size()) throw std::runtime_error("invalid integer: " + value);
    return result;
}

// Parse a finite real-valued option, rejecting trailing characters.
static double real(const std::string& value) {
    std::size_t used = 0;
    double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result)) throw std::runtime_error("invalid number: " + value);
    return result;
}

// Turn CLI arguments into checked settings for diagnostics or annealing trials.
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
        else if (key == "--selection") o.selection = value;
        else if (key == "--sweeps-per-temperature") o.sweeps_per_temperature = integer(value);
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
        (o.schedule != "linear" && o.schedule != "geometric" && o.schedule != "fixed" && o.schedule != "paper25") ||
        (o.selection != "sequential" && o.selection != "random"))
        throw std::runtime_error("invalid trials, temperature, or schedule");
    if (o.schedule == "paper25") {
        if (!o.sweeps_per_temperature || o.sweeps_per_temperature > std::numeric_limits<std::uint64_t>::max()/25)
            throw std::runtime_error("invalid sweeps per temperature");
        o.sweeps = 25 * o.sweeps_per_temperature;
        o.start = std::pow(1.05, -1.0); o.end = std::pow(1.05, -25.0);
        o.checkpoint = o.sweeps_per_temperature;
    }
    if (!o.energy_bits.empty() && !o.delta_bits.empty()) throw std::runtime_error("choose one diagnostic mode");
    return o;
}

// Write the common sweep/update/time shape used by first_best and first_hit.
static void write_event(std::ostream& out, const Event& event) {
    out << "{\"sweep\":" << event.sweep << ",\"attempted_updates\":" << event.attempted_updates
        << ",\"elapsed_seconds\":" << event.elapsed_seconds << '}';
}

// Write one JSONL record with enough settings to reproduce this seeded trial.
static void write_trial(std::ostream& out, const Qubo& q, const Options& o,
                        std::uint64_t trial, const AnnealingResult& result) {
    out << "{\"trial\":" << trial << ",\"seed\":" << result.seed
        << ",\"input\":" << json_string(o.input)
        << ",\"benchmark\":" << json_string(benchmark_id(o.input))
        << ",\"variables\":" << q.n
        << ",\"schedule\":" << json_string(o.schedule)
        << ",\"selection\":" << json_string(o.selection)
        << ",\"sweeps_per_temperature\":" << (o.schedule == "paper25" ? o.sweeps_per_temperature : 0)
        << ",\"temperature_start\":" << o.start << ",\"temperature_end\":" << o.end
        << ",\"sweeps\":" << o.sweeps << ",\"checkpoint\":" << o.checkpoint
        << ",\"initial_bits\":";
    if (o.initial.empty()) out << "null";
    else out << json_string(o.initial);
    out << ",\"target_energy\":";
    if (o.has_target) out << o.target;
    else out << "null";
    out << ",\"attempted_updates\":" << result.attempted_updates
        << ",\"accepted_updates\":" << result.accepted_updates
        << ",\"best_energy\":" << result.best_energy
        << ",\"best_bits\":" << json_string(bit_string(result.best_state))
        << ",\"final_energy\":" << result.final_energy
        << ",\"final_bits\":" << json_string(bit_string(result.final_state))
        << ",\"first_best\":";
    write_event(out, result.first_best);
    out << ",\"progress\":[";
    for (std::size_t i = 0; i < result.progress.size(); ++i) {
        if (i) out << ',';
        const auto& point = result.progress[i];
        out << "{\"sweep\":" << point.sweep << ",\"attempted_updates\":" << point.attempted_updates
            << ",\"best_energy\":" << point.best_energy
            << ",\"elapsed_seconds\":" << point.elapsed_seconds;
        if (!point.bits.empty())
            out << ",\"current_energy\":" << point.current_energy
                << ",\"bits\":" << json_string(point.bits)
                << ",\"beta\":" << std::pow(1.05, double(point.sweep / o.sweeps_per_temperature));
        out << '}';
    }
    out << "],\"target_reached\":" << (result.target_reached ? "true" : "false")
        << ",\"first_hit\":";
    if (result.target_reached) write_event(out, result.first_hit);
    else out << "null";
    out << ",\"runtime_seconds\":"
        << result.runtime_seconds << "}\n";
}

// Load the fixed QUBO once, handle optional diagnostics, then run and log trials.
int run_experiment(int argc, char** argv) {
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
            auto field = local_fields(q, x);
            out << "{\"delta\":" << flip_delta(x, field, o.delta_index) << "}\n";
            return 0;
        }
        if (!o.initial.empty()) parse_bits(o.initial, q.n);
        AnnealingParameters parameters;
        parameters.schedule = o.schedule;
        parameters.selection = o.selection;
        parameters.sweeps = o.sweeps;
        parameters.sweeps_per_temperature = o.sweeps_per_temperature;
        parameters.checkpoint = o.checkpoint;
        parameters.start = o.start;
        parameters.end = o.end;
        parameters.target = o.target;
        parameters.has_target = o.has_target;
        if (!o.initial.empty()) parameters.initial_state = parse_bits(o.initial, q.n);
        for (std::uint64_t trial = 0; trial < o.trials; ++trial) {
            const auto result = simulated_annealing(q, parameters, o.seed + trial);
            write_trial(out, q, o, trial, result);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "qubo_anneal: " << error.what() << '\n';
        return 2;
    }
}
