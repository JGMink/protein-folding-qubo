#ifndef QUBO_MODEL_H
#define QUBO_MODEL_H

#include <string>
#include <utility>
#include <vector>

// The problem: E(x) = offset + sum_i diagonal[i]*x_i + sum_{i<j} Q_ij*x_i*x_j, with each x_i in {0,1}.
// Read-only after loading, so replicas (or MPI ranks, each loading its own copy) can share it.
// (`std::vector` is a heap array that knows its length; `neighbors` is a ragged array of (j, Q_ij) pairs.)
struct Qubo {
    int n = 0;                       // number of binary variables (bits)
    double offset = 0;               // constant added to every energy
    std::vector<double> diagonal;    // Q_ii, the linear term of each bit
    std::vector<std::vector<std::pair<int, double>>> neighbors;  // sparse Q_ij; each pair is in both bits' lists
};

// Read a QUBO1 text file into a Qubo. Throws std::runtime_error on any malformed input.
Qubo load_qubo(const std::string& path);

// Full energy computed from scratch; used at the start and for final reporting, not per flip.
// `bits` holds one byte (0 or 1) per variable; `const&` is a read-only reference, like a const pointer.
double energy(const Qubo& q, const std::vector<unsigned char>& bits);

// field[i] = diagonal[i] + sum of Q_ij over neighbors j that are currently 1: the energy cost of switching bit i on.
std::vector<double> local_fields(const Qubo& q, const std::vector<unsigned char>& bits);

// Energy change if bit `index` flips, in O(1) from the cached fields (0->1 adds the field, 1->0 subtracts it).
double flip_delta(const std::vector<unsigned char>& bits, const std::vector<double>& fields, int index);

// After bit `index` flips (change = +1 for 0->1, -1 for 1->0), shift each neighbor's field by weight * change.
// The caller still toggles the bit. `std::vector<double>&` is a non-const reference the function writes through.
void update_fields_after_flip(const Qubo& q, std::vector<double>& fields, int index, int change);

// ---- Optional fold-validity check (protein-folding QUBOs only) ----
// Needs the lattice from the file's meta line. Variable index = residue * sites + site; a bit is 1 when that residue sits there.
struct FoldLayout {
    int residues = 0;                       // chain length
    int sites = 0;                          // lattice sites; residues * sites must equal Qubo::n
    std::vector<std::vector<int>> next_to;  // next_to[s] = the sites that share a lattice edge with s
};

// Read the lattice from the meta line of a QUBO1 file with n variables (needs "sequence" and "lattice_edges",
// as written compactly by qubo_bench). Throws std::runtime_error if either is missing or does not fit n.
FoldLayout load_fold_layout(const std::string& path, int n);

// True when `bits` is a legal fold (same rules as qubo_bench.model.validate): one site per residue, no shared
// sites, consecutive residues on adjacent sites. Costs O(n), so use it for rare checks, not on every flip.
bool is_valid_fold(const FoldLayout& layout, const std::vector<unsigned char>& bits);

#endif
