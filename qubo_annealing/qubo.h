#ifndef QUBO_MODEL_H
#define QUBO_MODEL_H

#include <string>
#include <utility>
#include <vector>

struct Qubo {
    int n = 0;
    double offset = 0;
    std::vector<double> diagonal;
    std::vector<std::vector<std::pair<int, double>>> neighbors;
};

Qubo load_qubo(const std::string& path);
double energy(const Qubo& q, const std::vector<unsigned char>& bits);
std::vector<double> local_fields(const Qubo& q, const std::vector<unsigned char>& bits);
double flip_delta(const std::vector<unsigned char>& bits, const std::vector<double>& fields, int index);
void update_fields_after_flip(const Qubo& q, std::vector<double>& fields, int index, int change);

#endif
