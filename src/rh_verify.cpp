#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

// Euler-Mascheroni constant
// ? How many digits of precision are needed?
const double EULER_CONSTANT = 0.57721566490153286060651209008240243;

// Compute iota(eta)
double iota(double eta) {
    double term1 = 1.0 / (1.0 + eta * eta) + 2.0 / (4.0 + eta * eta);
    double term2 = 12.0 / (9.0 + 4.0 * eta * eta);
    return (term1 < term2) ? term1 : term2;
}

// Load intervals from file
vector<pair<double, double>> load_intervals(const string& filename) {
    vector<pair<double, double>> intervals;

    // Check if the file exists and can be opened
    ifstream infile(filename.c_str());
    if (!infile) {
        cerr << "Error: Could not open intervals file: " << filename << endl;
        exit(1);
    }

    // Read the file line by line
    string line;
    while (getline(infile, line)) {        
        // Parse each line into gamma_minus and gamma_plus
        istringstream iss(line);
        double gamma_minus, gamma_plus;
        if (!(iss >> gamma_minus >> gamma_plus)) {
            continue;  // Skip invalid lines (not having two doubles)
        }
        // Store the interval
        intervals.push_back(make_pair(gamma_minus, gamma_plus));
    }
    return intervals;
}

// Compute C(Z)
double C_Z(const vector<pair<double, double>>& intervals) {
    double sum = 0.0;

    // Loop over all disjoint subintervals 
    for (size_t i = 0; i < intervals.size(); ++i) {
        double gamma_minus = intervals[i].first;
        double gamma_plus  = intervals[i].second;

        // Separate by type
        if (fabs(gamma_minus + gamma_plus) < 1e-8) {    // ? Is this tolerance good enough?
            // Type 2: symmetric [-gamma0, gamma0]
            double gamma0 = fabs(gamma_plus);           // ? I don't like this. I'll change it later :) - Anw, from how we construct the file, how can we get this?
            sum += 6.0 / (9.0 + 4.0 * gamma0 * gamma0);
        } else {
            // Type 1: general [gamma_minus, gamma_plus]
            sum += 12.0 / (9.0 + 4.0 * gamma_plus * gamma_plus);
        }
    }
    return sum;
}

// Load precomputed Kronecker symbols from file
map<int, int> load_kronecker(const string& filename) {
    map<int, int> chi_d_arr;

    // Check if the file exists and can be opened
    ifstream infile(filename.c_str());
    if (!infile) {
        cerr << "Error: Could not open Kronecker file: " << filename << endl;
        exit(1);
    }

    // Read the file line by line
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        int n = 0;
        int chi_d_val = 0;
        
        // Parse each line into n and value
        if (!(iss >> n >> chi_d_val)) {
            // Skip invalid lines (not having two integers)
            continue;
        }
        // Store the Kronecker symbol
        chi_d_arr[n] = chi_d_val;       
    }
    return chi_d_arr;
}

// Load precomputed von Mangoldt function values from file
map<int, double> load_von_mangoldt(const string& filename) {
    map<int, double> lambda_arr;

    // Check if the file exists and can be opened
    ifstream infile(filename.c_str());
    if (!infile) {
        cerr << "Error: Could not open von Mangoldt file:" << filename << endl;
        exit(1);
    }

    // Read the file line by line
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        int n = 0;
        double lambda_val = 0.0;

        // Parse each line into n and value
        if (!(iss >> n >> lambda_val)) {
            // Skip invalid lines (not having two values)
            continue;
        }
        // Store the von Mangoldt function value
        lambda_arr[n] = lambda_val;
    }
    return lambda_arr;
}

double logarithmic_derivative(int delta, int K, const map<int, int>& chi_d_arr, const map<int, double>& lambda_arr) {
    double sum = 0.0;

    // Loop over all k from 1 to K
    for (int k = 1; k <= K; ++k) {
        // Compute the value of Lambda_L
        double lambda_L = lambda_arr.at(k) * chi_d_arr.at(k);

        // Add the contribution to the sum
        sum -= lambda_L / (pow(k, 1.0 - delta));
    }
    return sum;
}


int main(int argc, char* argv[]) {
    // Check number of arguments
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <d> <eta> <K>\n";
        return 1;
    } 

    // Parse command line arguments
    int d = atoi(argv[1]);          // Degree of the Dirichlet character
    double eta = atof(argv[2]);     // Height of interest to check RH
    int K = atoi(argv[3]);          // Number of nontrivial zeros to consider

    // Load the data
    vector<pair<double, double>> intervals = load_intervals("data/intervals.txt");
    map<int, int> chi_d_arr = load_kronecker("data/kronecker.txt");
    map<int, double> lambda_arr = load_von_mangoldt("data/von_mangoldt.txt");

    // RH inequality test
    double lhs = 2.0 * iota(eta) + C_Z(intervals);
    double rhs = 0.5 * log(static_cast<double>(d) / (M_PI * exp(EULER_CONSTANT))) \
                 + logarithmic_derivative(-1, K, chi_d_arr, lambda_arr);

    if (lhs > rhs) {
        cout << "Condition satisfied: RH holds for all the nontrivial zeros of L(s, chi_d) up to height eta.\n";
    } else {
        cout << "Condition not satisfied: cannot conclude RH holds for all the nontrivial zeros of L(s, chi_d) up to height eta.\n";
    }

    return 0;
}