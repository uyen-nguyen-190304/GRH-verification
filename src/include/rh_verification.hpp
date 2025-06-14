#ifndef RH_VERIFICATION
#define RH_VERIFICATION

#include <vector>
#include <cstdint>

namespace grh {

// Compute iota(eta)
double iota(double eta);

// Compute C(Z)
double C_Z(const std::vector<std::pair<double, double>>& intervals);

// Compute logarithmic derivative
double log_derivative(const std::vector<int8_t>& chi_arr, 
                      const std::vector<double>& lambda_arr, 
                      int K);

// RH Inequality test
bool rh_verify(long long d, int K, double eta,
               const std::vector<std::pair<double, double>>& intervals,
               const std::vector<int8_t>& chi_arr,
               const std::vector<double>& lambda_arr,
               /* out */ int& N_used,
               /* out */ double& lhs,
               /* out */ double& rhs);

} // namespace grh

#endif /* RH_VERIFICATION */
