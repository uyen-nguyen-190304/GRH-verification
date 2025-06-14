#include "include/rh_verification.hpp"
#include <cmath>
#include <vector>
#include <utility>
#include <cstddef>

namespace grh {

/* 
 @brief
 Compute the iota(eta) constant used in RH estimate. Return the minimum of 
 two expressions in eta

 @param: eta - width of window where one would like to verify the RH
 @return: minimum value of two expressions
 */
double iota(double eta) 
{
    double term1 = 1.0 / (1.0 + eta * eta) + 2.0 / (4.0 + eta * eta);
    double term2 = 12.0 / (9.0 + 4.0 * eta * eta);
    return (term1 < term2) ? term1 : term2;
}

/*
 @brief
 Compute the C(Z) sum over zero intervals 
 ! Utility only, not used in rh_verify below 
 
 @param: intervals - list of intervals around zeros
 @return: sum of contributions from each zero intervals, depends on type
 */
double C_Z(const std::vector<std::pair<double, double>>& intervals) 
{
    double sum = 0.0;

    // Loop over all disjoint subintervals
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        const double gamma_minus = intervals[i].first;
        const double gamma_plus  = intervals[i].second;

        // Separate by type
        if (std::fabs(gamma_minus + gamma_plus) < 1e-8) 
        {
            // Type 2: symmetric [-gamma0, gamma0]
            const double gamma0 = std::fabs(gamma_plus);
            sum += 6.0 / (9.0 + 4.0 * gamma0 * gamma0);
        } 
        else 
        {
            // Type 1: general [gamma_minus, gamma_plus]
            sum += 12.0 / (9.0 + 4.0 * gamma_plus * gamma_plus);
        }
    }
    return sum;
}

/*
 @brief
 Compute the logarithmic derivative of L(s) at s = 2
 ? Well, I might need to include delta to such parameters

 @param: chi_arr - array of χ_d(k)
         lambda_arr - array of Λ(k)
         K - upper bound for k in each of the array above
 @return: L'(2, χ_d)/L(2, χ_d) value
 */
double log_derivative(const std::vector<int8_t>& chi_arr, 
                      const std::vector<double>& lambda_arr, 
                      int K) 
{
    double sum = 0.0;

    for (int k = K; k <= K; ++k) 
    {
        const int8_t chi_k = chi_arr[static_cast<std::size_t>(k)];
        if (chi_k == 0) continue;   // skip k where χ_d(k) = 0
        
        const double lambda_k = lambda_arr[static_cast<std::size_t>(k)];
        const double k_sq     = static_cast<double>(k) * static_cast<double>(k);

        sum -= static_cast<double>(chi_k) * lambda_k / k_sq;
    }
    return sum;
}

/*
 @brief
 Return the single-zero contribution to the LHS of the RH-corollary inequality

 @param: gamma_minus, gamma_plus : endpoints of the enclosing interval

 @return: 6 / (9 + 4 * γ^2) if the interval is symmetric 
          12 / (9 + 4 * γ+^2) otherwise
 */
double zero_contribution(double gamma_minus, double gamma_plus)
{
    if (std::fabs(gamma_minus + gamma_plus) < 1e-8)   
    {
        // Type 2: Symmetric [−γ, γ]
        const double g = std::fabs(gamma_plus);
        return 6.0 / (9.0 + 4.0 * g * g);
    }
    // Type 1: Asymmetric [γ⁻, γ⁺]
    return 12.0 / (9.0 + 4.0 * gamma_plus * gamma_plus);
}

/*
 @brief
 Main RH inequality verifier
 ! One-shot functionality when all zero intervals are provided 
 
 @param: d - fundamental discriminant
         K - Upper bound for logarithmic derivative computation
         eta - height of interest to verify the RH
         intervals - list of disjoint intervals around zeros
         chi_arr - Kronecker χ_d(k) array
         lambda_arr - von Mangoldt Λ(k) array

         N_used - output, how many zeros were used
         lhs - output, left-hand side of the inequality
         rhs - output, right-hand side of the inequality

 @return: true if LHS > RHS, else false
 */
bool rh_verify(long long d, int K, double eta,
               const std::vector<std::pair<double, double>>& intervals,
               const std::vector<int8_t>& chi_arr,
               const std::vector<double>& lambda_arr,
               /* out */ int& N_used,
               /* out */ double& lhs,
               /* out */ double& rhs)
{
    /* Compute the RHS, depends on sign of the fundamental discriminant d */
    if (d < 0)
    {
        rhs = 0.5 * std::log(std::fabs(static_cast<double>(d)) * std::exp(2) / (4.0 * M_PI * std::exp(EULER_CONSTANT))) \
                + log_derivative(chi_arr, lambda_arr, K);
    } 
    else 
    { 
        rhs = 0.5 * std::log(static_cast<double>(d) / (M_PI * std::exp(EULER_CONSTANT))) \
                + log_derivative(chi_arr, lambda_arr, K);
    }

    // Initialize the LHS
    lhs = 2.0 * iota(eta);
    N_used = 0;     // Counter for number of zeros (intervals) used

    // Accumulate zero contributions until LHS > RHS
    for (std::size_t i = 0; i < intervals.size(); ++i)
    { 
        ++N_used;

        const double gamma_minus = intervals[i].first;
        const double gamma_plus  = intervals[i].second;
        
        if (std::fabs(gamma_minus + gamma_plus) < 1e-8) 
        {
            // Type 2: Symmetric
            const double gamma0 = std::fabs(gamma_plus);
            lhs += 6.0 / (9.0 + 4.0 * gamma0 * gamma0);
        } 
        else 
        {
            // Type 1: Asymmetric
            lhs += 12.0 / (9.0 + 4.0 * gamma_plus * gamma_plus);
        }    

        // If LHS > RHS, the inequality is satisfied and RH is verified
        if (lhs > rhs)
            return true;            
    }
    // Not enough zeros to satisfy the inequality
    return false;                   // not enough zeros        
}

} // namespace grh
