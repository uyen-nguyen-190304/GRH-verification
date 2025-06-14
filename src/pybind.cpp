#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  
#include "include/rh_verification.hpp"

namespace py = pybind11;
using namespace grh;

/* Wraper for Python tuple return */
static py::tuple py_rh_verify(long long d, int K, double eta,
                              const std::vector<std::pair<double, double>>& intervals,
                              const std::vector<int8_t>& chi_arr,
                              const std::vector<double>& lambda_arr)
{
      int N_used = 0;     /* Number of zeros used to verify RH */
      double lhs = 0.0;   /* LHS expression of the inequality */
      double rhs = 0.0;   /* RHS expression of the inequality*/

      bool verified = rh_verify(d, K, eta, intervals, 
                                chi_arr, lambda_arr,
                                N_used, lhs, rhs);

      return py::make_tuple(verified, N_used, lhs, rhs);
}

/* Create Python module */
PYBIND11_MODULE(grhverify, m)
{
      m.doc() = "C++ computation for verifying RH numerically";

      // Euler-Mascheroni constant
      m.attr("EULER_CONSTANT") = py::float_(EULER_CONSTANT);

      // iota(eta)
      m.def("iota", &iota,
            py::arg("eta"),
            "Return iota(eta) = min{ 1/(1 + eta^2) + 2 / (4 + eta^2), 12 / (9 + 4 * eta^2) }");

      // C(Z)
      m.def("C_Z", &C_Z,
            py::arg("intervals"),
            "Compute C(Z) = sum_{gamma-, gamma+} (12 / (9 + 4 * gamma+^2)) + sum_{-gamma0, gamma0} (6 / (9 + 4 * gamma0^2)");

      // Logarithmic derivative L'(2, χ_d)/L(2, χ_d)
      m.def("log_derivative", &log_derivative,
            py::arg("chi_arr"), py::arg("lambda_arr"), py::arg("K"),
            "Compute L'(2, χ_d)/L(2, χ_d) from Kronecker and von Mangoldt arrays");

      // Helper function for single zero contribution to the lhs
      m.def("zero_contribution", &zero_contribution,
      py::arg("gamma_minus"), py::arg("gamma_plus"),
      "Contribution of one zero interval to the inequality");
      
      // Main RH verification wrapper
      m.def("rh_verify", &py_rh_verify,
            py::arg("d"), py::arg("K"), py::arg("eta"),
            py::arg("intervals"), py::arg("chi_arr"), py::arg("lambda_arr"),
            R"pbdoc(
Verify the RH inequality for a given discriminant d.

Returns
-------
(success, N_used, lhs, rhs) : tuple
    verified : bool   - True if inequality satisfied
    N_used   : int    - how many zeros (intervals) were used
    lhs/rhs  : float  - evaluated sides of the inequality
)pbdoc");
}
