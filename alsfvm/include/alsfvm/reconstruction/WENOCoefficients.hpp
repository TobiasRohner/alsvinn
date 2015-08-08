#pragma once
#include "alsfvm/types.hpp"
#include <array>
namespace alsfvm { namespace reconstruction { 

	template<int k>
    class WENOCoefficients {
    public:
		static real coefficients[];

		template<size_t index>
		static real computeBeta(const std::array<real, 2 * k - 1>& stencil);

		
    };

	template<>
	template<size_t index>
	 real WENOCoefficients<2>::computeBeta(const std::array<real, 3>& V) {
		if (index == 0) {
			real beta = V[2] - V[1];
			return beta*beta;
		}
		else if (index == 1) {
			real beta = V[1] - V[0];
			return beta*beta;
		}

		static_assert(index < 2, "Only up to index 1 for order 2 in WENO");
		return 0;
	}

	template<>
	template<size_t index>
	real WENOCoefficients<3>::computeBeta(const std::array<real, 5>& V) {
		if (index == 0) {
			return  13.0 / 12.0 * pow(V[2] - 2 * V[3] + V[4], 2) + 1 / 4.0*pow(3 * V[2] - 4 * V[3] + V[4], 2);
		}
		else if (index == 1) {
			return  13.0 / 12.0 * pow(V[1] - 2 * V[2] + V[3], 2) + 1 / 4.0*pow(V[1] -  V[3], 2);
		}
		else if (index == 2) {
			return  13.0 / 12.0 * pow(V[0] - 2 * V[1] + V[2], 2) + 1 / 4.0*pow(V[0] - 4 * V[1] + 3*V[2], 2);
		}
		static_assert(index < 3, "Only up to index 1 for order 2 in WENO");
		return 0;
	}





} // namespace alsfvm
} // namespace reconstruction