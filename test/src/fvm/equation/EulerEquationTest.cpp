#include "gtest/gtest.h"
#include "alsfvm/equation/euler/Euler.hpp"

using namespace alsfvm::equation::euler;
using namespace alsfvm;
TEST(EulerEquationTest, FluxTest) {
	// This test checks that the point flux is correction setup

	// First we check that we get the correct output if everything is one

    AllVariables<3> input(1, rvec3{ 1, 1, 1 }, 1, 1, rvec3{ 1, 1, 1 });

    EulerParameters parameters;
    Euler<3> equation(parameters);
	{
        ConservedVariables<3> output(0, rvec3{ 0, 0, 0 }, 0);

        equation.computePointFlux < 0 >(input, output);

		ASSERT_EQ(output.E, 2);
		ASSERT_EQ(output.m, rvec3(2, 1, 1));
		ASSERT_EQ(output.rho, 1);
	}

		
	{
        ConservedVariables<3> output(0, rvec3{ 0, 0, 0 }, 0);
        equation.computePointFlux < 1 >(input, output);
		ASSERT_EQ(output.E, 2);
		ASSERT_EQ(output.m, rvec3(1, 2, 1));
		ASSERT_EQ(output.rho, 1);
	}

	{
        ConservedVariables<3> output(0, rvec3{ 0, 0, 0 }, 0);

        equation.computePointFlux < 2 >(input, output);

		ASSERT_EQ(output.E, 2);
		ASSERT_EQ(output.m, rvec3(1, 1, 2));
		ASSERT_EQ(output.rho, 1);
	}
}