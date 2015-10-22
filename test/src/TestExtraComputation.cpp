#include <gtest/gtest.h>

#include "alsfvm/equation/CellComputerFactory.hpp"
#include "alsfvm/volume/VolumeFactory.hpp"
#include "alsfvm/volume/volume_foreach.hpp"
#include "alsfvm/equation/euler/AllVariables.hpp"
#include "alsfvm/equation/euler/Euler.hpp"


using namespace alsfvm;
using namespace alsfvm::memory;
using namespace alsfvm::equation;
using namespace alsfvm::volume;

class TestExtraComputation : public ::testing::Test {
public:
    alsfvm::shared_ptr<DeviceConfiguration> deviceConfiguration;
    const std::string equation;
    const std::string platform;
    const size_t nx;
    const size_t ny;
    const size_t nz;
    alsfvm::shared_ptr<MemoryFactory> memoryFactory;
    VolumeFactory volumeFactory;
    alsfvm::shared_ptr<Volume> conservedVolume;
    alsfvm::shared_ptr<Volume> extraVolume;
    alsfvm::shared_ptr<simulator::SimulatorParameters> simulatorParameters;
    CellComputerFactory cellComputerFactory;
    TestExtraComputation()
        : deviceConfiguration(new DeviceConfiguration("cpu")),
          equation("euler"),
          platform("cpu"),
          nx(10), ny(10), nz(10),
          memoryFactory(new MemoryFactory(deviceConfiguration)),
          volumeFactory("euler", memoryFactory),
          conservedVolume(volumeFactory.createConservedVolume(nx, ny, nz)),
          extraVolume(volumeFactory.createExtraVolume(nx, ny, nz)),
          simulatorParameters(new simulator::SimulatorParameters("euler", "cpu")),
          cellComputerFactory(simulatorParameters, deviceConfiguration)
    {

    }
};

TEST_F(TestExtraComputation, ConstructTest) {
    auto cellComputer =cellComputerFactory.createComputer();
}

TEST_F(TestExtraComputation, CheckExtraCalculation) {
	auto cellComputer = cellComputerFactory.createComputer();

	// Fill up volume
	transform_volume<euler::ConservedVariables, euler::ConservedVariables>
		(*conservedVolume, *conservedVolume, [](const euler::ConservedVariables& in) -> euler::ConservedVariables {
		return euler::ConservedVariables(0.5, 1, 1, 1, 4.4);
	});

	cellComputer->computeExtraVariables(*conservedVolume, *extraVolume);

    auto& eulerParameters = static_cast<euler::EulerParameters&> (simulatorParameters->getEquationParameters());
    for_each_cell<euler::ExtraVariables>(*extraVolume, [&](const euler::ExtraVariables& in, size_t index) {
		ASSERT_EQ(in.u, rvec3(2, 2, 2));
        ASSERT_EQ(in.p, (eulerParameters.getGamma() - 1)*(4.4 - 0.5 * 3 / 0.5));
	});

}

TEST_F(TestExtraComputation, CheckMaximumWaveSpeed) {
	auto cellComputer = cellComputerFactory.createComputer();
    auto& eulerParameters = static_cast<euler::EulerParameters&> (simulatorParameters->getEquationParameters());
	// Fill up volume
	transform_volume<euler::ConservedVariables, euler::ConservedVariables>
		(*conservedVolume, *conservedVolume, [](const euler::ConservedVariables& in) -> euler::ConservedVariables {
		return euler::ConservedVariables(0.5, 1, 1, 1, 4.4);
	});

    const real gamma = eulerParameters.getGamma();
    cellComputer->computeExtraVariables(*conservedVolume, *extraVolume);
    {
        real maxWaveSpeed = cellComputer->computeMaxWaveSpeed(*conservedVolume, *extraVolume,0);

        ASSERT_EQ(maxWaveSpeed, 2 + sqrt(gamma * (gamma - 1)*(4.4 - 0.5 * 3 / 0.5) / 0.5));
    }

    {
        real maxWaveSpeed = cellComputer->computeMaxWaveSpeed(*conservedVolume, *extraVolume, 1);

        ASSERT_EQ(maxWaveSpeed, 2 + sqrt(gamma * (gamma - 1)*(4.4 - 0.5 * 3 / 0.5) / 0.5));
    }

    {
        real maxWaveSpeed = cellComputer->computeMaxWaveSpeed(*conservedVolume, *extraVolume, 2);

        ASSERT_EQ(maxWaveSpeed, 2 + sqrt(gamma * (gamma - 1)*(4.4 - 0.5 * 3 / 0.5) / 0.5));
    }
}


TEST_F(TestExtraComputation, CheckConstraints) {
	auto cellComputer = cellComputerFactory.createComputer();

	// Fill up volume
	transform_volume<euler::ConservedVariables, euler::ConservedVariables>
		(*conservedVolume, *conservedVolume, [](const euler::ConservedVariables& in) -> euler::ConservedVariables {
		return euler::ConservedVariables(0.5, 1, 1, 1, 4.4);
	});

	cellComputer->computeExtraVariables(*conservedVolume, *extraVolume);

	// This should be fine
	ASSERT_TRUE(cellComputer->obeysConstraints(*conservedVolume, *extraVolume));

	// Now we fill it with something that cancels the constraints
	conservedVolume->getScalarMemoryArea("rho")->getPointer()[4] = -0.4;

	ASSERT_FALSE(cellComputer->obeysConstraints(*conservedVolume, *extraVolume));

	// and some extra
	extraVolume->getScalarMemoryArea("p")->getPointer()[8] = -0.4;
	ASSERT_FALSE(cellComputer->obeysConstraints(*conservedVolume, *extraVolume));

	// And fix the first one, then we should still get something false
    conservedVolume->getScalarMemoryArea("rho")->getPointer()[4] = 2;
	ASSERT_FALSE(cellComputer->obeysConstraints(*conservedVolume, *extraVolume));

    // check that it can be fixed again
    extraVolume->getScalarMemoryArea("p")->getPointer()[8] = 0.4;
    ASSERT_TRUE(cellComputer->obeysConstraints(*conservedVolume, *extraVolume));


    // Add an inf
    extraVolume->getScalarMemoryArea("p")->getPointer()[8] = INFINITY;
    ASSERT_FALSE(cellComputer->obeysConstraints(*conservedVolume, *extraVolume));

}
