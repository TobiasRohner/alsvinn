/* Copyright (c) 2018 ETH Zurich, Kjetil Olsen Lye
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "alsfvm/simulator/Simulator.hpp"
#include "alsuq/run/SimulatorCreator.hpp"
#include "alsfvm/init/Parameters.hpp"
#include "alsuq/mpi/Configuration.hpp"
#include <mpi.h>
#include "alsuq/types.hpp"

namespace alsuq {
namespace run {
//!
//! \brief The FiniteVolumeSimulatorCreator class creates a new instance of the FVM simulator
//!
class FiniteVolumeSimulatorCreator : public SimulatorCreator {
public:
    FiniteVolumeSimulatorCreator(const std::string& configurationFile,
	size_t numberOfSamples,
        mpi::ConfigurationPtr mpiConfigurationSpatial,
        mpi::ConfigurationPtr mpiConfigurationStatistical,
        mpi::ConfigurationPtr mpiConfigurationWorld,
        ivec3 multiSpatial
    );

    alsfvm::shared_ptr<alsfvm::simulator::AbstractSimulator>
    createSimulator(const alsfvm::init::Parameters& initialDataParameters,
        size_t sampleNumber) override;

private:
    mpi::ConfigurationPtr mpiConfigurationSpatial;
    mpi::ConfigurationPtr mpiConfigurationStatistical;
    mpi::ConfigurationPtr mpiConfigurationWorld;

    ivec3 multiSpatial;

    //! Gathers all the current samples from all current mpi procs
    //! and creates a list of names of the samples now being computed
    std::vector<std::string> makeGroupNames(size_t sampleNumber);

    bool firstCall{true};
    const std::string filename;



};
} // namespace run
} // namespace alsuq
