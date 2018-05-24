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

#include "alsfvm/integrator/RungeKutta3.hpp"

namespace alsfvm {
namespace integrator {


RungeKutta3::RungeKutta3(alsfvm::shared_ptr<System> system)
    : system(system) {

}


///
/// Returns the number of substeps this integrator uses.
/// Since this is third order RK, we need three subtimesteps
///
/// \returns 3
///
size_t RungeKutta3::getNumberOfSubsteps() const {
    return 3;
}

///
/// Performs one substep and stores the result to output.
///
/// \param inputConserved should have the output from the previous invocations
///        in this substep, if this is the first invocation, then this will have one element,
///        second timestep 2 elements, etc.
/// \param spatialCellSizes should be the cell size in each direction
/// \param dt is the timestep
/// \param output where to write the output
/// \note the next invocation to performSubstep will get as input the previuosly calculated outputs
///
real RungeKutta3::performSubstep(
    std::vector<alsfvm::shared_ptr<volume::Volume> >& inputConserved,
    rvec3 spatialCellSizes, real dt, real cfl,
    volume::Volume& output, size_t substep,
    const simulator::TimestepInformation& timestepInformation) {
    // We compute U + dt * F(U)


    // Compute F(U)
    rvec3 waveSpeeds(0, 0, 0);
    (*system)(*inputConserved[substep], waveSpeeds, true, output);

    if (substep == 0) {
        dt = computeTimestep(waveSpeeds, spatialCellSizes, cfl, timestepInformation);
    }

    rvec3 cellScaling(dt / spatialCellSizes.x,
        dt / spatialCellSizes.y,
        dt / spatialCellSizes.z);

    output *= cellScaling.x;

    output += *inputConserved[substep];



    if (substep == 1) {
        output *= 1. / 3.;
        output += *inputConserved[0];
        output *= 3. / 4.;
    } else if (substep == 2) {
        output *= 2;
        output += *inputConserved[0];
        output *= 1. / 3.;
    }

    return dt;
}
}
}


