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

#include "alsfvm/io/NetCDFMPIWriter.hpp"
#include <pnetcdf.h>
#include "alsutils/log.hpp"
#include "alsutils/mpi/to_mpi_offset.hpp"
#include "alsfvm/io/parallel_netcdf_write_report.hpp"
#include "alsfvm/io/parallel_netcdf_write_attributes.hpp"
#include <boost/filesystem.hpp>
#include "alsutils/timer/Timer.hpp"
#include "alsfvm/io/parallel_netcdf_utils.hpp"

#include <fstream>

namespace alsfvm {
namespace io {

NetCDFMPIWriter::NetCDFMPIWriter(const std::string& filename,
    const grid::Grid &grid,
    size_t num_samples,
    const std::vector<real> &timesteps,
    const std::vector<std::string>& groupNames,
    size_t groupIndex, bool newFile,
    MPI_Comm mpiCommunicator, MPI_Info mpiInfo)
    : NetCDFWriter(filename, newFile, grid, num_samples, timesteps),
      groupNames(groupNames),
      groupIndex(groupIndex),
      newFile(newFile),
      mpiCommunicator(mpiCommunicator),
      mpiInfo(mpiInfo) {
}

void NetCDFMPIWriter::write(const volume::Volume& conservedVariables,
    const grid::Grid& grid,
    const simulator::TimestepInformation& timestepInformation) {

    ALSVINN_TIME_BLOCK(alsvinn, fvm, io, netcdf);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const size_t sampleIdx = std::stoll(&groupNames[rank][7]);

    auto globalPosition = alsutils::mpi::to_mpi_offset(grid.getGlobalPosition());
    auto localSize = alsutils::mpi::to_mpi_offset(grid.getDimensions());
    if (grid.getActiveDimension() == 2) {
        std::swap(globalPosition[0], globalPosition[1]);
        std::swap(localSize[0], localSize[1]);
    }
    if (grid.getActiveDimension() == 3) {
        std::swap(globalPosition[0], globalPosition[2]);
        std::swap(localSize[0], localSize[2]);
    }

    size_t start[5];
    size_t count[5];
    start[0] = sampleIdx;
    count[0] = 1;
    start[1] = snapshotNumber;
    count[1] = 1;
    start[2] = globalPosition[0];
    count[2] = localSize[0];
    start[3] = globalPosition[1];
    count[3] = localSize[1];
    start[4] = globalPosition[2];
    count[4] = localSize[2];
    std::cout << "start = {" << start[0] << ", " << start[1] << ", " << start[2] << ", " << start[3] << ", " << start[4] << "}" << std::endl;
    std::cout << "count = {" << count[0] << ", " << count[1] << ", " << count[2] << ", " << count[3] << ", " << count[4] << "}" << std::endl;
    const size_t slice_size = localSize[0] * localSize[1] * localSize[2];
    std::vector<real> dataTmp(slice_size);
    conservedVariables.copyInternalCells(0, dataTmp.data(), dataTmp.size());
    NETCDF_SAFE_CALL(
	nc_put_vara(ncid_, varids_[0], start, count, dataTmp.data()));
    conservedVariables.copyInternalCells(conservedVariables.getNumberOfVariables() - 1, dataTmp.data(), dataTmp.size());
    NETCDF_SAFE_CALL(
	nc_put_vara(ncid_, varids_[1], start, count, dataTmp.data()));
    conservedVariables.copyInternalCells(1, dataTmp.data(), dataTmp.size());
    NETCDF_SAFE_CALL(
	nc_put_vara(ncid_, varids_[2], start, count, dataTmp.data()));
    conservedVariables.copyInternalCells(2, dataTmp.data(), dataTmp.size());
    NETCDF_SAFE_CALL(
	nc_put_vara(ncid_, varids_[3], start, count, dataTmp.data()));
    if (conservedVariables.getNumberOfVariables() == 5) {
      conservedVariables.copyInternalCells(3, dataTmp.data(), dataTmp.size());
      NETCDF_SAFE_CALL(
	  nc_put_vara(ncid_, varids_[4], start, count, dataTmp.data()));
    }
    ++snapshotNumber;
}

}
}
