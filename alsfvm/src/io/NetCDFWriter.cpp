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

#include "alsfvm/io/NetCDFWriter.hpp"
#include "alsfvm/io/netcdf_utils.hpp"
#include "alsfvm/io/io_utils.hpp"
#include "alsutils/log.hpp"
#include "alsfvm/io/netcdf_write_report.hpp"
#include "alsfvm/io/netcdf_write_attributes.hpp"
#include "alsutils/timer/Timer.hpp"
#include "alsutils/log.hpp"
#if ALSVINN_USE_MPI
#include <netcdf_par.h>
#endif


namespace alsfvm {
namespace io {

static constexpr nc_type netcdf_type(float) { return NC_FLOAT; }
static constexpr nc_type netcdf_type(double) { return NC_DOUBLE; }

NetCDFWriter::NetCDFWriter(const std::string& filename, bool newFile, const grid::Grid &grid, size_t num_samples, const std::vector<real> &timesteps) : basefileName(filename) {
  const ivec3 grid_dims = grid.getDimensions();
  const int dim = grid_dims.z > 1 ? 3 : 2;

  if (newFile) {
#if ALSVINN_USE_MPI
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1) {
      NETCDF_SAFE_CALL(nc_create_par(filename.c_str(), NC_CLOBBER | NC_NETCDF4, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid_));
    } else {
      NETCDF_SAFE_CALL(nc_create(filename.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid_));
    }
#else
    NETCDF_SAFE_CALL(nc_create(filename.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid_));
#endif
    NETCDF_SAFE_CALL(nc_set_fill(ncid_, NC_NOFILL, NULL));

    NETCDF_SAFE_CALL(nc_def_dim(ncid_, "member", num_samples, &dimids_[0]));
    NETCDF_SAFE_CALL(nc_def_dim(ncid_, "time", timesteps.size(), &dimids_[1]));
    NETCDF_SAFE_CALL(nc_def_dim(ncid_, "x", grid_dims.x, &dimids_[2]));
    NETCDF_SAFE_CALL(nc_def_dim(ncid_, "y", grid_dims.y, &dimids_[3]));
    if (dim == 3) {
      NETCDF_SAFE_CALL(nc_def_dim(ncid_, "z", grid_dims.z, &dimids_[4]));
    }

    int varid_member, varid_time, varid_x, varid_y, varid_z;
    NETCDF_SAFE_CALL(
	nc_def_var(ncid_, "member", NC_INT, 1, &dimids_[0], &varid_member));
    NETCDF_SAFE_CALL(nc_def_var(
	ncid_, "time", netcdf_type(real{}), 1, &dimids_[1], &varid_time));
    NETCDF_SAFE_CALL(
	nc_def_var(ncid_, "x", netcdf_type(real{}), 1, &dimids_[2], &varid_x));
    NETCDF_SAFE_CALL(
	nc_def_var(ncid_, "y", netcdf_type(real{}), 1, &dimids_[3], &varid_y));
    if (dim == 3) {
      NETCDF_SAFE_CALL(nc_def_var(
	  ncid_, "z", netcdf_type(real{}), 1, &dimids_[4], &varid_z));
    }
    size_t chunksizes[5];
    chunksizes[0] = 1;
    chunksizes[1] = 1;
    chunksizes[2] = grid_dims.x;
    chunksizes[3] = grid_dims.y;
    chunksizes[4] = grid_dims.z;
    NETCDF_SAFE_CALL(nc_def_var(
	ncid_, "rho", netcdf_type(real{}), 2 + dim, dimids_, &varids_[0]));
    NETCDF_SAFE_CALL(nc_def_var_chunking(ncid_, varids_[0], NC_CHUNKED, chunksizes));
    NETCDF_SAFE_CALL(nc_def_var(
	ncid_, "E", netcdf_type(real{}), 2 + dim, dimids_, &varids_[1]));
    NETCDF_SAFE_CALL(nc_def_var_chunking(ncid_, varids_[1], NC_CHUNKED, chunksizes));
    NETCDF_SAFE_CALL(nc_def_var(
	ncid_, "mx", netcdf_type(real{}), 2 + dim, dimids_, &varids_[2]));
    NETCDF_SAFE_CALL(nc_def_var_chunking(ncid_, varids_[2], NC_CHUNKED, chunksizes));
    NETCDF_SAFE_CALL(nc_def_var(
	ncid_, "my", netcdf_type(real{}), 2 + dim, dimids_, &varids_[3]));
    NETCDF_SAFE_CALL(nc_def_var_chunking(ncid_, varids_[3], NC_CHUNKED, chunksizes));
    if (dim == 3) {
      NETCDF_SAFE_CALL(nc_def_var(
	  ncid_, "mz", netcdf_type(real{}), 2 + dim, dimids_, &varids_[4]));
      NETCDF_SAFE_CALL(
	  nc_def_var_chunking(ncid_, varids_[4], NC_CHUNKED, chunksizes));
    }

    std::vector<int> member(num_samples);
    for (int i = 0; i < num_samples; ++i) {
      member[i] = i;
    }
    NETCDF_SAFE_CALL(nc_put_var(ncid_, varid_member, member.data()));
    NETCDF_SAFE_CALL(nc_put_var(ncid_, varid_time, timesteps.data()));
    const std::vector<rvec3> &midpoints = grid.getCellMidpoints();
    std::vector<real> x(grid_dims.x);
    for (int i = 0 ; i < grid_dims.x ; ++i ) {
      x[i] = midpoints[i].x;
    }
    std::vector<real> y(grid_dims.y);
    for (int i = 0 ; i < grid_dims.y ; ++i) {
      y[i] = midpoints[grid_dims.x * i].y;
    }
    std::vector<real> z(grid_dims.z);
    for (int i = 0 ; i < grid_dims.z ; ++i) {
      z[i] = midpoints[grid_dims.x * grid_dims.y * i].z;
    }
    NETCDF_SAFE_CALL(nc_put_var(ncid_, varid_x, x.data()));
    NETCDF_SAFE_CALL(nc_put_var(ncid_, varid_y, y.data()));
    if (dim == 3) {
      NETCDF_SAFE_CALL(nc_put_var(ncid_, varid_z, z.data()));
    }
  } else {
#if ALSVINN_USE_MPI
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1) {
      NETCDF_SAFE_CALL(nc_open_par(filename.c_str(), NC_WRITE, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid_));
    } else {
      NETCDF_SAFE_CALL(nc_open(filename.c_str(), NC_WRITE, &ncid_));
    }
#else
    NETCDF_SAFE_CALL(nc_open(filename.c_str(), NC_WRITE, &ncid_));
#endif
    NETCDF_SAFE_CALL(nc_inq_varid(ncid_, "rho", &varids_[0]));
    NETCDF_SAFE_CALL(nc_inq_varid(ncid_, "E", &varids_[1]));
    NETCDF_SAFE_CALL(nc_inq_varid(ncid_, "mx", &varids_[2]));
    NETCDF_SAFE_CALL(nc_inq_varid(ncid_, "my", &varids_[3]));
    if (dim == 3) {
      NETCDF_SAFE_CALL(nc_inq_varid(ncid_, "mz", &varids_[4]));
    }
  }
}

void NetCDFWriter::write(const volume::Volume& conservedVariables,
    const grid::Grid& grid,
    const simulator::TimestepInformation& timestepInformation) {
    ALSVINN_TIME_BLOCK(alsvinn, fvm, io, netcdf);

    const size_t sampleIdx = std::stoll(&conservedVariables.getName(0)[7]);

    size_t start[5];
    size_t count[5];
    start[0] = sampleIdx;
    count[0] = 1;
    start[1] = snapshotNumber;
    count[1] = 1;
    start[2] = 0;
    count[2] = conservedVariables.getNumberOfXCells();
    start[3] = 0;
    count[3] = conservedVariables.getNumberOfYCells();
    start[4] = 0;
    count[4] = conservedVariables.getNumberOfZCells();
    const size_t slice_size = conservedVariables.getNumberOfXCells() * conservedVariables.getNumberOfYCells() * conservedVariables.getNumberOfZCells();
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
