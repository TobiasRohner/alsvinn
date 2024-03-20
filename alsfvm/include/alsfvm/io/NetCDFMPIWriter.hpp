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
#include <mpi.h>
#include "alsfvm/io/NetCDFWriter.hpp"
namespace alsfvm {
namespace io {

//! Writes to the mpi version of netcdf.
//! @note Due to the new mpi of pNetCDF, this can not be combined
//! in any meaningful way with the old NetCDFWriter class, the code is pretty
//! much disjoint.
//!
//! @note can easily be read by the netcdf python package, see
//!       http://www.hydro.washington.edu/~jhamman/hydro-logic/blog/2013/10/12/plot-netcdf-data/
class NetCDFMPIWriter : public NetCDFWriter {
public:
    ///
    /// \brief NetCDFMPIWriter constructs a new NetCDFMPIWriter
    /// \param basefileName the basefilename to use (this could be eg.
    ///                     "some_simulation".
    /// \param groupNames names of groups to create in the file
    ///        (this is especially useful for MPI). If left blank (""), no prefix will be given
    /// \param groupIndex the groupIndex to write to
    ///
    /// \param newFile should we create (or overwrite) the file? If false,
    ///                the file will be opened and fail if it does not exist
    ///                or if the data does not match (ie. if the sizes mismatch,
    ///                or if the datasets are named differently)
    ///
    /// \param mpiCommunicator the given mpiCommunicator (used for pNETCDF)
    ///
    /// \param mpiInfo the mpiInfo (passed to pNetCDF)
    ///
    /// \note Timestep information will be added to the filename, as well as
    ///       proper extension (.h5).
    ///
    NetCDFMPIWriter(const std::string& filename,
	const grid::Grid &grid,
	size_t num_samples,
	const std::vector<real> &timesteps,
        const std::vector<std::string>& groupNames,
        size_t groupIndex,
        bool newFile,
        MPI_Comm mpiCommunicator,
        MPI_Info mpiInfo);

    //! We could inherit from this, hence virtual destructor.
    virtual ~NetCDFMPIWriter() {}


    ///
    /// \brief write writes the data to disk
    ///
    /// This writes the data in the format
    ///
    /// \code
    ///     <groupName>_<variable_name>
    /// \endcode
    ///
    /// we do not use any groups to store the file at the moment,
    /// this is to ensure maximal compatability with pNetCDF.
    ///
    /// \param conservedVariables the conservedVariables to write
    /// \param grid the grid that is used (describes the _whole_ domain)
    /// \param timestepInformation
    ///
    virtual void write(const volume::Volume& conservedVariables,
        const grid::Grid& grid,
        const simulator::TimestepInformation& timestepInformation) override;

protected:
private:
    const std::vector<std::string> groupNames;
    const size_t groupIndex;
    const bool newFile;
    MPI_Comm mpiCommunicator;
    MPI_Info mpiInfo;

};
} // namespace io
} // namespace alsfvm
