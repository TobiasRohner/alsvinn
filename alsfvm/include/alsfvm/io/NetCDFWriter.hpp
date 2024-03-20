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
#include "alsfvm/io/Writer.hpp"
#include <netcdf.h>
#include "alsfvm/io/netcdf_utils.hpp"

namespace alsfvm {
namespace io {

//! The netcdf writer writes to the netcdf format.
//! This is the recommended writer to use
//!
//! @note can easily be read by the netcdf python package, see
//!       http://www.hydro.washington.edu/~jhamman/hydro-logic/blog/2013/10/12/plot-netcdf-data/
class NetCDFWriter : public Writer {
public:
    //! Creates a new instance of the NetCDFWriter
    //!
    //! @param basefileName the base filename to use. Resulting filenames
    //!                     will be of the form
    //!                         basefileName_<timestep>.nc
    //!
    NetCDFWriter(const std::string& filename, bool newFile, const grid::Grid &grid, size_t num_samples, const std::vector<real> &timesteps);

    //! Since we inherit from this class, this is the safest option
    //! (the destructor is anyway empty)
    virtual ~NetCDFWriter() {}

    //! Write the volume to file
    //! This will create a variable for each volume
    //!
    //! There will be no additioanl grid information written,
    //! this is implicit in the netcdf format
    virtual void write(const volume::Volume& conservedVariables,
        const grid::Grid& grid,
        const simulator::TimestepInformation& timestepInformation) override;

protected:
  int ncid_;
  int dimids_[5];
  int varids_[5];

    size_t snapshotNumber{0};
    const std::string basefileName;
};
} // namespace io
} // namespace alsfvm
