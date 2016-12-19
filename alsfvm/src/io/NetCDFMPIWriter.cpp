#include "alsfvm/io/NetCDFMPIWriter.hpp"
#include <pnetcdf.h>
#include "alsutils/log.hpp"
namespace alsfvm { namespace io {

NetCDFMPIWriter::NetCDFMPIWriter(const std::__cxx11::string &basefileName,
                                 const std::vector<std::string> &groupNames,
                                 size_t groupIndex, bool newFile,
                                 MPI_Comm mpiCommunicator, MPI_Info mpiInfo)
    : NetCDFWriter(basefileName),
      groupNames(groupNames),
      groupIndex(groupIndex),
      newFile(newFile),
      mpiCommunicator(mpiCommunicator),
      mpiInfo(mpiInfo)
{

}

void NetCDFMPIWriter::write(const volume::Volume &conservedVariables, const volume::Volume &extraVariables, const grid::Grid &grid, const simulator::TimestepInformation &timestepInformation)
{
    netcdf_raw_ptr file;
    auto filename = getFilename();

    if (newFile) {

        NETCDF_SAFE_CALl(ncmpi_create(mpiCommunicator, filename.c_str(), NC_NETCDF4|NC_MPIIO,
                                   mpiInfo, &file));
    }
    else {
        NETCDF_SAFE_CALl(ncmpi_open(mpiCommunicator, filename.c_str(), NC_WRITE|NC_NETCDF4|NC_MPIIO,
                            mpiInfo, &file));
        NETCDF_SAFE_CALl(ncmpi_redef(file));
    }
    writeToFile(file, conservedVariables, extraVariables,
                grid, timestepInformation, newFile);


    NETCDF_SAFE_CALl(ncmpi_close(file));
}

NetCDFMPIWriter::dimension_vector NetCDFMPIWriter::createDimensions(netcdf_raw_ptr baseGroup, const volume::Volume &volume, bool newFile)
{
    std::array<netcdf_raw_ptr,3> dimensions;
    netcdf_raw_ptr xdim, ydim, zdim;

    if (newFile) {
        NETCDF_SAFE_CALl(ncmpi_def_dim(baseGroup, "x", volume.getNumberOfXCells(),
                                       &xdim));
        NETCDF_SAFE_CALl(ncmpi_def_dim(baseGroup, "y", volume.getNumberOfYCells(),
                                       &ydim));
        NETCDF_SAFE_CALl(ncmpi_def_dim(baseGroup, "z", volume.getNumberOfZCells(),
                                       &zdim));
    } else {
        NETCDF_SAFE_CALl(ncmpi_inq_dimid(baseGroup, "x", &xdim));
        NETCDF_SAFE_CALl(ncmpi_inq_dimid(baseGroup, "y", &ydim));
        NETCDF_SAFE_CALl(ncmpi_inq_dimid(baseGroup, "z", &zdim));
    }

    dimensions[0] = xdim;
    dimensions[1] = ydim;
    dimensions[2] = zdim;

    return dimensions;
}

std::vector<netcdf_raw_ptr>
NetCDFMPIWriter::makeDataset(netcdf_raw_ptr baseGroup,
                             const volume::Volume &volume,
                             std::array<netcdf_raw_ptr, 3> dimensions)
{
    std::vector<netcdf_raw_ptr> datasets;


    for (const auto& groupName : groupNames) {

        for (size_t memoryIndex = 0; memoryIndex < volume.getNumberOfVariables();
             ++memoryIndex) {
            netcdf_raw_ptr dataset;

            auto memoryName = groupName + "_" + volume.getName(memoryIndex) ;


            NETCDF_SAFE_CALl(ncmpi_def_var(baseGroup, memoryName.c_str(), NC_DOUBLE, 3,
                                           dimensions.data(), &dataset));

            if (groupName == groupNames[groupIndex]) {
                datasets.push_back(dataset);
            }
        }

    }

    return datasets;
}

void NetCDFMPIWriter::writeToFile(netcdf_raw_ptr file,
                                  const volume::Volume &conservedVariables,
                                  const volume::Volume &extraVariables,
                                  const grid::Grid &grid,
                                  const simulator::TimestepInformation &timestepInformation,
                                  bool newFile)
{


    auto dimensions = createDimensions(file, conservedVariables, newFile);


    auto datasetsConserved = makeDataset(file, conservedVariables, dimensions);
    auto datasetsExtra = makeDataset(file, extraVariables, dimensions);

    NETCDF_SAFE_CALl(ncmpi_enddef(file));
    writeVolume(file, conservedVariables, dimensions, datasetsConserved);
    writeVolume(file, extraVariables, dimensions, datasetsExtra);

}

void NetCDFMPIWriter::writeMemory(netcdf_raw_ptr baseGroup, netcdf_raw_ptr dataset, const volume::Volume &volume, size_t memoryIndex)
{
    std::vector<real> dataTmp(volume.getNumberOfXCells() * volume.getNumberOfYCells() * volume.getNumberOfZCells());

    volume.copyInternalCells(memoryIndex, dataTmp.data(), dataTmp.size());

    std::vector<double> data(dataTmp.size());
    std::copy(dataTmp.begin(), dataTmp.end(), data.begin());
    NETCDF_SAFE_CALl(ncmpi_put_var_double_all(baseGroup, dataset, data.data()));
}

void NetCDFMPIWriter::writeVolume(netcdf_raw_ptr baseGroup, const volume::Volume &volume, std::array<netcdf_raw_ptr, 3> dimensions,
                                  const std::vector<netcdf_raw_ptr>& datasets)
{
    for (size_t memoryIndex = 0; memoryIndex < volume.getNumberOfVariables();
                 ++memoryIndex) {
        auto dataset = datasets[memoryIndex];

        writeMemory(baseGroup, dataset, volume, memoryIndex);
    }
}


}
}