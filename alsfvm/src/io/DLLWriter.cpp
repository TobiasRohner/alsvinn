#include "alsfvm/io/DLLWriter.hpp"
#include <boost/algorithm/algorithm.hpp>
#include <boost/dll.hpp>
#include <boost/algorithm/string.hpp>
#include "alsutils/config.hpp"
#include "alsutils/cuda/get_current_gpu_id.hpp"

namespace alsfvm {
namespace io {

DLLWriter::DLLWriter(const std::string& basename,
    const Parameters& parameters,
    alsutils::mpi::ConfigurationPtr mpiConfigration) {


    auto filename = parameters.getString("library");
    auto createFunctionName = parameters.getString("create_function");
    auto makeParametersName = parameters.getString("make_parameters_function");
    auto setMpiCommName = parameters.getString("set_mpi_comm_function");

    if (boost::algorithm::to_lower_copy(makeParametersName) != "none") {
        auto makeParametersFunction = boost::dll::import_symbol <void* ()>(filename,
                makeParametersName);
        parametersStruct = makeParametersFunction();

        auto setParameterFunctionName = parameters.getString("set_parameter_function");

        auto setParameterFunction = boost::dll::import_symbol
            <void(void*, const char*, const char*)>(filename,
                setParameterFunctionName);

        auto deleteParametersFunctionName =
            parameters.getString("delete_parameters_function");

        deleteParametersFunction = boost::dll::import_symbol
            <void(void*)>(filename,
                deleteParametersFunctionName);

        for (auto key : parameters.getKeys()) {

            setParameterFunction(parametersStruct, key.c_str(),
                parameters.getString(key).c_str());
        }
    }

    auto createFunction =
        boost::dll::import_symbol<void* (const char*, const char*, void*)>(filename,
            createFunctionName);

    dllData = createFunction("alsvinn",
            (std::string("https://github.com/alsvinn/alsvinn git") +
                alsutils::getVersionControlID()).c_str(),
            parametersStruct
        );

    auto newTimestepFunctionName =
        parameters.getString("new_timestep_function");

    if (boost::algorithm::to_lower_copy(newTimestepFunctionName) != "none") {
        newTimestepFunction = boost::dll::import_symbol<void(DLLData, DLLData, real, int)>
            (filename,
                newTimestepFunctionName);
    }

    auto endTimestepFunctionName =
        parameters.getString("end_timestep_function");


    if (boost::algorithm::to_lower_copy(endTimestepFunctionName) != "none") {
        endTimestepFunction = boost::dll::import_symbol<void(DLLData, DLLData, real, int)>
            (filename,
                endTimestepFunctionName);

    }

    auto setMpiCommFunctionName = parameters.getString("set_mpi_comm_function");

    if (boost::algorithm::to_lower_copy(setMpiCommFunctionName) != "none") {
        auto setMpiCommFunction = boost::dll::import_symbol<void(DLLData, DLLData, MPI_Comm)>
            (filename, setMpiCommFunctionName);

        setMpiCommFunction(dllData, parametersStruct,
            mpiConfigration->getCommunicator());
    }

    auto needsDataOnHostFunctionName =
        parameters.getString("needs_data_on_host_function");

    if (boost::algorithm::to_lower_copy(needsDataOnHostFunctionName) != "none") {
        auto needsDataOnHostFunction = boost::dll::import_symbol<bool(void*, void*)>(filename,
                needsDataOnHostFunctionName);

        needsDataOnHost = needsDataOnHostFunction(dllData, parametersStruct);
    }

    const auto writeFunctionName = parameters.getString("write_function");
    writeFunction = boost::dll::import_symbol<write_function_t>(filename,
            writeFunctionName);
}

void DLLWriter::write(const volume::Volume& conservedVariables,
    const grid::Grid& grid,
    const simulator::TimestepInformation& timestepInformation) {



    if (newTimestepFunction) {
        newTimestepFunction(dllData, parametersStruct,
            timestepInformation.getCurrentTime(),
            timestepInformation.getNumberOfStepsPerformed());
    }

    const real ax = grid.getOrigin().x;
    const real ay = grid.getOrigin().y;
    const real az = grid.getOrigin().z;

    const real bx = grid.getTop().x;
    const real by = grid.getTop().y;
    const real bz = grid.getTop().z;

    const int ngx = int(conservedVariables.getNumberOfXGhostCells());
    const int ngy = int(conservedVariables.getNumberOfYGhostCells());
    const int ngz = int(conservedVariables.getNumberOfZGhostCells());


    const int nx = int(conservedVariables.getNumberOfXCells());
    const int ny = int(conservedVariables.getNumberOfYCells());
    const int nz = int(conservedVariables.getNumberOfZCells());


    for (size_t var = 0; var < conservedVariables.getNumberOfVariables(); ++var) {
        auto dataSmartPointer =  conservedVariables.getScalarMemoryArea(var);

        if (needsDataOnHost && !dataSmartPointer->isOnHost()) {
            dataSmartPointer = dataSmartPointer->getHostMemory();

        }

        int gpuID = -1;

        if (!dataSmartPointer->isOnHost()) {
            gpuID = alsutils::cuda::getCurrentGPUId();
        }

        writeFunction(dllData, parametersStruct,
            timestepInformation.getCurrentTime(),
            conservedVariables.getName(var).c_str(),
            dataSmartPointer->getPointer(),
            nx, ny, nz,
            ngx, ngy, ngz,
            ax, ay, az,
            bx, by, bz,
            gpuID);
    }

    if (endTimestepFunction) {
        endTimestepFunction(dllData, parametersStruct,
            timestepInformation.getCurrentTime(),
            timestepInformation.getNumberOfStepsPerformed());
    }

    // remember, signature is
    //real write_function(void* data,
    //  void* parameters,
    //  real time,
    //  const char* variable_name,
    //  const real* variable_data,
    //  int nx, int ny, int nz,
    //  int ngx, int ngy, int ngz,
    //  real ax, real ay, real az,
    //  real bx, real by, real bz,
    //  int gpu_number );




}

void DLLWriter::finalize(const grid::Grid& grid,
    const simulator::TimestepInformation& timestepInformation) {

    if (deleteFunction) {
        deleteFunction(dllData);
    }

    if (deleteParametersFunction) {
        deleteParametersFunction(parametersStruct);
    }
}

}
}
