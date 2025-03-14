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

#include "alsuq/config/Setup.hpp"
#include "alsuq/run/FiniteVolumeSimulatorCreator.hpp"
#include "alsuq/generator/GeneratorFactory.hpp"
#include "alsuq/distribution/DistributionFactory.hpp"
#include "alsuq/mpi/SimpleLoadBalancer.hpp"
#include <boost/property_tree/xml_parser.hpp>
#include <sstream>
#include "alsutils/io/TextFileCache.hpp"
#include "alsuq/stats/StatisticsFactory.hpp"
#include "alsfvm/io/WriterFactory.hpp"
#include "alsfvm/io/MpiWriterFactory.hpp"

#include "alsuq/stats/FixedIntervalStatistics.hpp"
#include "alsuq/stats/TimeIntegratedWriter.hpp"
#include <boost/algorithm/string.hpp>
#include "alsutils/log.hpp"
#include "alsutils/timer/Timer.hpp"
#include <boost/lexical_cast.hpp>

namespace alsuq {
namespace config {
// example:
// <samples>1024</samples>
// <generator>auto</generator>
// <parameters>
//   <parameter>
//     <name>a</name>
//     <length>40</length>
//     <type>uniform</type>
//   </parameter>
// </parameters>
// <stats>
//   <stat>
//     meanvar
//   </stat>
// </stats>

namespace {
template<class T>
vec3<T> parseVector(const std::string& vectorAsString) {
    std::vector<std::string> splitString;
    boost::split(splitString, vectorAsString, boost::is_any_of("\t "));

    T a = boost::lexical_cast<T>(splitString[0]);
    T b = boost::lexical_cast<T>(splitString[1]);
    T c = boost::lexical_cast<T>(splitString[2]);

    return vec3<T>(a, b, c);

}
}


std::shared_ptr<run::Runner> Setup::makeRunner(const std::string& inputFilename,
    alsutils::mpi::ConfigurationPtr mpiConfigurationWorld,
    int multiSample, ivec3 multiSpatial) {
    ALSVINN_TIME_BLOCK(alsvinn, uq, init);
    auto& textCache = alsutils::io::TextFileCache::getInstance();
    auto configurationFileContent = textCache.loadTextFile(inputFilename);
    std::stringstream stream(configurationFileContent);

    ptree configurationBase;
    boost::property_tree::read_xml(stream, configurationBase);
    auto configuration = configurationBase.get_child("config");
    auto sampleGenerator = makeSampleGenerator(configuration);
    auto numberOfSamples = readNumberOfSamples(configuration);
    auto sampleStart = readSampleStart(configuration);

    std::cout << "numberOfSamples = " << numberOfSamples << std::endl;
    std::cout << "sampleStart = " << sampleStart << std::endl;

    ALSVINN_LOG(INFO, "sampleStart = " << sampleStart);

    std::vector<size_t> samples;
    samples.reserve(numberOfSamples);

    for (size_t i = 0; i < numberOfSamples; ++i) {
        samples.push_back(sampleStart + i);
    }


    mpi::SimpleLoadBalancer loadBalancer(samples);

    auto loadBalanceConfiguration = loadBalancer.loadBalance(multiSample,
            multiSpatial,
            *mpiConfigurationWorld);
    auto& samplesForProc = std::get<0>(loadBalanceConfiguration);
    auto statisticalConfiguration = std::get<1>(loadBalanceConfiguration);
    auto spatialConfiguration = std::get<2>(loadBalanceConfiguration);


    auto simulatorCreator = std::dynamic_pointer_cast<run::SimulatorCreator>
        (std::make_shared<run::FiniteVolumeSimulatorCreator>
            (inputFilename,
	     numberOfSamples,
                spatialConfiguration,
                statisticalConfiguration,
                mpiConfigurationWorld,
                multiSpatial));

    auto name = boost::algorithm::trim_copy(
            configuration.get<std::string>("fvm.name"));
    auto runner = std::make_shared<run::Runner>(simulatorCreator, sampleGenerator,
            samplesForProc,
            statisticalConfiguration, name);
    const auto grid = createGrid(configuration);
    auto statistics  = createStatistics(*grid, numberOfSamples, configuration, statisticalConfiguration,
            spatialConfiguration, mpiConfigurationWorld);
    runner->setStatistics(statistics);

    // We want to make sure everything is created before going further
    MPI_Barrier(mpiConfigurationWorld->getCommunicator());
    return runner;
}

std::shared_ptr<samples::SampleGenerator> Setup::makeSampleGenerator(
    const std::string& inputFilename) {
    auto& textCache = alsutils::io::TextFileCache::getInstance();
    auto configurationFileContent = textCache.loadTextFile(inputFilename);
    std::stringstream stream(configurationFileContent);
    ptree configurationBase;
    boost::property_tree::read_xml(stream, configurationBase);
    auto configuration = configurationBase.get_child("config");
    auto sampleGenerator = makeSampleGenerator(configuration);

    return sampleGenerator;

}

size_t Setup::readNumberOfSamples(const std::string& inputFilename) {
    auto& textCache = alsutils::io::TextFileCache::getInstance();
    auto configurationFileContent = textCache.loadTextFile(inputFilename);
    std::stringstream stream(configurationFileContent);
    ptree configurationBase;
    boost::property_tree::read_xml(stream, configurationBase);
    auto configuration = configurationBase.get_child("config");

    return readNumberOfSamples(configuration);
}

size_t Setup::readSampleStart(const std::string& inputFilename) {
    auto& textCache = alsutils::io::TextFileCache::getInstance();
    auto configurationFileContent = textCache.loadTextFile(inputFilename);
    std::stringstream stream(configurationFileContent);
    ptree configurationBase;
    boost::property_tree::read_xml(stream, configurationBase);
    auto configuration = configurationBase.get_child("config");

    return readSampleStart(configuration);
}

std::shared_ptr<samples::SampleGenerator> Setup::makeSampleGenerator(
    Setup::ptree& configuration) {
    auto numberOfSamples = readNumberOfSamples(configuration);

    samples::SampleGenerator::GeneratorDistributionMap generators;

    auto generatorName = configuration.get<std::string>("uq.generator");

    if (generatorName == "auto") {
        generatorName = "stlmersenne";
    }

    auto parametersNode = configuration.get_child("uq.parameters");
    generator::GeneratorFactory generatorFactory;
    distribution::DistributionFactory distributionFactory;

    for (auto parameterNode : parametersNode) {
        auto name = parameterNode.second.get<std::string>("name");
        auto length = parameterNode.second.get<size_t>("length");
        auto type = parameterNode.second.get<std::string>("type");

        distribution::Parameters parametersToDistribution(parameterNode.second);
        parametersToDistribution.addDoubleParameter("lower", 0);
        parametersToDistribution.addDoubleParameter("upper", 1);

        parametersToDistribution.addDoubleParameter("a", 0);
        parametersToDistribution.addDoubleParameter("b", 1);


        parametersToDistribution.addDoubleParameter("mean", 0);
        parametersToDistribution.addDoubleParameter("sd", 1);

        auto distribution = distributionFactory.createDistribution(type,
                length,
                numberOfSamples,
                parametersToDistribution);


        auto generator = generatorFactory.makeGenerator(generatorName, length,
                numberOfSamples);
        generators[name] = std::make_pair(length, std::make_pair(generator,
                    distribution));
    }



    return std::make_shared<samples::SampleGenerator> (generators);


}


// <stats>
//
//   <stat>
//   <name>
//   structure
//   </name>
//   <numberOfSaves>1</numberOfSaves>
//   <writer>
//
//     <type>netcdf</type>
//     <basename>kh_structure</basename>
//   </writer>
//   </stat>
//   <stat>
//   <name>
//     meanvar
//   </name>
//   <numberOfSaves>1</numberOfSaves>
//   <writer>
//     <type>netcdf</type>
//     <basename>kh_structure</basename>
//   </writer>
//   </stat>
// </stats>
std::vector<std::shared_ptr<stats::Statistics> > Setup::createStatistics(
    const alsfvm::grid::Grid &grid,
    size_t num_samples,
    Setup::ptree& configuration,
    mpi::ConfigurationPtr statisticalConfiguration,
    mpi::ConfigurationPtr spatialConfiguration,
    mpi::ConfigurationPtr worldConfiguration) {
    auto statisticsNodes = configuration.get_child("uq.stats");
    stats::StatisticsFactory statisticsFactory;
    std::shared_ptr<alsfvm::io::WriterFactory> writerFactory;


    writerFactory.reset(new alsfvm::io::MpiWriterFactory(spatialConfiguration));

    auto platform = configuration.get<std::string>("fvm.platform");
    std::vector<std::shared_ptr<stats::Statistics> > statisticsVector;

    for (auto& statisticsNode : statisticsNodes) {
        auto name = statisticsNode.second.get<std::string>("name");
        boost::trim(name);
        stats::StatisticsParameters parameters(statisticsNode.second);
        parameters.setMpiConfiguration(statisticalConfiguration);
        parameters.setNumberOfSamples(readNumberOfSamples(configuration));
        parameters.setPlatform(platform);
        auto statistics = statisticsFactory.makeStatistics(platform, name, parameters);




        // Make writer
        std::string type = statisticsNode.second.get<std::string>("writer.type");
        std::string basename =
            statisticsNode.second.get<std::string>("writer.basename");

        for (auto statisticsName : statistics->getStatisticsNames()) {

            auto outputname = basename + "_" + statisticsName;
	    const auto& writerNode = statisticsNode.second.get_child("writer");
	    size_t numberOfSaves = writerNode.get<size_t>("numberOfSaves");
	    bool writeInitialTimestep = false;
	    if (writerNode.find("writeInitialTimestep") != writerNode.not_found()) {
		writeInitialTimestep = writerNode.get<bool>("writeInitialTimestep");
	    }
	    std::vector<real> timesteps(numberOfSaves + writeInitialTimestep);
	    timesteps[0] = 0;
	    for (size_t t = 0 ; t < numberOfSaves ; ++t) {
	      timesteps[t + writeInitialTimestep] = static_cast<real>(t + 1) / numberOfSaves;
	    }
            auto baseWriter = writerFactory->createWriter(type, outputname, grid, num_samples, timesteps,
                    alsfvm::io::Parameters(statisticsNode.second.get_child("writer")));
            baseWriter->addAttributes("uqAttributes", configuration);
            statistics->addWriter(statisticsName, baseWriter);
        }

        if (statisticsNode.second.find("numberOfSaves") !=
            statisticsNode.second.not_found()) {

            bool writeInitialTimestep = false;

            if (statisticsNode.second.find("writeInitialTimestep") !=
                statisticsNode.second.not_found()) {
                writeInitialTimestep = statisticsNode.second.get<bool>("writeInitialTimestep");
            }


            auto numberOfSaves = statisticsNode.second.get<size_t>("numberOfSaves");
            ALSVINN_LOG(INFO, "statistics.numberOfSaves = " << numberOfSaves);
            real endTime = configuration.get<real>("fvm.endTime");
            real timeInterval = endTime / numberOfSaves;
            auto statisticsInterval =
                std::shared_ptr<stats::Statistics>(
                    new stats::FixedIntervalStatistics(statistics, timeInterval,
                        endTime, writeInitialTimestep));
            statisticsVector.push_back(statisticsInterval);
        } else if (statisticsNode.second.find("time") !=
            statisticsNode.second.not_found()) {
            const real time = statisticsNode.second.get<real>("time");
            const real radius = statisticsNode.second.get<real>("timeRadius");

            auto timeIntegrator = std::shared_ptr<stats::Statistics>(
                    new stats::TimeIntegratedWriter(statistics, time,
                        radius));
            statisticsVector.push_back(timeIntegrator);
        }


    }

    return statisticsVector;
}

size_t Setup::readNumberOfSamples(Setup::ptree& configuration) {
    return configuration.get<real>("uq.samples");
}

size_t Setup::readSampleStart(Setup::ptree& configuration) {
    auto& uq = configuration.get_child("uq");

    if (uq.find("sampleStart") != uq.not_found()) {
        ALSVINN_LOG(INFO, "sampleStart tag present");
        return uq.get<size_t>("sampleStart");
    }

    return 0;
}

alsfvm::shared_ptr<alsfvm::grid::Grid> Setup::createGrid(
    const Setup::ptree& configuration) {
    const ptree& gridNode =  configuration.get_child("fvm.grid");

    const std::string& lowerCornerString = gridNode.get<std::string>("lowerCorner");
    const std::string& upperCornerString = gridNode.get<std::string>("upperCorner");
    const std::string& dimensionString = gridNode.get<std::string>("dimension");


    auto lowerCorner = parseVector<real>(lowerCornerString);
    auto upperCorner = parseVector<real>(upperCornerString);
    auto dimension = parseVector<int>(dimensionString);

    std::array<alsfvm::boundary::Type, 6> boundaryConditions;

    return alsfvm::make_shared<alsfvm::grid::Grid>(lowerCorner, upperCorner, dimension,
            boundaryConditions);
}
}
}
