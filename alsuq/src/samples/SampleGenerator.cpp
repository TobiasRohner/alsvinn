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

#include "alsuq/samples/SampleGenerator.hpp"
#include "alsutils/error/Exception.hpp"

namespace alsuq {
namespace samples {

SampleGenerator::SampleGenerator(const
    SampleGenerator::GeneratorDistributionMap& generators)
    : generators(generators) {

}

std::vector<real> SampleGenerator::generate(const std::string& parameter,
    const size_t sampleIndex) {
    if (generators.find(parameter) == generators.end()) {
        THROW("Unknown parameter " << parameter);
    }

    if (sampleIndex < currentSample) {
        THROW("This shouldn't happen, we have requested a sample number lower"
            << " than what we have generated. "
            << std::endl << "sampleIndex = " << sampleIndex << std::endl
            << "currentSample = " << currentSample);
    }

    size_t dimension = generators[parameter].first;
    auto generator = generators[parameter].second.first;
    auto distribution =  generators[parameter].second.second;

    // now we throw away samples we do not need
    while (currentSample < sampleIndex) {
        for (size_t i = 0; i < dimension; ++i) {
            distribution->generate(*generator, i);
        }

        currentSample++;
    }





    std::vector<real> samples(dimension);

    for (size_t i = 0; i < dimension; ++i) {
        samples[i] = distribution->generate(*generator, i);
    }

    currentSample++;
    return samples;
}

std::vector<std::string> SampleGenerator::getParameterList() const {
    std::vector<std::string> parameters;

    for (auto pair : generators) {
        parameters.push_back(pair.first);
    }

    return parameters;
}

}
}
