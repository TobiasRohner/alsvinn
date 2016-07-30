#pragma once
#include "alsfvm/diffusion/DiffusionOperator.hpp"
#include "alsfvm/volume/VolumeFactory.hpp"
#include "alsfvm/simulator/SimulatorParameters.hpp"
#include "alsfvm/reconstruction/ReconstructionFactory.hpp"

namespace alsfvm { namespace diffusion { 


    /// Applies the Tecno diffusion to the operator. This will always take 
    /// the form
    ///
    /// \f[R\Lambda R^{-1} \langle\langle v\rangle \rangle\f]
    ///
    /// where \f$R\f$ is the matrix of eigenvalues of the flux jacobian, and 
    /// $f\Lambda\f$ is either the Rusanov or Roe matrix. See
    /// http://www.cscamm.umd.edu/tadmor/pub/TV+entropy/Fjordholm_Mishra_Tadmor_SINUM2012.pdf
    /// 
    /// The matrix \f$\Lambda\f$ is specified through the DiffusionMatrix template argument.
    template<class Equation, class DiffusionMatrix>
    class TecnoDiffusionCPU : public DiffusionOperator {
    public:

        TecnoDiffusionCPU(volume::VolumeFactory& volumeFactory,
            alsfvm::shared_ptr<reconstruction::Reconstruction>& reconstructionFactory,
            const alsfvm::shared_ptr<simulator::SimulatorParameters>& simulatorParameters,
            size_t nx, size_t ny, size_t nz);

        /// Applies numerical diffusion to the outputVolume given the data in conservedVolume.
        ///
        /// \note The numerical diffusion will be added to outputVolume, ie. the code will 
        /// essentially work like
        /// \code{.cpp}
        /// outputVolume += diffusion(conservedVolume);
        /// \endcode
        virtual void applyDiffusion(volume::Volume& outputVolume,
            const volume::Volume& conservedVolume);

    private:
        alsfvm::shared_ptr<reconstruction::Reconstruction> reconstruction;

        // Reconstructed values (these are basically R^{-1}v, 
        // where v is the entropy variables and R^{-1} is the inverse of the 
        // eigenvalues of the flux.
        alsfvm::shared_ptr<volume::Volume> left;
        alsfvm::shared_ptr<volume::Volume> right;

        alsfvm::shared_ptr<volume::Volume> entropyVariables;

        Equation equation;
    };
} // namespace diffusion
} // namespace alsfvm