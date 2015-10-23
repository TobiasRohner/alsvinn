#include "alsfvm/reconstruction/ReconstructionCUDA.hpp"
#include "alsfvm/reconstruction/WENOF2.hpp"
#include "alsfvm/reconstruction/WENO2.hpp"
#include "alsfvm/equation/euler/Euler.hpp"



namespace alsfvm {
    namespace reconstruction {
        namespace {


            template<class Equation, class ReconstructionType,  size_t dimension, bool xDir, bool yDir, bool zDir>
            __global__ void reconstructDevice(Equation eq, typename Equation::ConstViews input,
                typename Equation::Views left, typename Equation::Views right,
                
                size_t numberOfXCells, size_t numberOfYCells, size_t numberOfZCells) {

                const size_t index = threadIdx.x + blockDim.x * blockIdx.x;
                // We have
                // index = z * nx * ny + y * nx + x;
                const size_t xInternalFormat = index % numberOfXCells;
                const size_t yInternalFormat = (index / numberOfXCells) % numberOfYCells;
                const size_t zInternalFormat = (index) / (numberOfXCells * numberOfYCells);

                if (xInternalFormat >= numberOfXCells || yInternalFormat >= numberOfYCells || zInternalFormat >= numberOfZCells) {
                    return;
                }
                const size_t x = xInternalFormat + (1) * xDir;
                const size_t y = yInternalFormat + (1) * yDir;
                const size_t z = zInternalFormat + (1) * zDir;

                ReconstructionType::reconstruct(eq, input, x, y, z, left, right, xDir, yDir, zDir);


            }

            template<class Equation, class ReconstructionType, size_t dimension, bool xDir, bool yDir, bool zDir >
            void callReconstructionDevice(const Equation& equation, const volume::Volume& inputVariables,
                size_t direction,
                size_t indicatorVariable,
                volume::Volume& leftOut,
                volume::Volume& rightOut) {

                const size_t numberOfXCells = leftOut.getTotalNumberOfXCells() - 2 * (1) * xDir;
                const size_t numberOfYCells = leftOut.getTotalNumberOfYCells() - 2 * (1) * yDir;
                const size_t numberOfZCells = leftOut.getTotalNumberOfZCells() - 2 * (1) * zDir;

                const size_t totalSize = numberOfXCells * numberOfYCells * numberOfZCells;


                const size_t blockSize = 512;
                const size_t gridSize = (totalSize + blockSize - 1) / blockSize;



                typename Equation::Views viewLeft(leftOut);
                typename Equation::Views viewRight(rightOut);
                typename Equation::ConstViews viewInput(inputVariables);

                reconstructDevice<Equation, ReconstructionType, dimension, xDir, yDir, zDir> << <gridSize, blockSize >> >(equation, viewInput,
                    viewLeft, viewRight,
                    numberOfXCells, numberOfYCells, numberOfZCells);

            }

            template<size_t dimension, class Equation, class ReconstructionType>
            void performReconstructionDevice(const Equation& equation, const volume::Volume& inputVariables,
                size_t direction,
                size_t indicatorVariable,
                volume::Volume& leftOut,
                volume::Volume& rightOut) {
                assert(direction < 3);
                switch (direction) {
                case 0:
                    callReconstructionDevice<Equation, ReconstructionType, dimension, 1, 0, 0>(equation, inputVariables, direction, indicatorVariable, leftOut, rightOut);
                    break;

                case 1:
                    callReconstructionDevice<Equation, ReconstructionType, dimension, 0, 1, 0>(equation, inputVariables, direction, indicatorVariable, leftOut, rightOut);
                    break;

                case 2:
                    callReconstructionDevice<Equation, ReconstructionType, dimension, 0, 0, 1>(equation, inputVariables, direction, indicatorVariable, leftOut, rightOut);
                    break;
                }

            }

        }



        ///
        /// Performs reconstruction.
        /// \param[in] inputVariables the variables to reconstruct.
        /// \param[in] direction the direction:
        /// direction | explanation
        /// ----------|------------
        ///     0     |   x-direction
        ///     1     |   y-direction
        ///     2     |   z-direction
        ///
        /// \param[in] indicatorVariable the variable number to use for
        /// stencil selection. We will determine the stencil based on
        /// inputVariables->getScalarMemoryArea(indicatorVariable).
        ///
        /// \param[out] leftOut at the end, will contain the left interpolated values
        ///                     for all grid cells in the interior.
        ///
        /// \param[out] rightOut at the end, will contain the right interpolated values
        ///                     for all grid cells in the interior.
        ///
        template<class ReconstructionType, class Equation>
        void ReconstructionCUDA<ReconstructionType, Equation>::performReconstruction(const volume::Volume& inputVariables,
            size_t direction,
            size_t indicatorVariable,
            volume::Volume& leftOut,
            volume::Volume& rightOut)
        {
            size_t dimension = 1 + (leftOut.getNumberOfYCells() > 1) + (leftOut.getNumberOfZCells() > 1);

            switch (dimension) {
            case 1:
                performReconstructionDevice<1, Equation, ReconstructionType>(equation, inputVariables, direction, indicatorVariable, leftOut, rightOut);
                break;
            case 2:
                performReconstructionDevice<2, Equation, ReconstructionType>(equation, inputVariables, direction, indicatorVariable, leftOut, rightOut);
                break;
            case 3:
                performReconstructionDevice<3, Equation, ReconstructionType>(equation, inputVariables, direction, indicatorVariable, leftOut, rightOut);
                break;
            }

        }

        template<class ReconstructionType, class Equation>
        size_t ReconstructionCUDA<ReconstructionType, Equation>::getNumberOfGhostCells()  {
            return ReconstructionType::getNumberOfGhostCells();
        }

        template<class ReconstructionType, class Equation>
        ReconstructionCUDA<ReconstructionType, Equation>::ReconstructionCUDA(simulator::SimulatorParameters& parameters)
            : equation(static_cast<typename Equation::Parameters&>(parameters.getEquationParameters()))
        {


        }

        template class ReconstructionCUDA < WENOF2<equation::euler::Euler>, equation::euler::Euler>;
        template class ReconstructionCUDA < WENO2<equation::euler::Euler>, equation::euler::Euler>;

    }
}
