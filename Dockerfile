FROM nvidia/cuda:12.6.3-devel-ubuntu24.04

WORKDIR alsvinn

RUN apt-get update && \
    apt-get install --no-install-recommends -y git cmake wget libssl-dev libhdf5-mpi-dev python3-dev libfftw3-dev python3-numpy libboost-all-dev libnetcdf-mpi-dev libgtest-dev libpnetcdf-dev && \
    apt-get clean

#RUN wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.0.tar.gz \
#    && tar -xvf openmpi-5.0.0.tar.gz \
#    && cd openmpi-5.0.0 \
#    && ./configure --with-pmix \
#    && make -j$(nproc) \
#    && make install

COPY addons addons
COPY alsfvm alsfvm
COPY alsuq alsuq
COPY alsuqcli alsuqcli
COPY alsutils alsutils
COPY alsvinncli alsvinncli
COPY assets assets
COPY cmake cmake
COPY containers containers
COPY documentation documentation
COPY examples examples
COPY library_examples library_examples
COPY python python
COPY test test
COPY utils utils
COPY CMakeLists.txt .

RUN cmake -DCMAKE_BUILD_TYPE=Release -DALSVINN_USE_FLOAT=ON -DALSVINN_USE_CUDA=ON -DALSVINN_BUILD_DOXYGEN=OFF -DALSVINN_PYTHON_VERSION=3.12 -DNETCDF_DIR=/usr/lib/aarch64-linux-gnu/netcdf/mpi -DPNETCDF_DIR=/usr/lib/aarch64-linux-gnu/netcdf/pnetcdf -S . -B build && \
    cmake --build build -j$(nproc)

ENTRYPOINT ["/alsvinn/build/alsuqcli/alsuqcli"]
