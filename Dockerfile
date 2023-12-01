# Use an official Ubuntu base image with OpenMPI
FROM ubuntu:latest

# Install necessary packages
RUN apt-get update && apt-get install -y \
    openmpi-bin \
    automake \
    gfortran \
    wget \
    make \
    && rm -rf /var/lib/apt/lists/*

COPY flosic /usr/src/flosic

# Set the working directory
WORKDIR /usr/src/flosic

# Compile FLOSIC
RUN gfortran -o condcomp condcomp.f
RUN cp Makefile.mpi Makefile
RUN make

# Install FLOSIC
RUN mkdir -p /usr/local/bin \
    && cp bin/mpnrlmol.mpi /usr/local/bin/

# Set the default command or entrypoint, if necessary
CMD ["bash"]
