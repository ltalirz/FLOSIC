from spack import *

class Flosic23(MakefilePackage):
    """FLOSIC is a software package for performing Fermi-Löwdin orbital self-interaction corrections."""

    homepage = "https://github.com/FLOSIC/PublicRelease_2023"
    url      = "https://github.com/FLOSIC/PublicRelease_2023/archive/refs/tags/v1.0.tar.gz"
    git      = "https://github.com/FLOSIC/PublicRelease_2023.git"

    version('1.0', sha256='830ef06669b6f0e14b77c579dae64c4abdf399a57590c13bdabd1071ea984f48')

    variant('mpi', default=True, description='Enable MPI support')

    #depends_on('openblas')
    depends_on('mpi', when='+mpi')

    def edit(self, spec, prefix):
       # Add linker flags to the Makefile before compiling
        makefile_path = join_path('flosic', 'Makefile.mpi')
        #openblas_flags = self.spec['openblas'].libs
        #filter_file('-llapack -lblas', ' '.join(openblas_flags), makefile_path)

        makefile = FileFilter('flosic/Makefile.mpi')

        # Adjust include and library directories
        include_path = spec['mpi'].prefix.include
        #lib_path = spec['openblas'].prefix.lib

        # Adjust the Makefile to include the correct paths for MPI and OpenBLAS
        #makefile.filter('^IDIR = .*', 'IDIR = -I' + include_path)
        #makefile.filter('^LDIR = .*', 'LDIR = -L' + lib_path)

        # Link against the OpenBLAS library
        #openblas_flags = self.spec['openblas'].libs
        #makefile.filter('LIBS = -llapack -lblas', 'LIBS = {}'.format(openblas_flags))

        with working_dir('flosic'):
            makefile = FileFilter('Makefile.mpi')

            # Set MPI option
            if '+mpi' in spec:
                makefile.filter('-DMPI', '-DMPI')
            else:
                makefile.filter('-DMPI', ' ' )

            # Add any additional filters or changes needed here

            # Adjust compilers and flags if needed
            # Example:
            # makefile.filter('CC = gcc', 'CC = {}'.format(spack_cc))
            # makefile.filter('FC = mpif90', 'FC = {}'.format(spack_fc))
            makefile.filter('FFLAGS = -O3 -mcmodel=large', 'FFLAGS = -O3 -mcmodel=large -fallow-argument-mismatch')

            # Rename Makefile.fedora to Makefile
            move('Makefile.mpi', 'Makefile')
            mkdirp('bin')


    def build(self, spec, prefix):
        with working_dir('flosic'):
            fc = Executable(spack_fc)

            # Compile modules that generate .mod files first
            fc('-o', 'condcomp', 'condcomp.f')

            make(parallel=False)
 

    def install(self, spec, prefix):
        # Continue with the rest of the installation steps
        mkdirp(prefix.bin)
        install('flosic/bin/mpnrlmol.mpi', prefix.bin)

        # Handle other installation steps as needed
 
