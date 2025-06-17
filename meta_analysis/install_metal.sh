cd /well/lindgren/dpalmer
git clone git@github.com:statgen/METAL.git

# # Local
# # Ensure homebrew is installed
# brew install cmake

# BMRC
module load CMake

cd METAL
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make test
make install
