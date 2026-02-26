#!/bin/bash

git clone https://github.com/cosmolattice/cosmolattice.git CL1.0
cp -r ./analogues_CL1.0/* CL1.0/
mkdir CL1.0/build
cd CL1.0/build
cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DMPI=ON
make -j bench-phi4
