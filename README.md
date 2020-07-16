## mpcc

This repo contains an inprogress implementation of model predictive contouring control based on the ACADOS solver.

#setup
1. Clone this repo

2. Head to https://github.com/acados/acados and install acados following their instructions.
   Below is a brief summary.
   
2.2 git clone git@github.com:acados/acados.git

2.3 cd acados

2.4 git submodule update --recursive --init

2.5 mkdir -p build

    cd build
    
    cmake ..
    
    make install
    
2.6 install python dependencies:
    pip3 install casadi
    
2.7 make python bindings:

    cd <acados_root>/build
    
    cmake -DACADOS_WITH_QPOASES=ON ..
    
    make install -j4
    
    pip3 install <acados_root>/interfaces/acados_template
    
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"<acados_root>/lib"
    
    export ACADOS_SOURCE_DIR=<acados_root>
    
    Note: I ended up putting the last two commands into my .bashrc
    
    It should look something like this:
    
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/home/pw/acados/lib/"
    
    export ACADOS_SOURCE_DIR=/home/pw/acados
    
2.8 test installation:

    cd <acados_root>/examples/acados_python/getting_started
    
    python3 minimal_example_closed_loop.py    
    
3.  All set up!
    Currently i am working in mpcc/scripts/acadosmodel
    try running "python3 python_sim.py"
