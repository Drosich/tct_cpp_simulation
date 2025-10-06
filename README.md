# TCT pulse simulation
This is an attempt at making a simulation of Transient Photon Absorption 
Transient Current Technique (TPA-TCT) pulses from scratch.
It is intended to provide a somewhat accurate description of the pulses we
observe at the measurement campaigns.
This simply simulates a 2D charge injection into a semiconductor corresponding
to a laser providing a gaussian beam. The charges generated are then transported
with a simple 2D transport algorithm along the electric field generated inside
a PN junction diode. Diffusion will be also simulated.

# Description
...

# How to run
The code is compiled with CMake. You can use the CmakeLists file provided here
as usual

```bash
$ mkdir build
$ cmake ..
$ make
```

An single executable file will then be created

```bash
$ ./tct_sim
```

# Dependencies
This code uses CERN's ROOT framework for the visualization of the results. In orther for this to work ROOT has to be installed

# Authors
Diego Rosich (IFCA-CSIC-UC)