# TCT pulse simulation
This is an attempt at making a simulation of Transient Photon Absorption 
Transient Current Technique (TPA-TCT) pulses from scratch.
It is intended to provide a somewhat accurate description of the pulses we
observe at the measurement campaigns.
This simply simulates a 2D charge injection into a semiconductor corresponding
to a laser providing a gaussian beam. The charges generated are then transported
with a simple 2D transport algorithm along the electric field generated inside
a PN junction diode. Diffusion will be also simulated.

