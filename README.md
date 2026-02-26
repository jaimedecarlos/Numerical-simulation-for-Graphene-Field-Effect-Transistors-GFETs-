# Numerical-simulation-for-Graphene-Field-Effect-Transistors-GFETs-

Numerical simulation in MATLAB for a physical model of Graphene FETs to study static saturation points. Designed through the Finite Differences Method. 

Possibility to extend the model to transient models for frequency analysis. 

## Installation

Developed in MATLAB R2021a - academic use. 

No specific actions needed

## Usage

TFG__tarea3_v3.m has the simulation code. 

Inside directory "results/": Graphics.m has the graphics code for plotting specific voltage input. Simulation data can be extracted from variables.mat

### Instructions to plot specific static curves from variables.dat

The constructed data only accepts drain voltage from 0V to 1V with 0.05V intervals, and gate voltage from -2V to 2V with 0.5V intervals. 
- To fix $V_D$ (drain voltage), fix $i$ int variable in Graphics.m, from 1 to 21. 
- To fix $V_G$ (gate voltage), fix $j$ int variable in Graphics.m, from 1 to 9. 

To plot $V_D=1V$, $V_G=1V$, state i=21; j=7
To plot $V_D=0.5V$, $V_G=-2V$, state i=11; j=1
etc

<img src="https://github.com/user-attachments/assets/ce659696-8b8e-4db9-80dd-d90cfa6c4d1f" width="300" />
<img src="https://github.com/user-attachments/assets/0004ca6b-1763-490c-a1b9-98a0e9dca558" width="300" />
<img src="https://github.com/user-attachments/assets/c1ad1318-c1cf-40bc-8c62-23898d4d2d01" width="300" />
<img src="https://github.com/user-attachments/assets/a9416eb9-c28c-4463-83ec-47e549cd57d7" width="300" />

## Analytical work

Mathematical proof and construction explained in PDF file. 

This work effectively couples the dominant equations for the system: 
- Electrostatic poisson equation in 2D.
- Solid physics equation: current conservation equation derived from semiconductor physics principles. The result is an ODE in the y-axis. 

### Poisson equation 2D

Solved numerically through a system of linear equations, derived through the 2nd order Taylor polynomial. 

### Current conservation equation

Given to its heavily nonlinear form, this model was solved through the shooting method, via a 2nd order RK method. 
