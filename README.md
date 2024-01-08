# TO_integration_scheme

Code for my thesis.

I compared the time-step-size-consistency of various integrators in Trajectory Optimization (TO) problems.
Integrators of superior time-step-size-consistency are exepcted to perform well in sim-to-real transfer.
Experiments were conducted on the cartpole swing up problem, and the underactuated double pendulum problem.
The integrators of interest are the symplectic Euler integrator(SEI), the variational integrator(VI), the passive midpoint integrator(PMI), and the fourth-order Runge-Kutta method (RK4).

Results show that RK4 consistently outperforms the others, followed by PMI, VI and SEI respectively.
