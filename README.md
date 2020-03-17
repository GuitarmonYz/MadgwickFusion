# MadgwickFusion
This repo contains implementation of both the improved Madgwick algorithm proposed in 'Formulation of a new gradient descent MARG orientation algorithm: Case study on robot teleoperation' https://doi.org/10.1016/j.ymssp.2019.04.064 and the algorithm in the original paper. The only difference between the two is ref vector for magnetic field, hence the equations for gradient calcuation are different.

Unlike the heuristic approach in the paper where the step size $\mu$, which adaptly changing based on the quaternion changing rate, is canceled out. In this implementation, the step size is preserved and the step size amplifier $\alpha$ needs to be set by the user. This is to achieve more explicit behavior as stated in the paper: 

> The combination of sensors allows the gyro to track orientation during high freq motion while gyro drift is compensated during low frequency using gradient descent steps.

The code here will not compile by just downloading it. It's up to the user to implement define sensor data struct, and calculate $\Delta t$, since they are platform dependant. 

Hope this helps people who are interested in motion tracking. Some extention tweaks are provided, one can use transient shaper to filter on quaternion changing rate, and use extension for gradient decent defined in gradient_optimize.h.