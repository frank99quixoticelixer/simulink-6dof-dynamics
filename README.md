## Model Architecture Overview

This project implements a custom **six-degree-of-freedom (6DOF) rigid-body dynamics model** using Simulink and user-defined MATLAB functions. The architecture follows a standard separation between **dynamics** and **kinematics**, with explicit state integration and reference-frame transformations.

The model is structured into three main subsystems:

- **Translational Dynamics** (linear motion expressed in the body frame)
- **Rotational Dynamics** (angular motion expressed in the body frame)
- **Kinematics** (orientation propagation and frame transformations)

All state derivatives are integrated explicitly using Simulink integrator blocks, allowing full visibility of the underlying equations and data flow.

---

## Translational Dynamics
This function computes the time derivative of the body-frame linear velocity vector.

The translational dynamics are expressed in the body reference frame and include the Coriolis term resulting from angular motion:

$$
\dot{\mathbf{V}}_b = \frac{1}{m}\mathbf{F} - \boldsymbol{\omega}_b \times \mathbf{V}_b
$$


Where:

- F is the applied force vector in the body frame
- m is the vehicle mass
- V_b is the body-frame linear velocity
- ω_b is the body-frame angular velocity

A rotation matrix from the body frame to the inertial frame is also computed using Euler angles. This matrix is later used for kinematic propagation and for transforming vectors between the body and inertial reference frames.

```matlab
 function V_dot = Traslational_Dynamics(F, V_b, w_b, m, e)

    % Force inputs into a column vector
    F   = F(:);
    V_b = V_b(:);
    w_b = w_b(:);

    % Coriolis effect 
    cross_term = cross(w_b, V_b);


    
    phi = e(1);
    theta = e(2);
    psi = e(3);

    % --- Rotation Matrix ---
    cpsi = cos(psi); spsi = sin(psi);
    cth = cos(theta); sth = sin(theta);
    cphi = cos(phi); sphi = sin(phi);

    R_eb = [ ...
        cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
        cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
        -sth,       sphi*cth,                 cphi*cth ];

    Rot = R_eb

    % Traslational Dynamics
    V_dot = (F / m) - cross_term;

  
end
```
## Rotational Dynamics

The rotational dynamics of the rigid body expressed in the body frame are governed by Euler’s rotational equations:

$$
\dot{\boldsymbol{\omega}}_b = \mathbf{I}^{-1}
\left(
\mathbf{M} - \boldsymbol{\omega}_b \times (\mathbf{I}\boldsymbol{\omega}_b)
\right)
$$

Where:

- **M** is the applied moment vector
- **I** is the inertia matrix
- **ω_b** is the body-frame angular rate vector

The cross-product term accounts for gyroscopic coupling between the rotational axes.
```matlab
function w_dot = Rotational_Dynamics(M, w_b, I)

    % Force inputs into a column vector
    M   = M(:);
    w_b = w_b(:);

    % Dot procut for the Inertia term
    Iw_b = I * w_b;

    % Cross product
    cross_term = cross(w_b, Iw_b);

    % Rotational Dynamics op
    w_dot = I \ (M - cross_term);
end
```
## Kinematics

The kinematics block is responsible for:

- Euler angle propagation  
- Transformation of linear velocity from the body frame to the inertial frame  
- Rotation matrix computation  

Euler angle rates are computed using the standard transformation matrix:

$$
\dot{\mathbf{e}} = \mathbf{J}(\mathbf{e}) \, \boldsymbol{\omega}_b
$$

The inertial-frame linear velocity is obtained as:

$$
\mathbf{V}_e = \mathbf{R}_{eb} \, \mathbf{V}_b
$$

Where:

- **e** is the vector of Euler angles  
- **J(e)** is the Euler angle rate transformation matrix  
- **ω_b** is the body-frame angular rate vector  
- **V_b** is the body-frame linear velocity  
- **R_{eb}** is the rotation matrix from the body frame to the inertial frame, computed directly from the Euler angles
```matlab
function [e_dot, V_e, Rot] = Kinematics(w_b, V_b, e)

    % Force inputs into a column vectora
    w_b = w_b(:);
    V_b = V_b(:);
    e   = e(:);

    phi = e(1);
    theta = e(2);
    psi = e(3);

    % --- J Matrix ---
    J_e = [ ...
        1, sin(phi)*tan(theta),  cos(phi)*tan(theta);
        0,      cos(phi),       -sin(phi);
        0, sin(phi)/cos(theta),  cos(phi)/cos(theta) ];

    % Euler Derivation
    e_dot = J_e * w_b;   % → 3x1

    % --- Rotation Matrix ---
    cpsi = cos(psi); spsi = sin(psi);
    cth = cos(theta); sth = sin(theta);
    cphi = cos(phi); sphi = sin(phi);

    R_eb = [ ...
        cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
        cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
        -sth,       sphi*cth,                 cphi*cth ];

    Rot = R_eb


    
    % Inertial rates
    V_e = R_eb * V_b;   % → 3x1

end
```
## Euler Angle Wrapping

To avoid numerical drift and discontinuities, Euler angles are wrapped to the interval \([-\pi, \pi]\) using the identity:

$$
\mathbf{e} = \arctan\!\left(\frac{\sin(\mathbf{e})}{\cos(\mathbf{e})}\right)
$$


This operation is applied after state integration to ensure that the orientation states remain bounded and well-defined during long-duration simulations.

