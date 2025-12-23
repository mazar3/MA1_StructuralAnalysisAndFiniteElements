# Past exams questions

---

## Chapter 3-4 - FEM basics

### Multiple-choice questions
- In structural finite elements, what are the units of an element $k_{ij}$ in the stiffness matrix $K$? ($N/m^2$ / $N/m$ / $kg/m$ / $N$)
- The stress equilibrium at inter-elemental boundaries is always satisfied in FEM. (True / False)
- Isotropic linear elastic materials can be fully described by how many material parameters? (2 / 3)
- Gauss divergence theorem is expressed as: ($\int_V u_i n_i dV = \oint_S u_{i,i} dS$ / $\int_V u_{i,j} dV = \oint_S u_i n_j dS$)
- The shear modulus of elasticity, $G$, is always positive. (True / False)
- The units of the strain energy density $\int \tau_{ij} d\epsilon_{ij}$ are: ($N/m^3$ / $J/m^3$ / $J/m^2$ / $W/m^2$)
- The units of the body forces $b_i$ are: ($N$ / $kg/m^3$ / $N/m^3$ / $kg/m$)
- The variational formulation states that at the equilibrium, the total potential energy $\Pi$ is equal to 0. (True / False)
- A displacement field $u$ is said to be *kinematically admissible* if it vanishes on the surface $S_u$, i.e. $u = 0$ on $S_u$. (True / False)
- In linear elasticity, the strain energy of the theoretical solution is a lower bound of the finite element solutions. (True / False)
- The weak formulation of the linear elastic problem involves the second derivatives of the displacements. (True / False)

### Open-ended questions
- Starting from the weak form of the equilibrium equations:
  $$ \int_V \sigma_{ij}(u) w_{(i, j)}dV = \int_V f_i w_i dV + \int_{S_T} w_i T^{(n)}_i dS $$
  a) Derive the *discrete equilibrium equation* for a homogeneous, linear elastic isotropic material using the Galerkin method to explain what *discretization error in the finite element method is*.
  b) Explain the concept of *h-refinement* and how you use it in practice to achieve a *converged mesh*.
  All of the variables should be defined, and the main steps of the mathematical developments need to be explained. List all of the sources of errors and the assumptions introduced in the presented successive steps.
- What are the sources of errors in the Finite Element Method? List all of them. Perform the developments of the discretized weak form and highlight where these errors appear. Are the errors internally related?
- Derive the governing equations of solid mechanics (translational and rotational equilibrium) and the constitutive equations for linear elastic isotropic materials.
- Explain the difference between the concept of *p-refinement* and *h-refinement* in the Finite Element Method. What is the purpose of each? Illustrate on a freely chosen example what strategy you would use for mesh refinement in practice (sketching a drawing can be helpful).
- Restricting the validity of your developments to a linear structural behaviour, derive the expression of the stiffness matrix of a finite element starting from the discretized internal force vector: $\int_V [B]^T \{\sigma\} dV = \{f_{ext}\}$.
  Enumerate the simplifying assumptions required to result in a linear structural behaviour (small strains, small displacements, linear elastic material). Detail the steps to reach: $[K] = \int_V [B]^T [H] [B] dV$.
- Draw the main steps of the finite element code you implemented during the computer labworks in a flowchart (from the data reading to the result saving). Explain in a single sentence each step and highlight when the routine of a finite element is used in a given step.
- Derive the *Weak form* of the equilibrium equations starting from the *Strong form* ($b_i + \tau_{ij,j} = 0$). Explain the mathematical steps (use of virtual displacement, integration over volume, Gauss divergence theorem).
- Explain why the Finite Element solution generally converges towards the real solution from the "stiffer side" (i.e., the structure appears stiffer in FEA than in reality). How can a FE solution of a given problem be enhanced? Explain two approaches (h-refinement and p-refinement).

---

## Chapter 5 - Shape functions

### Multiple-choice questions
- The condition on a shape function for the correct representation of rigid body mode is: ($\sum_i N_i(\xi) = 1$ / $N_i(\xi_j) = \delta_{ij}$)
- The condition on shape functions for the local support condition is: ($\sum_i N_i(\xi) = 1$ / $N_i(\xi_j) = \delta_{ij}$)

### Open-ended questions
- Explain why the sum of the interpolation functions within the finite element must be equal to one.
- What are the requirements for the shape functions to ensure convergence of the Finite Element Method? (Discuss convergence criteria).
- What are the general properties of shape functions?
- Explain the continuity requirements in a finite element approach. Why do these requirements differ depending on the element type (e.g., bar vs beam)?
- Explain the difference between *Global approximation* and *Local approximation*. What is the effect of the order of the interpolation function on the goodness of the approximation?

---

## Chapter 6 - Isoparametric elements

### Multiple-choice questions
- A second-order isoparametric element (e.g. REM-8) can represent exactly a contour defined by straight lines. (Yes always / Yes but provided the mesh is fine enough / Yes but only in plane strain state / No)

### Open-ended questions
- Describe the isoparametric approach. Explain the concept of geometrical mapping to convert a reference element to the real element.
- In the context of isoparametric elements, demonstrate how to transform an integral defined in the physical $(x, y)$ axis into an integral over the natural coordinates ($-1, +1$).
- Explain the specific sources of errors related to shape functions. Discuss geometrical errors (mesh not coinciding with real geometry) and discretization errors (approximation of the function by a finite subset).

---

## Chapter 7 - Numerical integration

### Multiple-choice questions
- Hourglass nodes are non-zero energy modes for zero displacements. (True / False)
- Hourglass modes are zero energy modes for non-zero displacements. (True / False)

### Open-ended questions
- Explain, using equations, why the 4-noded tetrahedron and the 3-noded triangular finite elements do not require numerical integration to compute their stiffness matrix while the 8-noded quadrilateral finite element does.
- How do you perform numerical integration in FEM? Explain how to transform a continuous integral into a numerical summation, defining every term in the equation.
- Compare *Exact integration* and *Reduced integration*. What is the difference between the two methods? Why and when would you use one instead of the other?
- In the case of an isoparametric finite element, explain the main steps to compute the stiffness matrix when using numerical integration. Your explanation should include the equations used in the mathematical development. In which case does integration error appear?
- Explain the concept of *Spurious mechanism* (or Hourglass mode).
  a) How can reduced integration provoke this mechanism?
  b) What are the conditions to avoid singularity of the stiffness matrix?
  c) How does it propagate through the mesh?
- What is an integration error in the context of FEM?
- In Gauss-Legendre quadrature, how many integration points ($K$) are required to exactly integrate a polynomial of order $P$? (Explain the relation $P = 2K - 1$).
- Explain the concept of ill-conditioning or rank deficiency that can occur if the number of integration points is insufficient (e.g., when the number of unknowns is greater than the rank of the stiffness matrix).

---

## Chapter 8 - Bar and beam elements

### Multiple-choice questions
- In Timoshenko beam theory, the transverse sections which were initially planar: (remain planar after deformation / are not planar anymore after deformation)
- How many DOFs does a 2D Bernoulli Finite Element have? (6 / 12 / 8 / 9)
- In a 2D Bernoulli element, if $v$ refers to the transverse displacement, the curvature $1/R$ is given by: ($dv/dx$ / $1/(dv/dx)$ / $d^2v/dx^2$ / $1/(d^2v/dx^2)$)

### Open-ended questions
- Why does a secondary derivative ($d^2/dx^2$) come into play for beam elements formulation, and what does it imply for continuity requirements ($C^0$ vs $C^1$)?
- Describe the Finite Element formulation for beam elements.
- Perform the demonstration of the formulation for a 1D bar element.
- Describe the difference between *local* and *global* coordinates in the context of bar and truss elements.
- Write the matrix relationship between generalized strains ($\epsilon, 1/R$) and generalized stresses ($N, M$) for a Bernoulli beam element. Define the terms in the constitutive matrix $H$.

### Exercises
- A bar structural element with the following material and geometrical properties (Young's modulus $E$, density $\rho$, cross-sectional area $A$, length $L$) is subjected to body forces $f_x = \rho g$ and a concentrated load $F_x$ (applied at the end) with $g$ the gravitational acceleration along $x$.
  The analytical solution of the differential equation $\frac{d}{dx} \left( E \frac{du_x}{dx} \right) + f_x = 0$ in the volume of the bar $V$ is:
  $$ u_x = \frac{1}{AE} \left( -A\rho g \frac{x^2}{2} + (F_x + A\rho gL) x \right) $$
  Suppose that the structural bar is now discretized using a *single bar finite element*.
  a) Derive the finite element solution of the problem, including stiffness matrix, total externally applied forces and apply boundary conditions.
  b) Calculate the nodal displacement, the strain and stress field for finite element and compare them to the analytical solution. Visualize the results in graphs showing the displacement, the strain and the stress versus the element length, $x$.

---

## Chapter 9-10 - 2D and 3D finite elements

### Multiple-choice questions
- In the plane strain state, the axial strain $\epsilon_z$ (i.e. perpendicular to the thickness) is always equal to 0. (True / False)
- In the plane stress state, the axial strain $\epsilon_z$ (i.e. perpendicular to the thickness) is always equal to 0. (True / False)
- Which Finite Element is subject to shear locking in bending dominated problems? (4-noded quadrilateral / 6-noded triangle)
- How many degrees of freedom does a TRIM-6 element possess? (6 / 9 / 12 / 18)
- How many degrees of freedom does a TRIM-3 element possess? (3 / 6 / 12 / 24)
- What is the size of the shape function matrix $N(x,y)$ for a REM4 element? ($1 \times 4$ / $2 \times 4$ / $1 \times 8$ / $2 \times 8$)

### Open-ended questions
- How do we obtain the plane stress and plane strain matrices? Explain the difference.
- For a 3-noded triangular element (CST), derive the strain-displacement matrix ($B$) starting from the shape functions.
- Why are the strain and stresses constant in a Constant Strain Triangle (CST) element?
- Explain the phenomenon of *Shear Locking*.
  a) On which type of elements does it occur?
  b) Why does the Q4 (bilinear quadrilateral) element struggle with bending (pure bending example)?
  c) How does it affect convergence?
- When would you use 2D Finite Elements versus 3D Finite Elements? Justify your choice with a given example (e.g., a train wheel, a plate). Discuss the advantages (BCs, geometry) and disadvantages (computational effort) of 3D FE.
- Draw two sketches representative of plane stress and of plane strain conditions, respectively. State the mathematical conditions corresponding to plane strain and plane stress modelling assumptions.
  Derive from the 3D Hooke's elastic stiffness matrix ($H$) the 2D reduced material stiffness matrixes for both conditions.
- Which element converges faster: the Q4 (4-noded quadrilateral) or the T3 (3-noded triangle)? Discuss in terms of degrees of freedom (DOF) and accuracy.
- What are the differences between 3-noded (linear) and 4-noded (quadratic) finite elements (in 1D/2D context)? Do linear or quadratic elements converge faster and why?
- Explain which order of finite element (linear vs. quadratic) you would use in a 2D/3D structural simulation. Present the advantages and disadvantages of the available solutions (accuracy, computational cost, risk of locking).
- Explain how you would decide between *plane stress* or *plane strain* modelling assumptions for a given physical problem.

---

## Chapter 11 - Assembly and application of boundary conditions

### Open-ended questions
- Compare the *Direct method* and the *Penalty method* for applying boundary conditions.
  a) What is the difference between the two methods?
  b) When should you use one instead of the other?
  c) What are the advantages and disadvantages of each?
  d) Discuss the influence of the penalty parameter (too high vs too low).
- Regarding the *Penalty method*: How should the penalty parameter value, $Z$, be chosen? Give a recommendation for its value and justify your choice.
- Explain if different types of finite elements can be assembled in a single finite element problem (e.g., connecting a beam element to a solid element). Give examples to support your argumentation and discuss continuity requirements.

---

## Chapter 12 - System solving

### Multiple-choice questions
- The Cholesky decomposition takes advantages of which property of the stiffness matrix? (Large size / Symmetry)
- The computational effort for solving the discrete system of equilibrium equations is known a priori when using an iterative method. (True / False)
- Different direct solution methods of the discrete system of equilibrium equations yield the same result. (True / False)
- The computational effort/cost for solving the discrete system of equilibrium equations is known a priori when using direct methods. (True / False)
- Different iterative solution methods of the discrete system of equilibrium equations may give different results. (True / False)

### Open-ended questions
- Explain the concept of the frontal method for solving the discrete system of equilibrium equations. Illustrate using a figure how the approach works and list its advantages and disadvantages.
- Compare *Direct methods* and *Iterative methods* for solving the system of equations.
  a) Explain the concept of each (use a drawing if necessary).
  b) List the advantages and disadvantages of both approaches.

---

## Chapter 13 - Symmetry and modularity

### Multiple-choice questions
- A structural problem having a geometry with symmetry of revolution around an axis can always be modelled as axisymmetric. (True / False)
- A structural problem having a symmetry of the boundary conditions and loads can be modelled as symmetric. (True / False)
- Can the following structural problem be considered as axisymmetric? (Cooling tower in Drogenbos, loading = self-weight only, $z$ = vertical axis). (True / False)

### Open-ended questions
- What are the necessary conditions to have symmetry in a structural problem?
- If you have symmetry in your model, would you exploit it? Explain why or why not.
- In a symmetric problem, would you use a full model, a 1/2 model, or a smaller fraction? Explain your reasoning.
- Explain the concepts of **repeatability** and **modularity** in Finite Element Analysis. What is a **superelement**, how is it constructed, and what precautions should be taken regarding its boundaries?

---

## Chapter 14 - Verification and validation and error estimation

### Multiple-choice questions
- The concept of Validation corresponds to checking if the computational model is a reliable representation of: ("the physical model" / "the real life phenomena")

### Open-ended questions
- Elaborate on the difference between *Verification* and *Validation* in the context of finite element analysis.
- Why is the *Displacement Patch Test (DPT)* important in Finite Element Analysis?
- Explain how the *recovery-based error estimation* works. Specifically, explain the concept of approximating a smooth gradient (stress recovery).
- How does the error evolve when increasing the degrees of freedom (DOF) for linear versus quadratic elements? (Discuss convergence rates).
- Define *Convergence* in FEM. How do you perform a convergence study and what are the different steps involved?

---

## Chapter 15 - Plate theory

### Multiple-choice questions
- In Reissner-Mindlin plate theory, the sections which were initially plane: (remain plane after deformation / are not planar anymore after deformation)
- A membrane is a structure mostly subject to bending. (True / False)

### Open-ended questions
- Define the stress resultants (stress items) in plate theory.
- Establish the equilibrium equations of translation and rotation for a plate finite element.
- Why does the finite element formulation for Plate theory (specifically Kirchhoff plates) require $C^1$ continuity?
- Discuss the computational advantages of using Plate and Shell finite elements compared to full 3D solid elements. Under what geometric conditions are they applicable?
- Explain the advantages and limitations of using Beam and Plate finite elements compared to Solid finite elements. Explain how $C^1$ continuity intervenes in their construction and why it is required for these specific elements compared to standard solid elements.

---

## Chapter 16 - Structural dynamics

### Multiple-choice questions
- In structural dynamics, the units of the elements composing the damping matrix $C$ are: ($kg/s$ / $N/s$ / $N/m$ / $kg/m$)
- Viscous damping is proportional to: (the acceleration / the velocity)
- For safety reasons, bridge designers seek to achieve the lowest possible value of the first eigenfrequency. (True / False)

### Open-ended questions
- Define what an eigenmode and an eigenfrequency are for an undamped dynamic system both mathematically and as physical interpretation.
- What is a modal basis?
  a) What is the goal of using a modal basis?
  b) What is the starting system?
  c) Why do we perform structural dynamics analysis?
- Explain how to pass from the coupled system of equations of motion to uncoupled equations. (Hint: Modal superposition / Modal basis).

---

## Chapter 17 - Computational homogenization

### Multiple-choice questions
- In multiscale methodology, phenomenological constitutive laws are used in the: (Structural scale / Microstructural scale)

### Open-ended questions
- What is the purpose of *Computational Homogenization*? Give an example of an application where this method is used.
