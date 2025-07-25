# SmoothingMBA-QCQP-ℓ₁
This is a MATLAB package that implements the smoothing moving ball approximation (*s*MBA) method in [[1]] for solving the following quadratically constrained quadratic program with $\ell_1$ regularization (QCQP-ℓ₁):

$$\eqalign{
	\min\limits_{x\in ℝ^{n}} & x^\top Q_0 x + (q^{0})^{\top}x + \rho ‖x‖_1 \\
	\text{s.t.} & \frac{1}{2} x^\top Q_i x + (q^{i})^\top x - b_i \leq 0 \quad \forall\, i=1,\ldots,m,
}$$

where $Q_{i}\in 𝕊^{n}_{+}$, $q^{i}\in ℝ^{n}$, $i=0,1,\ldots,m$, $b\in ℝ^{m}$ and $\rho > 0$.
The performance of *s*MBA is compared with CVX (version 2.2) using the SDPT3 solver (version 4.0).
For further details, please refer to our paper in [[1]].

# Matlab source codes
- **demo_QCQP.m**\
A demo of the numerical experiments in [[1]].

- **QCQP.m**\
A function that generates the data of a convex QCQP.

- **sMBA.m**\
The implementation of sMBA for solving QCQP-ℓ₁.

- **CaseSg.m** and **SubP_alpha.m**\
The implementations of the root-finding scheme described in [[2], Appendix A], which we use to solve the subproblem of *s*MBA (i.e., [[1], (3.1)]). The codes are hosted by Pong T. K. and are available at [https://www.polyu.edu.hk/ama/profile/pong/MBA_l1vl2/.](https://www.polyu.edu.hk/ama/profile/pong/MBA_l1vl2/)

# References
[1]: https://arxiv.org/pdf/2505.12314 "J. Xu, T. K. Pong and N. S. Sze. A smoothing moving balls approximation method for a class of conic-constrained difference-of-convex optimization problems. Preprint (2025)."
\[1\] [J. Xu, T. K. Pong and N. S. Sze. A smoothing moving balls approximation method for a class of conic-constrained difference-of-convex optimization problems. Preprint (2025).](https://arxiv.org/pdf/2505.12314)

[2]: https://epubs.siam.org/doi/abs/10.1137/20M1314057 "P. Yu, T. K. Pong and Z. Lu. Convergence rate analysis of a sequential convex programming 
method with line search for a class of constrained difference-of-convex optimization problems. *SIAM J. Optim.* 31, pp. 2024--2054 (2021)."
\[2\] [P. Yu, T. K. Pong and Z. Lu. Convergence rate analysis of a sequential convex programming method with line search for a class of constrained difference-of-convex optimization problems. *SIAM J. Optim.* 31, pp. 2024--2054 (2021).](https://epubs.siam.org/doi/abs/10.1137/20M1314057)



