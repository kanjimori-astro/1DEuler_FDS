# 1DEuler_FDS
The code solves the 1D Euler equations with the flux difference splitting method.

# Introduction
One-dimensional compressible flow is described by the Euler equation

$$\frac{\partial \vec{Q}}{\partial t}+\frac{\partial \vec{E}}{\partial t}=0,$$

where the $$\vec{Q}=(\rho, \rho u, e)^T$$ is the conserved variables and $$\vec{E}=(\rho u, p+\rho u^2, (e+p)u)^T$$ is the flux. The time evolution of the equation can be determined as

$$Q^{n+1}_ j=Q^n_j+\frac{\Delta t}{\Delta x}(\tilde{E}^n_{j+1/2}-\tilde{E}^n_{j-1/2}),$$

where $$\tilde{E}^n_{j+1/2}$$ is the numerical flux. In the flux difference splitting method, the numerical flux is given by

$$\tilde{E}^n_{j+1/2}=\frac{1}{2}(E_{j+1}+E_j)-|A|_ {j+1/2}(Q_{j+1}-Q_j)),$$

where $$|A|_ {j+1/2}=R_{j+1/2}|\Lambda|_ {j+1/2}R^{-1}_{j+1/2}$$, $$\Lambda$$ is the diagonalized flux Jacobian matrix, and $R$ is defined in Section 5.3 in Ref. [1]. The quantities at the cell edges are estimated with Roe's average.

# Test

I solve the Sod shock tube as a test problem. The spatial grid number is 500 and Courant number is 0.2.
![FDS](https://github.com/user-attachments/assets/06bd66f2-08bd-4ca6-bb17-f09cd1fee847)

# References
[1] 藤井孝藏 (1994), 『流体力学の数値計算法』, 東京大学出版会.
