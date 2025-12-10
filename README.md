# Introduction
This repository is dedicated to implementing two QRD-based versions of the EX-RLS algorithm: one using Givens rotations and the other using Householder reflectors.
For evaluating the algorithms, a system identification setup was adopted, where the target system to be identified is a Rayleigh fading channel, as illustrated in the figure below.

<img width="470" height="301" alt="image" src="https://github.com/user-attachments/assets/1dcc1b16-8e05-4ec5-9d50-a9137617d8a8" />

# QR-EX-RLS

The two implemented algorithms are presented below.

Givens-based EX-RLS:

<img width="813" height="761" alt="image" src="https://github.com/user-attachments/assets/45330f9c-3ac7-405f-8273-7228888c7094" />

Householder-based EX-RLS:

<img width="923" height="588" alt="image" src="https://github.com/user-attachments/assets/2f64d392-c612-4332-b963-21997ed05b92" />

For the Cholesky decomposition stage, one may use MATLAB’s built-in chol function or the custom implementation shown below:

<img width="972" height="667" alt="image" src="https://github.com/user-attachments/assets/10584f12-8c57-4366-a12c-07c90ce7cf61" />


# How to use

Download this repository and leave all the codes in the same folder. 
For DPD applications use isaacmacario2/DPD repository (https://github.com/isaacmacario2/DPD).

# Main references

J. A. Apolinário Jr., "QRD-RLS Adaptive Filtering". New York: Springer, 2009.

W.Liu, J. C. Príncipe and S. Haykin, "Kernel Adaptive Filtering: A comprehensive introduction". 1. ed. New Jersey: Wiley, 2010.

A. Sayed, \textit{Adaptive Filters}. 1. ed. New Jersey: Wiley, 2008.

S. Haykin, “Adaptive Filter Theory”, 1st ed., New Jersey: Prentice Hall, 1986.

P. S. R. Diniz, "Adaptive Filtering: Algorithms and practical implementation". 4. ed. New York: Springer, 2013.

P. S. R. Diniz, M. L. R. Campos,  W. A. Martins, M. V. S. Lima and J. A. Apolinário Jr., "Online Learning and Adaptive Filters". New York: Wiley, 2023.

William Ford, “Numerical Linear Algebra with Applications,” Elsevier, London, 2015.

Gene H. Golub and Charles F. Van Loan, “Matrix Computations,” The Johns Hopkins University Press, Baltimore, Maryland, 4th edition, 2013.

# How to Cite
@article{QREXRLS,
  author       = {Gouveia, I. M. S. and Apolin\'{A}rio Jr., J. A. and Saunders Filho, C. A. B. and Ramos, A. L. L.},
  title        = {Stable EX-RLS Algorithms: a QR Decomposition Approach},
  journal      =, 
  volume       =,
  number       =,
  pages        =, 
  doi          = 
}
