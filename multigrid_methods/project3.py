""" 
Creator: Letizia D'Achille 

Functions:

	full_weighting_1d(vh, v2h)
	full_weighting_2d(vh, v2h)
		Implementation of the full weighting restriction operator in 1d/2d
	prolongation_1d(vh, v2h)
	prolongation_2d(vh, v2h)
		Implementation of the prolongation operator in 1d/2d
	jacobi_step_1d(uh, fh, omega)
	jacobi_step_2d(uh, fh, omega)
		Performs one step of the weighted Jacobi method in 1d/2d
	two_grid_correction_step_1d(uh, fh, omega)
	two_grid_correction_step_2d(uh, fh, omega)
		Performs one step of the two-grid correction scheme in 1d/2d
	w_cycle_step_1d(uh, fh, omega, alpha1, alpha2)
	w_cycle_step_2d(uh, fh, omega, alpha1, alpha2)
		Performs one step of the W-cycle method in 1d/2d
	full_mg_1d(uh, fh, omega, alpha1, alpha2, nu)
	full_mg_2d(uh, fh, omega, alpha1, alpha2, nu)
		Performs one step of the full multigrid method in 1d/2d
"""

import numpy as np
import numpy.linalg as la

def full_weighting_1d(vh, v2h):
	"""
	Implementation of the full weighting restriction operator in 1d

		Parameters:
			vh (1d array): An array defined over the fine grid
			v2h (1d array): An array defined over the coarse grid
	"""
	n2 = len(v2h) - 1
	v2h[1:n2] = (vh[1:2*n2-2:2] + 2*vh[2:2*n2-1:2] + vh[3:2*n2:2]) / 4

def full_weighting_2d(vh, v2h):
	"""
	Implementation of the full weighting restriction operator in 2d

		Parameters:
			vh (2d array): An array defined over the fine grid
			v2h (2d array): An array defined over the coarse grid
	"""
	n2 = len(v2h) - 1
	v2h[1:n2,1:n2] = (vh[1:2*n2-2:2, 1:2*n2-2:2] + vh[1:2*n2-2:2, 3:2*n2:2] + vh[3:2*n2:2, 1:2*n2-2:2] + vh[3:2*n2:2, 3:2*n2:2] + 2*(vh[2:2*n2-1:2, 1:2*n2-2:2] + vh[2:2*n2-1:2, 3:2*n2:2] + vh[1:2*n2-2:2, 2:2*n2-1:2] + vh[3:2*n2:2, 2:2*n2-1:2]) + 4*vh[2:2*n2-1:2, 2:2*n2-1:2]) / 16

def prolongation_1d(vh, v2h):
	"""
	Implementation of the full prolongation operator in 1d

		Parameters:
			vh (1d array): An array defined over the fine grid
			v2h (1d array): An array defined over the coarse grid
	"""
	n2 = len(v2h) - 1
	vh[0:2*n2-1:2] = v2h[0:n2]
	vh[1:2*n2:2] = (v2h[0:n2] + v2h[1:n2+1]) / 2

def prolongation_2d(vh, v2h):
	"""
	Implementation of the full prolongation operator in 2d

		Parameters:
			vh (2d array): An array defined over the fine grid
			v2h (2d array): An array defined over the coarse grid
	"""
	n2 = len(v2h) - 1
	vh[0:2*n2-1:2, 0:2*n2-1:2] = v2h[0:n2,0:n2]
	vh[1:2*n2:2, 0:2*n2-1:2] = (v2h[0:n2,0:n2] + v2h[1:n2+1,0:n2]) / 2
	vh[0:2*n2-1:2, 1:2*n2:2] = (v2h[0:n2,0:n2] + v2h[0:n2,1:n2+1]) / 2
	vh[1:2*n2:2, 1:2*n2:2] = (v2h[0:n2,0:n2] + v2h[1:n2+1,0:n2] + v2h[0:n2,1:n2+1] + v2h[1:n2+1,1:n2+1]) / 4

def jacobi_step_1d(uh, fh, omega):
	"""
	Performs one step of the weighted Jacobi method in 1d

		Parameters:
			uh (1d array): Initial guess of the solution
			fh (1d array): Right hand side of the linear system
			omega (float): Weight applied to Jacobi method
		
		Returns:
			res (float): Pseudo-residual of one step of the method
	"""
	n = len(fh) - 1
	h = 1/n
	x = np.linspace(0,1,n+1)
	c = 2 - abs(x)

	uh1 = uh.copy()
	uh[1:n] = (1 - omega) * uh1[1:n] + omega * (h**2 * fh[1:n] + uh1[0:n-1] + uh1[2:n+1]) / (2 + c[1:n] * h**2)
	res = la.norm(uh - uh1, np.inf)
	return res

def two_grid_correction_step_1d(uh, fh, omega):
	"""
	Performs one step of the two-grid correction scheme in 1d
	Relaxing is implemented using Jacobi method

		Parameters:
			uh (1d array): Initial guess of the solution
			fh (1d array): Right hand side of the linear system
			omega (float): Weight applied to Jacobi method
		
		Returns:
			res (float): Pseudo-residual of one step of the method
	"""
	n = len(fh) - 1
	h = 1/n
	x = np.linspace(0,1,n+1)
	c = 2 - abs(x)
	
	uh1 = uh.copy()
	jacobi_step_1d(uh1, fh, omega)
	rh = uh1.copy()
	rh[1:n] = fh[1:n] + (uh1[0:n-1] + uh1[2:n+1] - (2 + c[1:n] * h**2) * uh1[1:n]) / h**2
	n2 = n//2
	r2h = np.zeros(n2 + 1)
	full_weighting_1d(rh, r2h)
	e2h = np.zeros_like(r2h)
	for _ in range(10):
		jacobi_step_1d(e2h, r2h, omega)
	eh = np.zeros_like(x)
	prolongation_1d(eh, e2h)
	uh1 = uh1 + eh
	jacobi_step_1d(uh1, fh, omega)
	res = la.norm(uh - uh1, np.inf)
	uh[:] = uh1[:]
	return res

def w_cycle_step_1d(uh, fh, omega, alpha1, alpha2):
	"""
	Performs one step of the W-cycle method in 1d
	Relaxing is implemented using Jacobi method

		Parameters:
			uh (1d array): Initial guess of the solution
			fh (1d array): Right hand side of the linear system
			omega (float): Weight applied to Jacobi method
			alpha1 (int): Number of pre-smoothing steps
			alpha2 (int): Number of post-smoothing steps
		
		Returns:
			res (float): Pseudo-residual of the final post-smoothing step
	"""
	n = len(fh) - 1
	h = 1/n
	x = np.linspace(0,1,n+1)
	c = 2 - abs(x)
	
	if (n == 2):
		uh[1] = fh[1] * h**2 / (2 + c[1] * h**2)
		return 0
	
	uh1 = uh.copy()
	for _ in range(alpha1):
		jacobi_step_1d(uh1, fh, omega)
	rh = uh1.copy()
	rh[1:n] = fh[1:n] + (uh1[0:n-1] + uh1[2:n+1] - (2 + c[1:n] * h**2) * uh1[1:n]) / h**2
	n2 = n//2
	r2h = np.zeros(n2 + 1)
	full_weighting_1d(rh, r2h)
	e2h = np.zeros_like(r2h)
	for _ in range(2):
		w_cycle_step_1d(e2h, r2h, omega, alpha1, alpha2)
	eh = np.zeros_like(x)
	prolongation_1d(eh, e2h)
	uh1 = uh1 + eh
	for _ in range(alpha2):
		res = jacobi_step_1d(uh1, fh, omega)
	uh[:] = uh1[:]
	return res

def full_mg_1d(uh, fh, omega, alpha1, alpha2, nu):
	"""
	Performs one step of the full multigrid method in 1d
	Relaxing is implemented using Jacobi method

		Parameters:
			uh (1d array): Initial guess of the solution
			fh (1d array): Right hand side of the linear system
			omega (float): Weight applied to Jacobi method
			alpha1 (int): Number of pre-smoothing steps
			alpha2 (int): Number of post-smoothing steps
			nu (int): Number of W-cycle steps
		
		Returns:
			res (float): Result of the final step of W-cycle method
	"""
	n = len(fh) - 1
	h = 1/n
	x = np.linspace(0,1,n+1)
	c = 2 - abs(x)

	if (n == 2):
		uh[1] = fh[1] * h**2 / (2 + c[1] * h**2)
		return 0
	
	n2 = n//2
	f2h = np.zeros(n2 + 1)
	full_weighting_1d(fh, f2h)
	u2h = np.zeros_like(f2h)
	full_mg_1d(u2h, f2h, omega, alpha1, alpha2, nu)
	uh1 = np.zeros_like(x)
	prolongation_1d(uh1, u2h)
	for _ in range(nu):
		res = w_cycle_step_1d(uh1, fh, omega, alpha1, alpha2)
	uh[:] = uh1[:]
	return res

def jacobi_step_2d(uh, fh, omega):
	"""
	Performs one step of the weighted Jacobi method in 2d

		Parameters:
			uh (2d array): Initial guess of the solution
			fh (2d array): Right hand side of the linear system
			omega (float): Weight applied to Jacobi method
		
		Returns:
			res (float): Pseudo-residual of one step of the method
	"""
	n = len(fh) - 1
	h = 1/n
	x = np.linspace(0, 1, n+1)
	y = np.linspace(0, 1, n+1)
	x, y = np.meshgrid(x, y)
	c = 2 - np.sqrt(x**2 + y**2)
	
	uh1 = uh.copy()
	uh[1:n,1:n] = (1 - omega) * uh1[1:n,1:n] + omega * (h**2 * fh[1:n,1:n] + uh1[0:n-1,1:n] + uh1[2:n+1,1:n] + uh1[1:n,0:n-1] + uh1[1:n,2:n+1]) / (4 + c[1:n,1:n] * h**2)
	res = la.norm((uh - uh1).flat, np.inf)
	return res

def two_grid_correction_step_2d(uh, fh, omega):
	"""
	Performs one step of the two-grid correction scheme in 2d
	Relaxing is implemented using Jacobi method

		Parameters:
			uh (2d array): Initial guess of the solution
			fh (2d array): Right hand side of the linear system
			omega (float): Weight applied to Jacobi method
		
		Returns:
			res (float): Pseudo-residual of one step of the method
	"""
	n = len(fh) - 1
	h = 1/n
	x = np.linspace(0, 1, n+1)
	y = np.linspace(0, 1, n+1)
	x, y = np.meshgrid(x, y)
	c = 2 - np.sqrt(x**2 + y**2)

	uh1 = uh.copy()
	jacobi_step_2d(uh1, fh, omega)
	rh = uh1.copy()
	rh[1:n,1:n] = fh[1:n,1:n] + (uh1[0:n-1,1:n] + uh1[2:n+1,1:n] + uh1[1:n,0:n-1] + uh1[1:n,2:n+1] - (4 + c[1:n,1:n] * h**2) * uh1[1:n,1:n]) / h**2
	n2 = n//2
	r2h = np.zeros((n2 + 1,n2 + 1))
	full_weighting_2d(rh, r2h)
	e2h = np.zeros_like(r2h)
	for _ in range(10):
		jacobi_step_2d(e2h, r2h, omega)
	eh = np.zeros_like(x)
	prolongation_2d(eh, e2h)
	uh1 = uh1 + eh
	jacobi_step_2d(uh1, fh, omega)
	res = la.norm((uh - uh1).flat, np.inf)
	uh[:] = uh1[:]
	return res

def w_cycle_step_2d(uh, fh, omega, alpha1, alpha2):
	"""
	Performs one step of the W-cycle method in 2d
	Relaxing is implemented using Jacobi method

		Parameters:
			uh (2d array): Initial guess of the solution
			fh (2d array): Right hand side of the linear system
			omega (float): Weight applied to Jacobi method
			alpha1 (int): Number of pre-smoothing steps
			alpha2 (int): Number of post-smoothing steps
		
		Returns:
			res (float): Pseudo-residual of the final post-smoothing step
	"""
	n = len(fh) - 1
	h = 1/n
	x = np.linspace(0, 1, n+1)
	y = np.linspace(0, 1, n+1)
	x, y = np.meshgrid(x, y)
	c = 2 - np.sqrt(x**2 + y**2)

	if (n == 2):
		uh[1,1] = fh[1,1] * h**2 / (4 + c[1,1] * h**2)
		return 0
	
	uh1 = uh.copy()
	for _ in range(alpha1):
		jacobi_step_2d(uh1, fh, omega)
	rh = uh1.copy()
	rh[1:n,1:n] = fh[1:n,1:n] + (uh1[0:n-1,1:n] + uh1[2:n+1,1:n] + uh1[1:n,0:n-1] + uh1[1:n,2:n+1] - (4 + c[1:n,1:n] * h**2) * uh1[1:n,1:n]) / h**2
	n2 = n//2
	r2h = np.zeros((n2 + 1,n2 + 1))
	full_weighting_2d(rh, r2h)
	e2h = np.zeros_like(r2h)
	for _ in range(2):
		w_cycle_step_2d(e2h, r2h, omega, alpha1, alpha2)
	eh = np.zeros_like(x)
	prolongation_2d(eh, e2h)
	uh1 = uh1 + eh
	for _ in range(alpha2):
		res = jacobi_step_2d(uh1, fh, omega)
	uh[:] = uh1[:]
	return res

def full_mg_2d(uh, fh, omega, alpha1, alpha2, nu):
	"""
	Performs one step of the full multigrid method in 2d
	Relaxing is implemented using Jacobi method

		Parameters:
			uh (2d array): Initial guess of the solution
			fh (2d array): Right hand side of the linear system
			omega (float): Weight applied to Jacobi method
			alpha1 (int): Number of pre-smoothing steps
			alpha2 (int): Number of post-smoothing steps
			nu (int): Number of W-cycle steps
		
		Returns:
			res (float): Result of the final step of W-cycle method
	"""
	n = len(fh) - 1
	h = 1/n
	x = np.linspace(0, 1, n+1)
	y = np.linspace(0, 1, n+1)
	x, y = np.meshgrid(x, y)
	c = 2 - np.sqrt(x**2 + y**2)

	if (n == 2):
		uh[1,1] = fh[1,1] * h**2 / (4 + c[1,1] * h**2)
		return 0
	
	n2 = n//2
	f2h = np.zeros((n2 + 1,n2 + 1))
	full_weighting_2d(fh, f2h)
	u2h = np.zeros_like(f2h)
	full_mg_2d(u2h, f2h, omega, alpha1, alpha2, nu)
	uh1 = np.zeros_like(x)
	prolongation_2d(uh1, u2h)
	for _ in range(nu):
		res = w_cycle_step_2d(uh1, fh, omega, alpha1, alpha2)
	uh[:] = uh1[:]
	return res