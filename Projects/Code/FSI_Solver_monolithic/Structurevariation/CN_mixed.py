from dolfin import Constant, inner, inv, dot, grad, det, Identity,\
solve, lhs, rhs, assemble, DirichletBC, div, sym, tr, norm, \
MPI, mpi_comm_world
#from semi_implicit import *

def F_(U):
	return Identity(len(U)) + grad(U)

def J_(U):
	return det(F_(U))

def E(U):
	return 0.5*(F_(U).T*F_(U) - Identity(len(U)))

def S(U,lamda_s,mu_s):
    I = Identity(len(U))
    return 2*mu_s*E(U) + lamda_s*tr(E(U))*I

def Piola1(U,lamda_s,mu_s):
	return F_(U)*S(U,lamda_s,mu_s)


def structure_setup(d_, v_, p_, phi, psi, gamma, dS, mu_f, n,\
            dx_s, dx_f, mu_s, rho_s, lamda_s, k, mesh_file, **semimp_namespace):
	delta = 1E10
	alfa = 0.01
	h_ = mesh_file.hmin()
	F_solid_linear = rho_s/k*inner(v_["n"] - v_["n-1"], psi)*dx_s
	# Stress tensor
	#F_solid += 0.5*inner(Piola1(d_["n"], lamda_s, mu_s), grad(phi))*dx_s
	#F_solid += 0.5*inner(Piola1(d_["n-1"], lamda_s, mu_s), grad(phi))*dx_s
	F_solid_nonlinear = inner(Piola1(0.5*(d_["n"] + d_["n-1"]), lamda_s, mu_s), grad(psi))*dx_s

	#Deformation relation to velocity
	F_solid_linear += delta*((1./k)*inner(d_["n"] - d_["n-1"],phi)*dx_s - inner(0.5*(v_["n"] + v_["n-1"]), phi)*dx_s)
	#F_w = delta*((1.0/k)*inner(d-d0,psi)*dx_s - inner(0.5*(u+u0),psi)*dx_s)

	return dict(F_solid_linear = F_solid_linear, F_solid_nonlinear = F_solid_nonlinear)
