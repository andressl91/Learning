from dolfin import *

def solver_setup(F_fluid_linear, F_fluid_nonlinear, \
                 F_solid_linear, F_solid_nonlinear, DVP, dvp_, **monolithic):

    F_lin = F_fluid_linear + F_solid_linear
    F_nonlin = F_fluid_nonlinear + F_solid_nonlinear
    F = F_lin + F_nonlin

    chi = TrialFunction(DVP)
    Jac = derivative(F, dvp_["n"], chi)

    A = None
    b = None

    return dict(F=F, Jac=Jac, A=A, b=b)


def newtonsolver(F, J_nonlinear, A_pre, A, b, bcs, \
                dvp_, up_sol, dvp_res, rtol, atol, max_it, T, t, **monolithic):
    Iter      = 0
    residual   = 1
    rel_res    = residual
    lmbda = 1

    while rel_res > rtol and residual > atol and Iter < max_it:
        if Iter % 10 == 0:
            A = assemble(Jac, tensor=A)
            A.ident_zeros()

        b = assemble(-F, tensor=b)

        [bc.apply(A, b, dvp_["n"].vector()) for bc in bcs]
        up_sol.solve(A, dvp_res.vector(), b)
        dvp_["n"].vector().axpy(lmbda, dvp_res.vector())
        [bc.apply(dvp_["n"].vector()) for bc in bcs]
        rel_res = norm(dvp_res, 'l2')
        residual = b.norm('l2')

        if MPI.rank(mpi_comm_world()) == 0:
            print "Newton iteration %d: r (atol) = %.3e (tol = %.3e), r (rel) = %.3e (tol = %.3e) " \
        % (Iter, residual, atol, rel_res, rtol)
        Iter += 1

    return dict(t=t)
