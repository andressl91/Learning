from dolfin import *

def solver_setup(F_fluid_linear, F_fluid_nonlinear, \
                 F_solid_linear, F_solid_nonlinear, DVP, dvp_, **monolithic):

    F_lin = F_fluid_linear + F_solid_linear
    F_nonlin = F_fluid_nonlinear + F_solid_nonlinear
    F = F_lin + F_nonlin

    chi = TrialFunction(DVP)
    Jac = derivative(F, dvp_["n"], chi)

    return dict(F=F, Jac=Jac)


def newtonsolver(F, Jac, bcs,\
                dvp_, up_sol, dvp_res, rtol, atol, max_it, T, t, **monolithic):
    Iter      = 0
    residual   = 1
    rel_res    = residual
    lmbda = 1
    timeAssemble = []
    fulltime = []
    timeAssemble2 = []
    solvetime = []


    while rel_res > rtol and residual > atol and Iter < max_it:
        #if Iter % 5 == 0:
        timefull0 = timeme.time()
        time0 = timeme.time()
        A = assemble(Jac, form_compiler_parameters = {"quadrature_degree": 4})
        timeAssemble.append(timeme.time() - time0)
        #A = assemble(J_nonlinear, tensor=A) #keep_diagonal = True
        A.axpy(1.0, A_pre, True)
        A.ident_zeros()

        time1 = timeme.time()
        b = assemble(-F, tensor=b)
        timeAssemble2.append(timeme.time() - time1)

        time2 = timeme.time()
        [bc.apply(A, b, dvp_["n"].vector()) for bc in bcs]
        up_sol.solve(A, dvp_res.vector(), b)
        solvetime.append(timeme.time() - time2)

        dvp_["n"].vector().axpy(lmbda, dvp_res.vector())
        [bc.apply(dvp_["n"].vector()) for bc in bcs]
        rel_res = norm(dvp_res, 'l2')
        residual = b.norm('l2')

        fulltime.append(timeme.time() - timefull0)


        if isnan(rel_res) or isnan(residual):
            print "type rel_res: ",type(rel_res)
            t = T*T

        if MPI.rank(mpi_comm_world()) == 0:
            print "Newton iteration %d: r (atol) = %.3e (tol = %.3e), r (rel) = %.3e (tol = %.3e) " \
        % (Iter, residual, atol, rel_res, rtol)
        Iter += 1
    

    return dict(t=t)
