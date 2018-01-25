"""
This is the main file for computing Fluid Structure Interaction problems from the master thesis
(https://www.duo.uio.no/handle/10852/57788).
FSI problems are solved in a monolithic environment, meaning all equations are solved at once.

The solver can easily be validated using the benchmark made by Stefan Turek and Jaroslav Hron, as computing domains in
the Mesh folder, inlet fluid velocities etc, are taken from that benchmark.
(https://www.researchgate.net/publication/226447172_Proposal_for_Numerical_Benchmarking_of_Fluid-Structure_Interaction_Between_an_Elastic_Object_and_Laminar_Incompressible_Flow)


Example:
    To run for instance the first FSI test fsi1 run the line:
    $ python fsi_monolithic_main.py -problem fsi1 -T 10 -dt 0.5

    This will run the fsi1 test with a simulation time of 10 seconds with a timestep of 0.5 seconds. The fsi1 problem
    can be found in the Problems folder, which is passed in to the general solver.

    One can also run the fsi2 test case with increased control over parameters:
    $ python fsi_monolithic_main.py -fluidvari thetaCN -solidvari thetaCN -extravari alfa -solver newtonsolver -problem fsi2
    -T 10 -dt 0.5 -theta 0.5

Its also always possible to get help by running python fsi_monolithic_main.py -h
Other FSI problems can be made by making a problem python file in the Problem folder,
and the subsequent mesh in the Mesh folder, passing in fluid and structure parameters, boundary conditions etc.

"""

from dolfin import *
from Utils.argpar import *

# Silence FEniCS output
set_log_active(False)

# Get user input
args = parse()

# Mesh refiner
exec("from Problems.%s import *" % args.problem)
if args.refiner is not None:
    for i in range(args.refiner):
        mesh = refine(mesh)

update_variables = {}

# Update argparser input
for key in args.__dict__:
    if args.__dict__[key] is not None:
        update_variables[key] = args.__dict__[key]

vars().update(update_variables)


# Import variationalform and solver
print(args.solver)
exec("from Fluidvariation.%s import *" % args.fluidvari)
exec("from Structurevariation.%s import *" % args.solidvari)
exec("from Extrapolation.%s import *" % args.extravari)
exec("from Newtonsolver.%s import *" % args.solver)


# Function Domains
D = VectorElement("CG", mesh_file.ufl_cell(), d_deg)
V = VectorElement("CG", mesh_file.ufl_cell(), v_deg)
P = FiniteElement("CG", mesh_file.ufl_cell(), p_deg)

# Define coefficients
k = Constant(dt)
n = FacetNormal(mesh_file)

# Sets up the functionspaces and functions in dictionaries
if args.extravari == "biharmonic":
    print("Biharmonic")
    DVP = FunctionSpace(mesh_file, MixedElement([D, V, P, D]))

    dvp_ = {}
    d_ = {}
    v_ = {}
    p_ = {}
    w_ = {}

    for time in ["n", "n-1", "n-2"]:
        dvp = Function(DVP)
        dvp_[time] = dvp
        d, v, p, w = split(dvp)

        d_[time] = d
        v_[time] = v
        p_[time] = p
        w_[time] = w

    phi, psi, gamma, beta = TestFunctions(DVP)

else:
    DVP = FunctionSpace(mesh_file, MixedElement([D, V, P]))

    # Create functions
    dvp_ = {}
    d_ = {}
    v_ = {}
    p_ = {}

    for time in ["n", "n-1", "n-2"]:
        dvp = Function(DVP)
        dvp_[time] = dvp
        d, v, p = split(dvp)

        d_[time] = d
        v_[time] = v
        p_[time] = p

    phi, psi, gamma = TestFunctions(DVP)

t = 0

# Check for available solvers:
for i in ["mumps", "superlu_dist", "default"]:
    if has_lu_solver_method(i):
        solver_method = i

up_sol = LUSolver(solver_method)
up_sol.parameters["same_nonzero_pattern"] = False
up_sol.parameters["reuse_factorization"] = False


vars().update(fluid_setup(**vars()))
vars().update(structure_setup(**vars()))
vars().update(extrapolate_setup(**vars()))
vars().update(solver_setup(**vars()))
vars().update(initiate(**vars()))
vars().update(create_bcs(**vars()))

atol = 1e-6
rtol = 1e-6
max_it = 100
lmbda = 1.0

dvp_res = Function(DVP)
chi = TrialFunction(DVP)

counter = 0

tic()
while t <= T + 1e-8:
    t += dt

    if MPI.rank(mpi_comm_world()) == 0:
        print("Solving for timestep %g" % t)

    pre_solve(**vars())
    vars().update(newtonsolver(**vars()))

    times = ["n-2", "n-1", "n"]
    for i, t_tmp in enumerate(times[:-1]):
        [t_tmp].vector().zero()
        dvp_[t_tmp].vector().axpy(1, dvp_[times[i+1]].vector())
    vars().update(after_solve(**vars()))
    counter += 1

simtime = toc()
print("Total Simulation time %g" % simtime)
post_process(**vars())
