import sys
from deflation import ForwardProblem
from dolfin import *
from math import log

set_log_level(ERROR)

"""10.1002/cpa.20005"""

args = [sys.argv[0]] + """
                       --petsc.snes_max_it 200
                       --petsc.snes_monitor
                       --petsc.snes_converged_reason
                       --petsc.snes_stol 0.0
                       --petsc.snes_rtol 0.0
                       --petsc.snes_type newtonls
                       --petsc.snes_linesearch_type basic

                       --petsc.ksp_type preonly
                       --petsc.ksp_monitor_cancel
                       --petsc.ksp_rtol 1.0e-12
                       --petsc.ksp_atol 1.0e-12

                       --petsc.inner_pc_type lu
                       --petsc.inner_pc_gamg_type agg
                       --petsc.inner_pc_gamg_verbose 10
                       --petsc.inner_pc_gamg_coarse_eq_limit 2000
                       --petsc.inner_pc_gamg_agg_nsmooths 4
                       --petsc.inner_pc_gamg_threshold 0.04
                       --petsc.inner_pc_gamg_sym_graph true
                       --petsc.inner_mg_coarse_pc_type redundant
                       --petsc.inner_mg_coarse_sub_pc_type lu
                       --petsc.inner_mg_coarse_sub_pc_factor_shift_type NONZERO
                       --petsc.inner_mg_coarse_sub_pc_factor_shift_amount 1.0e-12
                       --petsc.inner_mg_levels_pc_type sor
                       """.split()
parameters.parse(argv=args)

mesh = UnitSquareMesh(100, 100, "crossed")
size = MPI.size(mpi_comm_world())
if size > 1:
  nrefine = log(size, 4)
  if int(nrefine) != nrefine:
    print "Need to have processors a power of 4, as each refinement multiplies work by 4"
    assert False

  for i in range(int(nrefine) + 3):
    mesh = refine(mesh, redistribute=False)

V = FunctionSpace(mesh, "CG", 1)
Vdim = V.dim()
if MPI.rank(mpi_comm_world()) == 0:
    print "Degrees of freedom: ", Vdim
    print "Degrees of freedom per core: ", Vdim/float(size)

delta = Constant(0.04)

u = Function(V)

v = TestFunction(V)
F = delta * inner(grad(v), grad(u))*dx + 1.0/delta * inner(v, u**3 - u)*dx

bcs = [DirichletBC(V, +1.0, "x[0] == 0.0 || x[0] == 1"),
       DirichletBC(V, -1.0, "x[1] == 0.0 || x[1] == 1")]

def empty_vector(model):
  # should be able to do return Vector(model.size()) but it doesn't work in parallel
  # dolfin's Vector API is terrible.
  b = Vector(model)
  b.zero()
  return b

class DeflationOperator(object):
    def __init__(self, y, power, shift):
        self.power = power
        self.shift = shift
        self.y = y
        self.roots = []

    def deflate(self, root):
        self.roots.append(root.copy(deepcopy=True))

    def normsq(self, y, root):
        return inner(y - root, y - root)*dx

    def evaluate(self):
        m = 1.0
        for normsq in [assemble(self.normsq(self.y, root)) for root in self.roots]:
            factor = normsq**(-self.power/2.0) + self.shift
            m *= factor

        return m

    def derivative(self):
        if len(self.roots) == 0:
            deta = empty_vector(self.y.vector())
            deta.zero()
            return deta

        p = self.power
        factors  = []
        dfactors = []
        dnormsqs = []
        normsqs  = []

        for root in self.roots:
            form = self.normsq(self.y, root)
            normsqs.append(assemble(form))
            dnormsqs.append(assemble(derivative(form, self.y)))

        for normsq in normsqs:
            factor = normsq**(-p/2.0) + self.shift
            dfactor = (-p/2.0) * normsq**((-p/2.0) - 1.0)

            factors.append(factor)
            dfactors.append(dfactor)

        eta = product(factors)

        deta = empty_vector(self.y.vector())
        deta.zero()

        for (solution, factor, dfactor, dnormsq) in zip(self.roots, factors, dfactors, dnormsqs):
            deta.axpy((eta/factor)*dfactor, dnormsq)

        return deta

def newton(F, u, bcs, atol=1.0e-8, deflation=None):

    [bc.apply(u.vector()) for bc in bcs]

    def norm(F, u, bcs): # defines convergence in nonlinear iteration
        b = assemble(F)
        [bc.apply(b, u.vector()) for bc in bcs]
        return b.norm("l2")

    def printnorm(i, n):
        print "%3d SNES Function norm %1.15e" % (i, n)

    du = Function(u.function_space())
    J = derivative(F, u, du)
    hbcs = [DirichletBC(bc) for bc in bcs]; [hbc.homogenize() for hbc in hbcs]
    i = 0

    n = norm(F, u, bcs)
    printnorm(i, n)

    while n > atol:
        du.vector().zero()

        solve(J + F == 0, du, hbcs)

        if deflation is not None:
            Eu = deflation.derivative().inner(du.vector())
            minv = 1.0 / deflation.evaluate()
            tau = (1 + minv*Eu/(1 - minv*Eu))
            du.assign(tau * du)

        u.assign(u + du)

        i = i + 1
        n = norm(F, u, bcs)
        printnorm(i, n)
       

roots = []
pvd = File("output/%s/roots.pvd" % float(delta))
natol = 1.0e-8
power = 2
shift = 1
d = DeflationOperator(u, power, shift)

for i in range(3):
  if i == 0:
      newton(F, u, bcs, natol)
  else:
      newton(F, u, bcs, natol, d)
  roots.append(u.copy(deepcopy=True))
  d.deflate(u)

  distances = []
  for otherroot in roots[:-1]:
    dist = assemble(inner(u - otherroot, u - otherroot)*dx)
    distances.append(dist)
  if len(distances) > 0:
    print "Distances to other roots: %s" % distances

  u.assign(Constant(0))
  pvd << roots[-1]
