import sys; sys.path.append('..')
from deflation import ForwardProblem
from petscsnessolver import PetscSnesSolver

from dolfin import *

import gc
gc.disable()
import resource

import glob
from utils import prevent_MPI_output

prevent_MPI_output()

parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"
parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"
parameters["allow_extrapolation"] = True

args = [sys.argv[0]] +    """
                       --petsc.snes_max_it 100
                       --petsc.snes_type newtonls
                       --petsc.snes_linesearch_type l2
                       --petsc.snes_linesearch_max_it 1
                       --petsc.snes_linesearch_maxstep 1.0
                       --petsc.snes_stol 0.0
                       --petsc.snes_atol 1.0e-9
                       --petsc.snes_rtol 0.0
                       --petsc.snes_monitor
                       --petsc.snes_converged_reason
                       --petsc.snes_linesearch_monitor

                       --petsc.ksp_monitor_cancel
                       --petsc.ksp_monitor_true_residualx
                       --petsc.ksp_type preonly
                       --petsc.ksp_atol 1.0e-12
                       --petsc.ksp_rtol 1.0e-10

                       --petsc.inner_pc_type lu
                        """.split()
parameters.parse(args)

omega  = Constant(0.2)
mu     = Constant(0.18)
mu_end = 0.401

def forward(Z):
  phi = Function(Z)
  (r, c) = split(phi)
  (r_test, c_test) = split(TestFunction(Z))

  mesh = Z.mesh()
  x = SpatialCoordinate(mesh)

  mag = r**2 + c**2

  F = (
        mu*inner(r, r_test)*dx
      - 0.5*inner(grad(r), grad(r_test))*dx
      - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner(r, r_test)*dx
      - mag*inner(r, r_test)*dx

      + mu*inner(c, c_test)*dx
      - 0.5*inner(grad(c), grad(c_test))*dx
      - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner(c, c_test)*dx
      - mag*inner(c, c_test)*dx
      )

  bcs = [DirichletBC(Z, (0.0, 0.0), "on_boundary")]

  return (F, phi, bcs)

mesh = RectangleMesh(Point(-10, -10), Point(10, 10), 25, 25, "crossed")
V = FunctionSpace(mesh, "CG", 1)
Z = VectorFunctionSpace(mesh, "CG", 1, dim=2)
(F, phi, bcs) = forward(Z)

power = 2.0
shift = 1.0
solver = PetscSnesSolver()

phiex = Function(Z, name="Solution")
zero  = Function(Z, name="Zero")
groupphi = Function(Z, name="Negation")
groupsub = Function(V, name="Argh")

roots = {}

def guesses(prev_roots):
    for (i, root) in enumerate(prev_roots):
        root.label = "prev-root-%d" % i

    out = prev_roots
    if len(out) == 0:
        out = [interpolate(Constant((1, 1)), Z)]
        out[0].label = "initial-guess-1"
    return out

checkpoint = False
if not checkpoint:
  prev_mu_val = "init"
  roots[prev_mu_val] = []
  start = 0.4
  step = 0.001
else:
  prev_mu_val = 1.2400
  step = 0.001
  roots[prev_mu_val] = []
  for i in range(len(glob.glob("output-%2.4f/root-*.xml.gz" % prev_mu_val))):
      roots[prev_mu_val].append(Function(Z, "output-%2.4f/root-%d.xml.gz" % (prev_mu_val, i)))

  start = prev_mu_val + step

def nroots(mu):
    if mu <= 0.6:
        return 6
    return float("inf")

class H1ForwardProblemMagnitude(ForwardProblem):
    def norm(self, u, v):
        (ur, uc) = split(u)
        (vr, vc) = split(v)
        magusq = ur**2 + uc**2
        magvsq = vr**2 + vc**2
        diff = magusq - magvsq
        return inner(diff, diff)*dx + inner(grad(diff), grad(diff))*dx 

mu_val = start
while mu_val <= mu_end + step:
  info("Considering mu = %s" % mu_val)
  roots[mu_val] = []
  mu.assign(mu_val)
  problem = H1ForwardProblemMagnitude(F, Z, phi, bcs, power=power, shift=shift)
  problem.deflate(zero)

  prev_roots = roots[prev_mu_val]
  out = File("output-%2.4f/roots.pvd" % mu_val)
  

  guess_count = {} # how many times have we tried a particular guess
  guesss = guesses(prev_roots)

  while len(guesss) > 0:
    print('to len to guesss ine %s'%len(guesss))
    guess = guesss.pop(0) # fetch the first guess
    if guess.label in guess_count:
      guess_count[guess.label] += 1
    else:
      guess_count[guess.label] = 0
    j = guess_count[guess.label]

    phi.assign(guess)
    info("  Assigning %s (sub-iteration %s)" % (guess.label, j))

    try:
      solver.solve(problem, phi.vector())

      # Let's check it's an actual root of the original problem.
      res = assemble(F)
      
      [bc.apply(res, phi.vector()) for bc in bcs]
      info("  |F(root)|: %s" % res.norm("l2"))
      assert res.norm("l2") < 1.0e-8

      roots[mu_val].append(phi.copy(deepcopy=True))
      problem.deflate(phi)

      info("  Solution attempt (%s, %s) for mu = %s successful" % (guess.label, j, mu_val))

      distances = []
      for otherroot in roots[mu_val][:-1]:
        d = assemble(inner(phi - otherroot, phi - otherroot)*dx)
        distances.append(d)
      if len(distances) > 0:
        info("Distances to other roots: %s" % distances)

      phiex.assign(roots[mu_val][-1])
      out << phiex
      File("output-%2.4f/root-%d.xml.gz" % (mu_val, len(roots[mu_val])-1)) << roots[mu_val][-1]
      j += 1

      if len(roots[mu_val]) == nroots(mu_val):
        break

      guesss.append(guess)

    except ValueError:
      info("  Solution attempt (%s, %s) for mu = %s unsuccessful" % (guess.label, j, mu_val))

    if len(roots[mu_val]) == nroots(mu_val):
      break

  if mu_val in roots:
    info("mu = %s: %d roots found." % (mu_val, len(roots[mu_val])))
    info("Memory used: %s" % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    del roots[prev_mu_val]
    prev_mu_val = mu_val
    mu_val += step

  gc.collect()

print 'type res is'
print type(res)