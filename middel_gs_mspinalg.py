"""
Created on Fri Jul 14 14:15:33 2017

This script uses the MSPIN algorithm with MSPIN approximation 1 for the Jacobian (see thesis).
The components involve solving the subproblems as standalond and sequentially in order to obtain
the solutions to the subproblem F_preconditioned = (g,h)^T. 
@author: Georgios
"""
from dolfin import *
from petsc4py import PETSc
import sys; sys.path.append('..')
from deflation import ForwardProblem
import gc
gc.disable()
import resource
import glob
from utils import prevent_MPI_output
import copy
import numpy as np

# original tolerances
 #                      --petsc.snes_atol 1.0e-5
  #                     --petsc.ksp_atol 1.0e-12
  #                     --petsc.ksp_rtol 1.0e-10 
  # natol = 1.0e-8
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
mu_end = 0.600
natol = 1.0e-8

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

    def normsq(self, u, v):
        (ur, uc) = split(u) 
        (vr, vc) = split(v)
        magusq = ur**2 + uc**2                                                                           
        magvsq = vr**2 + vc**2                                                                           
        diff = magusq - magvsq
        return inner(diff, diff)*dx + inner(grad(diff), grad(diff))*dx                           

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
     
     def norm(F, u, bcs):
         b = assemble(F)
         [bc.apply(b, u.vector()) for bc in bcs]
         return b.norm("l2")

     def printnorm(i, n):
         print "%3d SNES Function norm %1.15e" % (i, n)

     du = Function(u.function_space())
     J = derivative(F, u, du)
     hbcs = [DirichletBC(bc) for bc in bcs]; [hbc.homogenize() for hbc in hbcs]
     i = 0
     maxit = 100
     n = norm(F, u, bcs)
     printnorm(i, n)

     while n > atol and i < maxit:
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

     if i == maxit: raise ValueError
     return n 

def MSPIN_alg_sep(Ffirst, Fsecond, F_old, phi,gh,bcs, atol=1.0e-8, deflation = None):#,deflation2 = None):
    
    [bc.apply(phi.vector()) for bc in bcs]
    
    
    def norm(F, u, bcs):
        b = assemble(F)
        [bc.apply(b, u.vector()) for bc in bcs]
        return b.norm("l2")

    def printnorm(i, n):
        print "%3d SNES Function norm %1.15e" % (i, n)
   
    i = 0
    maxit = 100
    hbcs = [DirichletBC(bc) for bc in bcs]; [hbc.homogenize() for hbc in hbcs]
    n = norm(F_old, phi, bcs)
    #n0 = norm(F_old, phi_old, bcs)
    #print('this is the printed norm before solution')
    printnorm(i, n)
    #phi_old.assign(rc_0)
    #phi_mspin.assign(rc_0)
    # options
    linesearch = True   # using linesearch
    ls_monitor = True   # linesearch monitor
    snes_monitor = True # nonlinear solver monitor
    plot_path = True    # plot paths of approximated solutions
    
    dphi = TrialFunction(phi.function_space())
    J = derivative(F_old, phi, dphi)
    phi_test = TestFunction(phi.function_space())
    dummy = inner(as_vector([1.0, 1.0]), phi_test)*dx
    out = Function(phi.function_space(), name="out")
    #pv = Function(phi_old.function_space())
    phi_tmp= Function(phi.function_space())
        #==================================

   # print('======================================================\n')
   # print('                      MSPIN                           \n')
    #print('======================================================\n')
    #random_n = Function(phi_old.function_space(), name="Noise")
    #eps = 1.0e-4
    #random_n.vector()[:] = eps*np.random.random(random_n.vector().size())
    #solve(F_old==0,phi_old,hbcs)
    #newton(F_old, phi_old, bcs, natol,d) 
    #[bc.apply(phi_old.vector()) for bc in bcs]
    #phi_mspin.assign(phi_old)
    
    #pvd = File("mspin_solution.pvd")
    #pvd << random_n
    #phi_mspin.assign(gh_init)
    #pvd = File("mspin_solution.pvd")
    #gh.rename("gh", "gh")

        # gh_0 =project(Expression(("0.0","0.0"), de
    while n > atol and i < maxit:#n> atol and i < maxit:
        #out.vector().zero()
        #if i > 1:
         #   break
        #[bc.apply(phi_old.vector()) for bc in bcs]
        #phi_mspin.assign(phi_old+random_n)
        # gh_0 =project(Expression(("0.0","0.0"), degree = 2), Z) # initial guess for (g,h) which must be zero

        phi_tmp.assign(phi) # to store phi_olduntil I produce the matrix and then restore it              
                       
        #gh.assign(gh_init) # initial gues must always be zero for this one
        #[bc.apply(gh.vector()) for bc in bcs]
    
           # Step 1: obtain the solution to the subproblem
        # this part must actually be solved by INB up to a tolerance
        #solve(F_old==0,phi_old,hbcs)
       # File("check_phiold_mspin.pvd") << phi_old
        #phi_mspin.assign(phi_old)
        
        #solve(F_mspin == 0, gh ,hbcs)
        #problem1 = NonlinearVariationalProblem(Ffirst, gh, hbcs)
        #problem2 = NonlinearVariationalProblem(Fsecond, gh, hbcs)
        #solver1 = NonlinearVariationalSolver(problem1)
        #solver2 = NonlinearVariationalSolver(problem2)
        #prm1 = solver1.parameters
        #prm2 = solver2.parameters
        #prm1["newton_solver"]["absolute_tolerance"] = 1E-3
        #prm1["newton_solver"]["relative_tolerance"] = 1E-3
        #prm2["newton_solver"]["absolute_tolerance"] = 1E-3
        #prm2["newton_solver"]["relative_tolerance"] = 1E-3
        
        #solver1.solve()
        #solver2.solve()
        solve(Ffirst == 0, gh ,hbcs,solver_parameters={"newton_solver": {"absolute_tolerance": 1.0e-3}})
        solve(Fsecond == 0, gh ,hbcs,solver_parameters={"newton_solver": {"absolute_tolerance": 1.0e-3}})
       # pvd << gh
        #solve(F_old == 0, phi_,hbcs)
       # print('check stopping condition on gh whic has norm %s' % (dolfin.norm(gh,'L2')))
       # File("check_gh_mspin.pvd")<<gh

        dofs_0 = Z.sub(0).dofmap().dofs() # dofs associated with subspace 0
        is_0   = PETSc.IS().createGeneral(dofs_0) # a PETSc index set
        dofs_1 = Z.sub(1).dofmap().dofs() # dofs associated with subspace 1
        is_1   = PETSc.IS().createGeneral(dofs_1) # a PETSc index set                   

        gh_v = as_backend_type(gh.vector()).vec()
        gh_0 = gh_v.getSubVector(is_0)
        gh_1 = gh_v.getSubVector(is_1)
        
        pv = Function(phi.function_space())
        pv_v = as_backend_type(pv.vector()).vec()
        pv_0 = pv_v.getSubVector(is_0)     
        pv_1 = pv_v.getSubVector(is_1)
        
        phi_v = as_backend_type(phi.vector()).vec()
        phi_0 = phi_v.getSubVector(is_0)
        phi_1 = phi_v.getSubVector(is_1)
        

        pv_0.axpy(+1.0,phi_0)
        pv_0.axpy(-1.0,gh_0)
        pv_1.axpy(+1.0,phi_1)
        
        pv_v.restoreSubVector(is_0, pv_0)
        pv_v.restoreSubVector(is_1, pv_1)
        
        
        phi.assign(pv)

        A = PETScMatrix()
        asm = SystemAssembler(J, dummy, hbcs)
        asm.assemble(A)



        A = A.mat()
        A_00 = A.getSubMatrix(is_0, is_0) # might need A.getSubMatrix, createSubMatrix
        A_01 = A.getSubMatrix(is_0, is_1) # depends on your version
        A_10 = A.getSubMatrix(is_1, is_0)
        A_11 = A.getSubMatrix(is_1, is_1)

        # Compute L = [A_00 0 ; A_10 A_11] . gh --- I don't know if this is the
        # right thing to compute or not, but it shows how to use PETSc

    
        Lgh = gh_v.duplicate()
        Lgh_0 = Lgh.getSubVector(is_0)
        Lgh_1 = Lgh.getSubVector(is_1)
#
        A_00.mult(gh_0, Lgh_0)
#        tmp = Lgh_0.duplicate()
#        A_01.mult(gh_1, tmp)
#        Lgh_0.axpy(+1.0, tmp)
        A_10.mult(gh_0, Lgh_1)
        tmp = Lgh_1.duplicate()
        A_11.mult(gh_1, tmp)
        Lgh_1.axpy(+1.0, tmp)

        Lgh.restoreSubVector(is_0, Lgh_0)
        Lgh.restoreSubVector(is_1, Lgh_1)
        gh_v.restoreSubVector(is_0, gh_0)
        gh_v.restoreSubVector(is_1, gh_1)

        Lgh.scale(-1.0)
        # Step 2
        # OK. Now that we have the RHS, set up a linear system and solve it.
        # Note from GEORGE: I have taken these lines out of the while loop: don't know 
        # whtehre this is correct: I only left the solve part in the while lloop
        #out = Function(Z)
        out_v = as_backend_type(out.vector()).vec()
        # now restore phi_old.
        
        #

        ksp = PETSc.KSP().create() # object that solves linear systems
        ksp.setOperators(A) # matrix will be A
        ksp.setType("preonly")
        ksp.pc.setType("lu")   # use LU factorisation, by default
        ksp.setOptionsPrefix("mspin_")
        ksp.setFromOptions()
        
        # Solve the linear system

        ksp.solve(Lgh, out_v)
       # pvd << out
        
        linesearch = False
        # Step 3 : the linesearch
        #if linesearch:
        #    # need preconditioned system
            #f = 0.5 * (dolfin.norm(gh,'L2'))*(dolfin.norm(gh,'L2'))
        #    f = 0.5 * np.linalg.norm(gh_v)*np.linalg.norm(gh_v)
        #    dg = -2.0 * f
        #    stp = 1.0
            #stp = ls_mspin( phi_old, phi_mspin, gh, out, hbcs, dg, f, stp, F_mspin,gh_init)#ls(x_soln, d, dg, f, stp, options, fun_MSPIN)
       
        stp = float(1)
        if deflation is not None:
            Eu = deflation.derivative().inner(out.vector())
            minv = 1.0 / deflation.evaluate()
            tau = (1 + minv*Eu/(1 - minv*Eu))
            out.assign(tau * out)
        
        #phi_old.assign(phi_oldtmp) # restoring phi_old
        phi.assign(phi_tmp + stp * out)#assign(phi_old + stp*out)
        [bc.apply(phi.vector()) for bc in bcs]
       
        # gh_0 =project(Expression(("0.0","0.0"), de
        i = i + 1
       
        n= norm(F_old,phi,bcs)
        printnorm(i, n)
    
    if i == maxit: raise ValueError
    return n   

 
def forward_MSPIN_total2(Z):
  phi = Function(Z)
  gh  = Function(Z)
  (r, c) = split(phi)
  (g, h) = split(gh)
  (r_test, c_test) = split(TestFunction(Z))

  mesh = Z.mesh()
  x = SpatialCoordinate(mesh)
 
  mag = r**2 + c**2
  mag1 = ((r-g)**2 + c**2)
  mag2 = ((r-g)**2 + (c-h)**2)
  F_old = (
        mu*inner(r, r_test)*dx
      - 0.5*inner(grad(r), grad(r_test))*dx
      - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner(r, r_test)*dx
      - mag*inner(r, r_test)*dx

      + mu*inner(c, c_test)*dx
      - 0.5*inner(grad(c), grad(c_test))*dx
      - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner(c, c_test)*dx
      - mag*inner(c, c_test)*dx
      )
  F1 = (
        mu*inner((r-g), r_test)*dx
        - 0.5*inner(grad((r-g)), grad(r_test))*dx
        - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner((r-g), r_test)*dx
        - mag1*inner((r-g), r_test)*dx
       )
  F2 = inner(h,c_test)*dx
  F3 = (mu*inner((c-h), c_test)*dx
        - 0.5*inner(grad(c-h), grad(c_test))*dx
        - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner((c-h), c_test)*dx
        - mag2*inner((c-h), c_test)*dx
        )
        
  
  Ffirst = F1 +F2
  Fsecond = F1 +F3
  
  F_mspin = (
        mu*inner((r-g), r_test)*dx
        - 0.5*inner(grad((r-g)), grad(r_test))*dx
        - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner((r-g), r_test)*dx
        - mag1*inner((r-g), r_test)*dx
    
        +mu*inner((c-h), c_test)*dx
        - 0.5*inner(grad(c-h), grad(c_test))*dx
        - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner((c-h), c_test)*dx
        - mag2*inner((c-h), c_test)*dx
        )
  
  bcs = [DirichletBC(Z, (0.0, 0.0), "on_boundary")]

  return (F_old,F_mspin, Ffirst, Fsecond, phi, gh,bcs)#(F, phi, gh, bcs)    
  
def forward_MSPIN_total(Z):
  phi = Function(Z)
  gh  = Function(Z)
  (r, c) = split(phi)
  (g, h) = split(gh)
  (r_test, c_test) = split(TestFunction(Z))

  mesh = Z.mesh()
  x = SpatialCoordinate(mesh)
 
  mag = r**2 + c**2
  mag1 = ((r-g)**2 + c**2)
  mag2 = ((r-g)**2 + (c-h)**2)
  F_old = (
        mu*inner(r, r_test)*dx
      - 0.5*inner(grad(r), grad(r_test))*dx
      - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner(r, r_test)*dx
      - mag*inner(r, r_test)*dx

      + mu*inner(c, c_test)*dx
      - 0.5*inner(grad(c), grad(c_test))*dx
      - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner(c, c_test)*dx
      - mag*inner(c, c_test)*dx
      )
  
  F_mspin = (
        mu*inner((r-g), r_test)*dx
        - 0.5*inner(grad((r-g)), grad(r_test))*dx
        - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner((r-g), r_test)*dx
        - mag1*inner((r-g), r_test)*dx
    
        +mu*inner((c-h), c_test)*dx
        - 0.5*inner(grad(c-h), grad(c_test))*dx
        - 0.5*omega**2 * (x[0]**2 + x[1]**2) * inner((c-h), c_test)*dx
        - mag2*inner((c-h), c_test)*dx
        )
  
  bcs = [DirichletBC(Z, (0.0, 0.0), "on_boundary")]

  return (F_old,F_mspin, phi, gh,bcs)#(F, phi, gh, bcs)


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
#(F, phi, bcs) = forward(Z)
#(F_old, F_mspin, phi, gh, bcs) = forward_MSPIN_total(Z)
(F_old,F_mspin, Ffirst, Fsecond, phi, gh,bcs) = forward_MSPIN_total2(Z)
power = 2.0
shift = 1.0

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
  for i in range(len(glob.glob("outputgs-%2.4f/root-*.xml.gz" % prev_mu_val))):
      roots[prev_mu_val].append(Function(Z, "outputgs-%2.4f/root-%d.xml.gz" % (prev_mu_val, i)))

  start = prev_mu_val + step

def nroots(mu):
    if mu <= 0.6:
        return 6
    return float("inf")

#class H1ForwardProblemMagnitude(ForwardProblem):
 #   def norm(self, u, v):
  #      (ur, uc) = split(u)
   #     (vr, vc) = split(v)
    #    magusq = ur**2 + uc**2
     #   magvsq = vr**2 + vc**2
      #  diff = magusq - magvsq
       # return inner(diff, diff)*dx + inner(grad(diff), grad(diff))*dx 

mu_val = start
correct_mu = []
while mu_val <= mu_end + step:
  info("Considering mu = %s" % mu_val)
  roots[mu_val] = []
  mu.assign(mu_val)
  d = DeflationOperator(phi,power,shift)
  d.deflate(zero)

  prev_roots = roots[prev_mu_val]
  out = File("outputgs-%2.4f/roots.pvd" % mu_val)
  

  guess_count = {} # how many times have we tried a particular guess
  guesss = guesses(prev_roots)
  m=0
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
      #solver.solve(problem, phi.vector())
     # if m == 0:
      #    n=newton(F, phi, bcs, natol) 
    #  else:
      #n=newton(F, phi, bcs, natol,d)
      n = MSPIN_alg_sep(Ffirst, Fsecond, F_old, phi,gh,bcs, natol, d)
      #n=MSPIN_alg_total(F_mspin, F_old, phi, gh,  bcs, natol, d)

      #Let's check it's an actual root of the original problem.
      res = assemble(F_old)
      [bc.apply(res, phi.vector()) for bc in bcs]
      info("  |F(root)|: %s" % res.norm("l2"))
      assert res.norm("l2") < 1.0e-8 #1.0e-8
      m+=1
      roots[mu_val].append(phi.copy(deepcopy=True))
      d.deflate(phi)
      if mu_val in roots.keys():
          correct_mu.append(mu_val)

      info("  Solution attempt (%s, %s) for mu = %s successful" % (guess.label, j, mu_val))

      distances = []
      for otherroot in roots[mu_val][:-1]:
        dist = assemble(inner(phi - otherroot, phi - otherroot)*dx)
        distances.append(dist)
      if len(distances) > 0:
        info("Distances to other roots: %s" % distances)

      phiex.assign(roots[mu_val][-1])
      out << phiex
      File("outputgs-%2.4f/root-%d.xml.gz" % (mu_val, len(roots[mu_val])-1)) << roots[mu_val][-1]
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
print correct_mu
