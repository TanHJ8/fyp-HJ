import numpy as np
import os
from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
import ufl

# 1. DEFINE THE DOMAIN
domain = mesh.create_unit_square(MPI.COMM_WORLD, 64, 64, mesh.CellType.quadrilateral)

# 2. DEFINE FUNCTION SPACE
V = fem.functionspace(domain, ("Lagrange", 1))
import numpy as np
# 1. DEFINE THE DOMAIN
domain = mesh.create_unit_square(MPI.COMM_WORLD, 64, 64, mesh.CellType.quadrilateral)

# 2. DEFINE FUNCTION SPACE
V = fem.functionspace(domain, ("Lagrange", 1))
import numpy as np
import os
from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
import ufl

# 1. DEFINE THE DOMAIN
domain = mesh.create_unit_square(MPI.COMM_WORLD, 64, 64, mesh.CellType.quadrilateral)

# 2. DEFINE FUNCTION SPACE
V = fem.functionspace(domain, ("Lagrange", 1))

# 3. PARAMETERS
t = 0
T_final = 2.0
num_steps = 50docker run -ti --rm -v $(pwd):/home/shared -w /home/shared dolfinx/dolfinx:v0.9.
0~
dt = T_final / num_steps
rho = 1.0; c = 1.0; k = 0.1 

# 4. BOUNDARY CONDITIONS
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)
boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)

u_D = fem.Function(V)
u_D.x.array[:] = 1.0
bc = fem.dirichletbc(u_D, boundary_dofs)

# 5. INITIAL CONDITION
T_n = fem.Function(V)
T_n.name = "T_n"
T_n.x.array[:] = 0.0

# 6. VARIATIONAL PROBLEM
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
F = rho * c * (u - T_n) / dt * v * ufl.dx + k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
a = ufl.lhs(F)
L = ufl.rhs(F)

problem = LinearProblem(a, L, bcs=[bc], 
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})

# 7. TIME LOOP
print("Starting simulation...")
with io.XDMFFile(domain.comm, "heat_simulation.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    
    for i in range(num_steps):
        t += dt
        
        # The Correct Solver Line
        T_h = problem.solve()
        
        # Write to file
        xdmf.write_function(T_h, t)
        
        # Update for next step
        T_n.x.array[:] = T_h.x.array
        
        if (i+1) % 10 == 0:
            print(f"Time step {i+1}/{num_steps} completed (t={t:.2f})")

print("Done! Check the Files folder.")
# 3. PARAMETERS
t = 0
T_final = 2.0
num_steps = 50docker run -ti --rm -v $(pwd):/home/shared -w /home/shared dolfinx/dolfinx:v0.9.
0~
dt = T_final / num_steps
rho = 1.0; c = 1.0; k = 0.1 

# 4. BOUNDARY CONDITIONS
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)
boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)

u_D = fem.Function(V)
u_D.x.array[:] = 1.0
bc = fem.dirichletbc(u_D, boundary_dofs)

# 5. INITIAL CONDITION
T_n = fem.Function(V)
T_n.name = "T_n"
T_n.x.array[:] = 0.0

# 6. VARIATIONAL PROBLEM
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
F = rho * c * (u - T_n) / dt * v * ufl.dx + k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
a = ufl.lhs(F)
L = ufl.rhs(F)

problem = LinearProblem(a, L, bcs=[bc], 
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})

# 7. TIME LOOP
print("Starting simulation...")
with io.XDMFFile(domain.comm, "heat_simulation.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    
    for i in range(num_steps):
        t += dt
        
        # The Correct Solver Line
        T_h = problem.solve()
        
        # Write to file
        xdmf.write_function(T_h, t)
        
        # Update for next step
        T_n.x.array[:] = T_h.x.array
        
        if (i+1) % 10 == 0:
            print(f"Time step {i+1}/{num_steps} completed (t={t:.2f})")

print("Done! Check the Files folder.")
