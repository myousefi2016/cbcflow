
# Function space conventions:
Q = scalar pressure space
U = scalar velocity space
V = vector velocity space
W = mixed velocity*pressure space

# To make functions:
Qc = Q  |  W.sub(1).collapsed()
Uc = U  |  V.sub(0).collapsed()  |  W.sub(0).sub(0).collapsed()
Vc = V  |  W.sub(0).collapsed()

# To make DBCs:
Qr = Q | W.sub(1)
Ur = [U  |  V.sub(d)  |  W.sub(0).sub(d)  for d in dims]
Vr = V  |  W.sub(0)




controls = problem.controls(Uc, Qc)

icu, icp = problem.initial_conditions(Uc, Qc, controls)

bcs = problem.boundary_conditions(Ur, Qr, t, controls)

bcu, bcp = self.transform_bcs_v1(bcs)
bcu, bcp, apbc, Lpbc = self.transform_bcs_v2(bcs)
bcup, Lupbc = self.transform_bcs_v3(bcs)

problem.update_controls(spaces, controls)
problem.update_boundary_conditions(spaces, bcs)




# Creating spaces for problem
spaces = Spaces(mesh, self.params.u_degree, self.params.p_degree, "normal")

# Initial data fetching from problem:
observations = problem.observations(spaces, t)
controls = problem.controls(spaces)
ics = problem.initial_conditions(spaces, controls)
bcs = problem.boundary_conditions(spaces, u, p, t, controls)

# Alternative application of initial condition for different schemes:
assign_ics_split(u0, p0, spaces, ics)
assign_ics_segregated(u0, p0, spaces, ics)
assign_ics_coupled(up0, spaces, ics)

# Alternative transformations for different schemes:
bcu, bcp, apbc, Lpbc = self.bcs_for_split_penalty(bcs)
bcu, bcp, apbc, Lpbc = self.bcs_for_segregated_penalty(bcs)

bcu, bcp = make_split_bcs(problem, spaces, bcs)
bcu, bcp = make_segregated_bcs(problem, spaces, bcs)
bcu, Lbc = make_mixed_bcs(spaces, bcs)

# TODO:
bcu, bcp, apbc, Lpbc = make_split_penalty_bcs(problem, spaces, bcs)
bcu, bcp, apbc, Lpbc = make_segregated_penalty_bcs(problem, spaces, bcs)


# Time loop update of problem data:
problem.update_observations(spaces, t, observations)
problem.update_boundary_conditions(spaces, t, bcs)
