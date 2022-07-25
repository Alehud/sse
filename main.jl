include("src/SSE.jl")
using .SSE

# Construct the lattice (i.e., the geometry of interaction terms, or `bond_map`)
Lx = 2
Ly = 2
bond_map = square_lattice_NN(Lx=Lx, Ly=Ly)
# println("bond_map: "); display(bond_map); println("\n");
Nb = size(bond_map)[2]     # total number of bonds
N = Lx * Ly             # total number of d.o.f.

# Construct interactions (i.e. associate a local Hamiltonian to every bond on the lattice)
dof_max = 2     # each d.o.f. is represented by an integer from 1 to `dof_max`.
# Local Hamiltonian
Δ = 0.5
J = 1.0
h = 0.5
hb = h / (2 * 2J)
ε = ((1-Δ)/2 - hb)/2
H_XXZ = [ε 0 0 0
         0 Δ/2+hb+ε 0.5 0
         0 0.5 Δ/2+hb+ε 0
         0 0 0 2hb+ε]
# println("H_XXZ: "); display(H_XXZ); println("\n");

H = Hamiltonian(dof_max, H_XXZ)
# println(H.ham)
# println(H.dof_max)
# println(H.n_legs)
# println(H.labels)
# println(H.labels_dof)

hams = [H]
ham_bonds = [[1:Nb;]]
inter = Interaction(dof_max, hams, ham_bonds, bond_map)
# println(inter)

# println("Compatible with the directed loop update scheme: ", is_compatible_with_dir_loops(inter))
solve_directed_loop_equations(inter)








# beta = 0.1  # inverse temperature

# # First, we perform 'steps_eq' MC sweeps to equilibrate the system. Then we perform 'nbins' bins with 'steps_bin' MC sweeps in each, and we do measurements after each bin.
# nbins = 100
# steps_bin = 200
# steps_eq = 100

# bond_map = construct_square_lattice(Lx=Lx, Ly=Ly)
# display(bond_map)

# # Initialization
# spins = Array{Bool}(undef, N)
# initialize!(spins)               # randomly initialize the spin state, 0 (false) - spin down, 1 (true) - spin up
# println("\n", spins)
# vals = Vals(10, 0)               # mutable struct Vals is defined in Updates.jl                    
# # vals.M is the highest power in the Taylor expansion, i.e., the number of layers in the "time" direction.
# # vals.n is the number of non-identity operators on the lattice
# op_type = zeros(Int8, vals.M)    # type of operator in each "time" layer, 0 - no operator (identity), 1 - diagonal operator, 2 - off-diagonal operator
# op_ind = fill(-1, vals.M)        # bond index of the operator (taken from `bond_map`), from 1 to Nb, -1 - no operator (identity)
# first_leg = fill(-1, N)          # label of the first (closest to m=1 "time" layer) leg of a non-identity operator on each site
# last_leg = fill(-1, N)           # label of the last (closest to m=M "time" layer) leg of a non-identity operator on each site
# vertex_list = fill(-1, 4vals.M)       # linked vertex list, see comment to linked_vertex_list()











