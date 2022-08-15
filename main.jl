include("src/SSE.jl")
using .SSE

## Construct the lattice (i.e., the geometry of interaction terms, or `bond_map`)
Lx = 16
Ly = 16
bond_map = square_lattice_NN(Lx=Lx, Ly=Ly)
# println("bond_map: "); display(bond_map); println("\n");


## Construct interactions (i.e. associate a local Hamiltonian to every bond on the lattice)
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
inter = Interaction(H, bond_map)
# println(inter.dof_max)
# println(inter.n_bonds)
# println(inter.n_dofs)
# println("bond_map: "); display(inter.bond_map); println("\n");
# println(inter.hams, "\n")
# println(inter.hams_bonds)


## Check compatibility with directed loops update scheme
# println("Compatible with the directed loop update scheme: ", is_compatible_with_dir_loops(inter), "\n")


## Solve directed loop equations
P = solve_directed_loop_equations(inter)

##
beta = 16  # inverse temperature

# First, we perform `steps_eq` MC sweeps to equilibrate the system. Then we perform `nbins` bins with `steps_bin` MC sweeps in each, and we do measurements after each bin.
nbins = 100
steps_bin = 200
steps_eq = 5000

## Initialize values
# vals = Vals(5, 5, [1,1,2,1,2,1], Int8[2,1,2,1,1], [2,1,2,4,4], Int8[1,1,1,1,1], [4, 2, 3, 2, 2], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
# vals = Vals(5, [1,2,1,2], 2, 4)
# vals = Vals(5, 4, 2, 4)
vals = Vals(M=10, inter=inter)

## Equilibration
open("Mn_saturation.dat", "w") do f
    for i in 1:steps_eq
        diagonal_update!(vals=vals, inter=inter, beta=beta)
        linked_vertex_list!(vals=vals, inter=inter)
        off_diagonal_update!(; vals=vals, inter=inter, mode="directed-loops", P=P, loop_flips=10, new_state_of_entr_leg="shifted", shifted_by=1, max_loop_length=200)
        adjust_cut_off!(vals=vals, inter=inter)
        println(f, vals.n, "\t", vals.M)
        flush(f)
    end
end







