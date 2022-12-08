include("src/SSE.jl")
using .SSE

## Construct the lattice (i.e., the geometry of interaction terms, or `bond_map`)
Lx = 4
Ly = 4
bond_map = square_lattice_NN(Lx=Lx, Ly=Ly)
println("bond_map: "); display(bond_map); println("\n");


## Construct interactions (i.e. associate a local Hamiltonian to every bond on the lattice)
dof_max = 2     # each d.o.f. is represented by an integer from 1 to `dof_max`.
# Local Hamiltonian
Δ = 1.0
J = 1.0
h = 0
hb = h / (2 * 2J)
ε = ((1-Δ)/2 - hb)/2
H_XXZ = [ε 0 0 0
         0 Δ/2+hb+ε 0.5 0
         0 0.5 Δ/2+hb+ε 0
         0 0 0 2hb+ε]
println("H_XXZ: "); display(H_XXZ); println("\n");
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

# First, we perform `steps_eq` MC sweeps to equilibrate the system. Then we perform `nbins` bins with `steps_bin` MC sweeps in each, and we do measurements after each bin.
nbins = 200
steps_bin = 200
steps_eq = 3000

## Temperature sweep
for T in 0.1:0.1:2
    println("T = ", T, " -------------------------------------------------------------------------------------------------")
    beta = 1/T

    ## Initialize values
    vals = Vals(M=10, inter=inter)

    ## Equilibration
    println("Equilibration:")
    open("Mn_saturation_T=$T.dat", "w") do f
        for i in 1:steps_eq
            if i % 500 == 0
                println(i, " / ", steps_eq)
                flush(stdout)
            end
            diagonal_update!(vals=vals, inter=inter, beta=beta)
            linked_vertex_list!(vals=vals, inter=inter)
            off_diagonal_update!(; vals=vals, inter=inter, mode="directed-loops", P=P, loop_flips=10, new_state_of_entr_leg="shifted", shifted_by=1, max_loop_length=200)
            adjust_cut_off!(vals=vals, inter=inter)
            println(f, vals.n, "\t", vals.M)
            flush(f)
        end
    end
    println("Equilibration finished.")

    # ## Main cycle, collect measurements
    # println("Measurement: ")
    # bipartition_mask = Int8[(-1)^(i-1) for i=1:inter.n_dofs]
    # open("measurements_T=$T.dat", "w") do f
    #     for i in 1:nbins
    #         if i % 10 == 0
    #             println("Bin ", i)
    #             flush(stdout)
    #         end
    #         obs = Observables()
    #         for j in 1:steps_bin
    #             diagonal_update!(vals=vals, inter=inter, beta=beta)
    #             linked_vertex_list!(vals=vals, inter=inter)
    #             off_diagonal_update!(; vals=vals, inter=inter, mode="directed-loops", P=P, loop_flips=10, new_state_of_entr_leg="shifted", shifted_by=1, max_loop_length=200)
    #             accumulate_bin_values!(; vals=vals, inter=inter, bipartition_mask=bipartition_mask, obs=obs)
    #         end
    #         measure!(vals=vals, inter=inter, steps_bin=steps_bin, beta=beta, obs=obs)
    #         println(f, obs.E, "\t", obs.C, "\t", obs.mag_stag², "\t", obs.χ)
    #         flush(f)
    #     end
    # end

    ## Main cycle, collect measurements
    println("Measurement C-like: ")
    open("measurements_T=$T.dat", "w") do f
        for i in 1:nbins
            if i % 100 == 0
                println("Bin ", i)
                flush(stdout)
            end
            obs = ObservablesC()
            for j in 1:steps_bin
                diagonal_update!(vals=vals, inter=inter, beta=beta)
                linked_vertex_list!(vals=vals, inter=inter)
                off_diagonal_update!(; vals=vals, inter=inter, mode="directed-loops", P=P, loop_flips=10, new_state_of_entr_leg="shifted", shifted_by=1, max_loop_length=200)
                measureC!(vals=vals, inter=inter, steps_bin=steps_bin, beta=beta, obs=obs)
            end
            obs.enrg1 = obs.enrg1 / steps_bin
            obs.enrg2 = obs.enrg2 / steps_bin
            obs.amag2 = obs.amag2 / steps_bin
            obs.ususc = obs.ususc / steps_bin
            obs.enrg2 = (obs.enrg2 - obs.enrg1*(obs.enrg1 + 1.0)) / inter.n_dofs
            obs.enrg1 = -(obs.enrg1 / (beta * inter.n_dofs) - 0.25*inter.n_bonds/inter.n_dofs)
            obs.amag2 = 3.0*obs.amag2 / (inter.n_dofs)^2
            obs.ususc = beta * obs.ususc / inter.n_dofs
            println(f, obs.enrg1, "\t", obs.enrg2, "\t", obs.amag2, "\t", obs.ususc)
            flush(f)
        end
    end
end






