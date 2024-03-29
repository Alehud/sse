export Interaction, Hamiltonian


"""
Hamiltonian (one of the terms in the model's Hamiltonian).
# Fields
- `dof_max::Int32`: each d.o.f. is represented by an integer from 1 to `dof_max`
- `n_legs::Int32`: number of legs from one side of the vertex (i.e., on how many degrees of freedom the term acts)
- `ham::Array{<:Number, 2}`: Hamiltonian matrix
- `labels::Vector{CartesianIndex{2}}`: cartesian coordinates of nonzero elements of `ham`; the i'th nonzero element can be accessed as `ham[labels[i][1], labels[i][2]]` or `ham[labels[i]]`.
- `labels_dof::Vector{Tuple{Tuple{Vararg{Int32}}, Tuple{Vararg{Int32}}}}`: values of degrees of freedom at each nonzero matrix element

"""
struct Hamiltonian
    dof_max::Int32
    n_legs::Int32
    ham::Array{<:Number, 2}
    labels::Vector{CartesianIndex{2}}
    labels_dof::Vector{Tuple{Vector{<:Integer}, Vector{<:Integer}}}
    
    function Hamiltonian(dof_max, ham)
        ss = size(ham)
        if length(ss) ≠ 2
            raise(error("`ham` is not a 2D array."))
        elseif ss[1] ≠ ss[2]
            raise(error("`ham` is not a square matrix."))
        elseif !ishermitian(ham)
            raise(error("`ham` is not Hermitian."))
        elseif any(diag(ham) .< 0)
            raise(error("`ham` contains negative diagonal elements. Add an appropriate constant to the Hamiltonian."))
        elseif any(ham .< 0)
            raise(error("`ham` contains negative off-diagonal elements. The probabilities in QMC cannot be properly defined."))
        elseif !isinteger(log(dof_max, ss[1]))
            raise(error("Dimensions of `ham` are incompatible with `dof_max`."))
        else
            n_legs = round(Int, log(dof_max, ss[1]))
            labels = findall(ham .≠ 0)
            labels_dof = Vector{Tuple{Vector{<:Integer}, Vector{<:Integer}}}(undef, length(labels))
            for (i, lbl) in enumerate(labels)
                labels_dof[i] = (int2base(lbl[1] - 1, dof_max, n_legs, offset_zero=true), int2base(lbl[2] - 1, dof_max, n_legs, offset_zero=true))
            end
            new(dof_max, n_legs, ham, labels, labels_dof)
        end
    end
end


"""
List of Hamiltonians (which can be different terms in the model's Hamiltonian) and bonds associated to each of them.
# Fields
- `dof_max::Int32`: each d.o.f. is represented by an integer from 1 to `dof_max`
- `hams::Vector{Array{<:Number, 2}}`: vector containing different types of local Hamiltonians
- `hams_bonds::Vector{Vector{<:Integer}}`: vector of vectors with bonds corresponding to every Hamiltonian from `ham` (has to be the same length as `ham`)
- `bond_map::Array{<:Integer, 2}`: spins associated to each bond (see `construct_lattice()` function from ConstructLattice.jl)
- `n_bonds::Integer`: total number of bonds
- `n_dofs::Integer`: total number of degrees of freedom
- `n_legs_max::Integer`: maximal number of legs (from one side) in a single vertex

"""
struct Interaction
    dof_max::Integer
    hams::Vector{Hamiltonian}
    hams_bonds::Vector{Vector{<:Integer}}
    bond_map::Array{<:Integer, 2}
    n_bonds::Integer
    n_dofs::Integer
    n_legs_max::Integer
    # bipartition_map::Vector{<:Integer}
    
    function Interaction(dof_max::Integer, hams::Vector{Hamiltonian}, hams_bonds::Vector{Vector{<:Integer}}, bond_map::Array{<:Integer, 2})
        if length(hams) ≠ length(hams_bonds)
            raise(error("Length of `hams` is not equal to the length of `hams_bonds`."))
        else
            for (H, idx) in zip(hams, hams_bonds)
                if H.dof_max ≠ dof_max
                    raise(error("`dof_max` of one of the Hamiltonians from `hams` is not equal to the `dof_max` of the interaction.`"))
                elseif any(bond_map[1:H.n_legs, idx] .== 0) || any(bond_map[(H.n_legs+1):end, idx] .≠ 0)
                    raise(error("Size of `hams` is not compatible with the number of d.o.f. in `bond_map` columns that correspond to `hams_bonds`."))
                end
            end
        end

        n_legs_max = 0
        for H in inter.hams
            if H.n_legs > n_legs_max
                n_legs_max = H.n_legs
            end
        end

        n_dofs = maximum(bond_map)
        n_bonds = size(bond_map)[2]

        # # Assign ±1 to the two sublattices of the bipartite lattice. Assign all zeros if the lattice is not bipartite.
        # bipartition_map = zeros(Int8, n_dofs)
        # bipartition_map[1] = 1
        # flag = false
        # for i in 1:n_bonds
        #     for site in bond_map[2:end, i]
        #         if bipartition_map[site] == 0
        #             bipartition_map[site] = bipartition_map[bond_map[1, i]] * (-1)
        #         elseif bipartition_map[site] == bipartition_map[bond_map[1, i]]
        #             flag = true
        #             break
        #         end
        #     end
        #     if flag
        #         bipartition_map = zeros(Int8, n_dofs)
        #         break
        #     end
        # end

        new(dof_max, hams, hams_bonds, bond_map, n_bonds, n_dofs, n_legs_max)
    end

    function Interaction(H::Hamiltonian, bond_map::Array{<:Integer, 2})
        new(H.dof_max, [H], [[1:size(bond_map)[2];]], bond_map, size(bond_map)[2], maximum(bond_map), H.n_legs)
    end
end
