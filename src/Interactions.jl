export Interaction, Hamiltonian


"""
Hamiltonian (one of the terms in the model's Hamiltonian).
# Fields
- `dof_max::Int32`: each d.o.f. is represented by an integer from 1 to `dof_max`
- `n_legs::Int32`: number of legs from one side of the vertex (i.e., on how many degrees of freedom the term acts)
- `ham::Array{<:Number}`: Hamiltonian matrix
- `labels::Vector{CartesianIndex{2}}`: cartesian coordinates of nonzero elements of `ham`; the i'th nonzero element can be accessed as `ham[labels[i][1], labels[i][2]]`
- `labels_dof::Vector{Tuple{Tuple{Vararg{Int32}}, Tuple{Vararg{Int32}}}}`: values of degrees of freedom at each nonzero matrix element

"""
struct Hamiltonian
    dof_max::Int32
    n_legs::Int32
    ham::Array{<:Number}
    labels::Vector{CartesianIndex{2}}
    labels_dof::Vector{Tuple{Tuple{Vararg{Int32}}, Tuple{Vararg{Int32}}}}
    
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
            labels_dof = Vector{Tuple{NTuple{n_legs, Int32}, NTuple{n_legs, Int32}}}(undef, length(labels))
            for (i, lbl) in enumerate(labels)
                labels_dof[i] = (Tuple(int2base(lbl[1] - 1, dof_max, n_legs)), Tuple(int2base(lbl[2] - 1, dof_max, n_legs)))
            end
            new(dof_max, n_legs, ham, labels, labels_dof)
        end
    end
end


"""
List of Hamiltonians (which can be different terms in the model's Hamiltonian) and bonds associated to each of them.
# Fields
- `dof_max::Int32`: each d.o.f. is represented by an integer from 1 to `dof_max`
- `hams::Vector{Array{<:Number}}`: vector containing different types of local Hamiltonians
- `ham_bonds::Vector{Vector{<:Integer}}`: vector of vectors with bonds corresponding to every Hamiltonian from `ham` (has to be the same length as `ham`)
- `bond_map::Array{<:Integer}`: spins associated to each bond (see `construct_lattice()` function from ConstructLattice.jl)

"""
struct Interaction
    dof_max::Int32
    hams::Vector{Hamiltonian}
    ham_bonds::Vector{Vector{<:Integer}}
    bond_map::Array{<:Integer}
    
    function Interaction(dof_max, hams, ham_bonds, bond_map)
        if length(hams) ≠ length(ham_bonds)
            raise(error("Length of `hams` is not equal to the length of `ham_bonds`."))
        else
            for (ham, idx) in zip(hams, ham_bonds)
                if ham.dof_max ≠ dof_max
                    raise(error("`dof_max` of one of the Hamiltonians from `hams` is not equal to the `dof_max` of the interaction.`"))
                elseif any(bond_map[1:ham.n_legs, idx] .== 0) || any(bond_map[(ham.n_legs+1):end, idx] .≠ 0)
                    raise(error("Size of `hams` is not compatible with the number of d.o.f. in `bond_map` columns that correspond to `ham_bonds`."))
                end
            end
        end
        new(dof_max, hams, ham_bonds, bond_map)
    end
end
