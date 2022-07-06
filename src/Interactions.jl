module Interactions

using LinearAlgebra

export Interaction


"""
Interaction terms for all bonds 
# Fields
- `dof_max::Int32`: each d.o.f. is represented by an integer from 1 to `dof_max`.
- `hams::Vector{Array{<:Number}}`: vector containing different types of local Hamiltonians
- `int_bonds::Vector{Vector{<:Integer}}`: vector of vectors with bonds corresponding to every Hamiltonian from `ham` (has to be the same length as `ham`)
- `bonds::Array{<:Integer}`: spins associated to each bond (see `construct_lattice()` function from ConstructLattice.jl)

"""
struct Interaction
    dof_max::Int32
    hams::Vector{Array{<:Number}}
    int_bonds::Vector{Vector{<:Integer}}
    bonds::Array{<:Integer}
    function Interaction(dof_max, hams, int_bonds, bonds)
        if length(hams) ≠ length(int_bonds)
            raise(error("Length of `hams` is not equal to the length of `int_bonds`."))
        else
            for (ham, idx) in zip(hams, int_bonds)
                if !ishermitian(ham)
                    raise(error("`ham` is not Hermitian."))
                else
                    if any(diag(ham) .< 0)
                        raise(error("`ham` contains negative diagonal elements. Add an appropriate constant to the Hamiltonian."))
                    elseif any(ham .< 0)
                        raise(error("`ham` contains negative off-diagonal elements. The probabilities in QMC cannot be properly defined."))
                    else
                        ss = size(ham)
                        bond_size = round(Int, log(dof_max, ss[1]))
                        if any(bonds[1:bond_size, idx] .== 0) || any(bonds[(bond_size+1):end, idx] .≠ 0)
                            raise(error("Size of `hams` is not compatible with the number of d.o.f. in `bonds` columns corresponding to `int_bonds`."))
                        end
                    end
                end
            end
        end
        new(dof_max, hams, int_bonds, bonds)
    end
end


end