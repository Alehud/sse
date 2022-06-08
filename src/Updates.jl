module Updates

using Random

export initialize!, Vals, equilibrate!, diagonal_update!, off_diagonal_update, linked_vertex_list, adjust_cut_off, measure

# A mutable struct is used, so that we can change M and n from inside of functions.
mutable struct Vals
    M::Int
    n::Int
end

function initialize!(dofs::Vector{Bool})
    for dof in dofs
        dof = bitrand()
    end
end

function equilibrate!()
    
end

function diagonal_update!(;spins::Vector{Bool}, bonds::Vector{<:Integer}, op_type::Vector{<:Integer}, 
                                            op_ind::Vector{<:Integer}, vals::Vals, Nb::Integer, beta::AbstractFloat)
    # Sweep over M "time" layers
    for m in 1:vals.M
        # Check the operator type
        if op_type[m] == 0
            # If there is no operator, insert one at a random bond with the acceptance probability
            bond = rand(1:Nb)   # choose a random bond
            # Accept the move with the acceptance probability (if prob>1, then accept automatically)
            prob = beta*Nb / (2*(vals.M - vals.n))
            if rand() < prob
                vals.n += 1
            end
        elseif op_type[m] == 1
        else
        end
    end
end

function off_diagonal_update()
end

function linked_vertex_list()
end

function adjust_cut_off()
end

function measure()
end





end