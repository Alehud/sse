export Vals, diagonal_update!, off_diagonal_update!, linked_vertex_list!, adjust_cut_off!, measure

# A mutable struct to keep all values that can change during the updates
@with_kw mutable struct Vals
    M::Integer                           # the highest power in the Taylor expansion, i.e., the number of layers in the "time" direction.
    n::Integer                           # the number of non-identity operators on the lattice
    dofs::Vector{<:Integer}              # current state of the degrees of freedom
    op_type::Vector{Int8}                # type of operator in each "time" layer, 0 - no operator (identity), 1 - diagonal operator, 2 - off-diagonal operator
    op_bond::Vector{<:Integer}           # bond index of the operator (taken from `bond_map`), from 1 to `n_bonds`, 0 - no operator (identity)
    op_ham_idx::Vector{Int8}             # Hamilonian index of the operator, i.e., which bond type it belongs to (e.g., if there is a single Hamiltonian, then all elements are 1), 0 - no operator                                 
    op_lbl::Vector{<:Integer}            # linear index (label) of the operator in each time layer, 0 - no operator
    first_leg::Vector{<:Integer}         # label of the first (closest to m=1 "time" layer) leg of a non-identity operator on each site
    last_leg::Vector{<:Integer}          # label of the last (closest to m=M "time" layer) leg of a non-identity operator on each site
    vertex_list::Vector{<:Integer}       # linked vertex list

    function Vals(M::Integer, n::Integer, dofs::Vector{<:Integer}, op_type::Vector{Int8}, op_bond::Vector{<:Integer}, op_ham_idx::Vector{Int8}, op_lbl::Vector{<:Integer}, 
                  first_leg::Vector{<:Integer}, last_leg::Vector{<:Integer}, vertex_list::Vector{<:Integer})
        # TODO: add consistency checks
        if M < 1
            raise(error("`M` has to be ≥ 1."))
        elseif n < 0
            raise(error("`n` has to be ≥ 0."))
        elseif n > M
            raise(error("`n` has to be ≤ M."))
        elseif count(op_type .≠ 0) ≠ n
            raise(error("The number of nonzero elements in `op_type` has to be `n`."))
        elseif count(op_bond .≠ 0) ≠ n
            raise(error("The number of nonzero elements in `op_bond` has to be `n`."))
        elseif count(op_ham_idx .≠ 0) ≠ n
            raise(error("The number of nonzero elements in `op_ham_idx` has to be `n`."))
        elseif count(op_lbl .≠ 0) ≠ n
            raise(error("The number of nonzero elements in `op_lbl` has to be `n`."))
        elseif any(dofs .< 1)
            raise(error("`dofs` have to be ≥ 1."))
        elseif length(op_type) ≠ M
            raise(error("Length of `op_type` has to be equal to `M`."))
        elseif length(op_bond) ≠ M
            raise(error("Length of `op_bond` has to be equal to `M`."))
        elseif length(op_ham_idx) ≠ M
            raise(error("Length of `op_ham_idx` has to be equal to `M`."))
        elseif length(op_lbl) ≠ M
            raise(error("Length of `op_lbl` has to be equal to `M`."))
        elseif !all((op_type .== 0) .|| (op_type .== 1) .|| (op_type .== 2))
            raise(error("`op_type` can only be 0, 1, or 2."))
        elseif any(op_bond .< 0)
            raise(error("`op_bond` have to be ≥ 0."))
        elseif any(first_leg .< 0)
            raise(error("`first_leg` have to be ≥ 0."))
        elseif length(first_leg) ≠ length(dofs)
            raise(error("Length of `first_leg` has to be equal to length of `dofs`."))
        elseif any(last_leg .< 0)
            raise(error("`last_leg` have to be ≥ 0."))
        elseif length(last_leg) ≠ length(dofs)
            raise(error("Length of `last_leg` has to be equal to length of `dofs`."))
        elseif any(vertex_list .< 0)
            raise(error("`vertex_list` have to be ≥ 0."))
        end
        for v in vertex_list
            if (v == 0) || (vertex_list[vertex_list[v]] == v)
                continue
            else
                raise(error("`vertex_list` is constructed incorrectly."))
            end
        end
        new(M, n, dofs, op_type, op_bond, op_ham_idx, op_lbl, first_leg, last_leg, vertex_list)
    end
end

function Vals(M::Integer, dofs::Vector{<:Integer}, dof_max::Integer, n_legs_max::Integer)
    if dof_max < 2
        raise(error("`dof_max` has to be ≥ 2"))
    elseif any(dofs .> dof_max)
        raise(error("`dofs` cannot be larger than `dof_max`."))
    elseif (n_legs_max < 2) || (n_legs_max % 2 ≠ 0)
        raise(error("`n_legs_max` has to be ≥ 2 and even."))
    end
    return Vals(M, 0, dofs, zeros(Int8, M), zeros(Int, M), zeros(Int8, M), zeros(Int, M), zeros(Int, length(dofs)), zeros(Int, length(dofs)), zeros(Int, n_legs_max*M))
end

function Vals(M::Integer, n_dofs::Integer, dof_max::Integer, n_legs_max::Integer)
    if n_dofs < 1
        raise(error("`n_dofs` has to be ≥ 1"))
    end
    return Vals(M, rand(1:dof_max, n_dofs), dof_max, n_legs_max)
end

function Vals(; M::Integer, inter::Interaction)
    return Vals(M, inter.n_dofs, inter.dof_max, 2*inter.n_legs_max)
end

function diagonal_update!(;vals::Vals, inter::Interaction, beta::Real)
    # Sweep over M "time" layers
    for m in 1:vals.M
        # Check the operator type
        if vals.op_type[m] == 0
            # If there is no operator, insert one at a random bond with the acceptance probability
            bond = rand(1:inter.n_bonds)   # choose a random bond
            # Determine which Hamiltonian this bond belongs to
            ham_idx = 0
            for (i, ham_bonds) in enumerate(inter.hams_bonds)
                if bond in ham_bonds
                    ham_idx = i
                    break
                end
            end
            state = vals.dofs[inter.bond_map[:,bond]]   # state of the dofs on the chosen bond
            # Check if the diagonal term is non-zero
            if (state, state) in inter.hams[ham_idx].labels_dof
                # Accept the move with the acceptance probability (if prob>1, then accept automatically)
                accept_prob = beta*inter.n_bonds / (2*(vals.M - vals.n))
                if rand() < accept_prob
                    vals.n += 1
                    vals.op_type[m] = 1
                    vals.op_bond[m] = bond
                    vals.op_ham_idx[m] = ham_idx
                    vals.op_lbl[m] = findfirst(x -> x == (state, state), inter.hams[ham_idx].labels_dof)
                end
            end
        elseif vals.op_type[m] == 1
            # If there is a diagonal operator, remove it.
            # Accept the move with the acceptance probability (if prob>1, then accept automatically)
            accept_prob = 2*(vals.M - vals.n + 1) / (beta * inter.n_bonds)
            if rand() < accept_prob
                vals.n -= 1
                vals.op_type[m] = 0
                vals.op_bond[m] = 0
                vals.op_ham_idx[m] = 0
                vals.op_lbl[m] = 0
            end
        elseif vals.op_type[m] == 2
            # If there is an off-diagonal operator, advance the dofs state
            vals.dofs[inter.bond_map[:,vals.op_bond[m]]] = inter.hams[vals.op_ham_idx[m]].labels_dof[vals.op_lbl[m]][2]
        else
            raise(error("Invalid value of `op_type`. Only 0 (no operator), 1 (diagonal operator), 2 (off-diagonal operator) are allowed."))
        end
    end
end

# TODO: check for interactions with several types of bonds
function linked_vertex_list!(; vals::Vals, inter::Interaction)
    vals.first_leg .= 0
    vals.last_leg .= 0
    vals.vertex_list .= 0
    for m in 1:vals.M
        if vals.op_type[m] == 0
            continue
        end
        v0 = 2*inter.n_legs_max*(m-1) + 1
        ss = inter.bond_map[:,vals.op_bond[m]]
        ss = ss[ss .≠ 0]
        vv = vals.last_leg[ss]
        for (i, (v, s)) in enumerate(zip(vv, ss))
            if v ≠ 0
                vals.vertex_list[v] = v0 + (i-1)
                vals.vertex_list[v0 + (i-1)] = v
            else
                vals.first_leg[s] = v0 + (i-1)
            end
            vals.last_leg[s] = v0 + inter.n_legs_max + (i-1)
        end
    end

    # Linking across periodic "time" boundary
    for s in 1:inter.n_dofs
        f = vals.first_leg[s]
        if f ≠ 0
            l = vals.last_leg[s]
            vals.vertex_list[f] = l
            vals.vertex_list[l] = f
        end
    end
end

function off_diagonal_update!(; vals::Vals, inter::Interaction, mode::String, P::Vector{Array{Float32, 5}}, loop_flips::Integer=1, new_state_of_entr_leg::String="random", shifted_by::Integer=0, max_loop_length::Integer=100)
    if mode ≠ "directed-loops"
        raise(error("Wrong mode of the off-diagonal update."))
    end
    for _ in 1:loop_flips
        # If there are no operators (i.e., n == 0)
        if vals.n == 0
            break
        end

        # Randomly choose initial leg
        v0 = 0
        while true
            v0 = rand(1:vals.M*2*inter.n_legs_max)
            if vals.vertex_list[v0] ≠ 0
                break
            end
        end

        # Choose a new state of the initial leg
        if new_state_of_entr_leg == "random"
            s0 = rand(1:inter.dof_max)
        elseif new_state_of_entr_leg == "shifted"
            m = (v0-1) ÷ (2*inter.n_legs_max) + 1                                                   # "time" layer
            label_dof = deepcopy(inter.hams[vals.op_ham_idx[m]].labels_dof[vals.op_lbl[m]])         # label_dof of the operator (i.e., in the form ([1,2,...], [2,1,...]))
            legₑ = (v0-1) % (2*inter.n_legs_max) + 1                                                # leg of the initial vₑ
            leg_state = label_dof[(legₑ-1)÷inter.n_legs_max + 1][(legₑ-1)%inter.n_legs_max + 1]     # state of the initial leg
            s0 = mod(leg_state - 1 + shifted_by, inter.dof_max) + 1
        else
            raise(error("`new_state_of_entr_leg` can be only \"random\" or \"shifted\"."))
        end

        # Set initial values
        vₑ = v0
        sₑ = s0

        # Save initial configurations. In case the loop becomes too long and is terminated, we will need to restore the original configuration of `vals`.
        dofs_backup = deepcopy(vals.dofs)
        op_lbl_backup = deepcopy(vals.op_lbl)
        op_type_backup = deepcopy(vals.op_type)
        
        # Construct the loop probabilistically according to P, follow it and flip it accordingly along the way.
        loop_length = 0
        while true
            # Get the exit leg and its state according to the probabilities from tensor P
            # P[ham index][label, entrance leg, new state of entrance leg, exit leg, new state of exit leg]
            probs = P[vals.op_ham_idx[m]][vals.op_lbl[m], legₑ, sₑ, :, :]
            r = rand()
            cum = Float32(0.0)
            legₓ = 0
            sₓ = 0
            for leg in 1:size(probs, 1), s in 1:size(probs, 2)
                if r - cum < probs[leg, s]
                    legₓ = leg
                    sₓ = s
                    break
                end
                cum += probs[leg, s]
            end

            # v of the exit leg
            vₓ = (m-1)*(2*inter.n_legs_max) + legₓ

            # Assign new states on entrance and exit legs to the label_dof
            label_dof[(legₑ-1)÷inter.n_legs_max + 1][(legₑ-1)%inter.n_legs_max + 1] = sₑ
            label_dof[(legₓ-1)÷inter.n_legs_max + 1][(legₓ-1)%inter.n_legs_max + 1] = sₓ

            # Update vals
            vals.op_lbl[m] = findfirst(==(label_dof), inter.hams[vals.op_ham_idx[m]].labels_dof)
            if label_dof[1] == label_dof[2]
                vals.op_type[m] = 1
            else
                vals.op_type[m] = 2
            end
            # If vₑ or/and vₓ are among the first legs, then update the stored dofs configuration.
            site = findfirst(==(vₑ), vals.first_leg)
            if site ≠ nothing
                vals.dofs[site] = sₑ
            end
            site = findfirst(==(vₓ), vals.first_leg)
            if site ≠ nothing
                vals.dofs[site] = sₓ
            end

            # If the loop is too long, terminate and restore the original configuration
            loop_length += 1
            if loop_length > max_loop_length
                vals.dofs = deepcopy(dofs_backup)
                vals.op_lbl = deepcopy(op_lbl_backup)
                vals.op_type = deepcopy(op_type_backup)
                break
            end
            # println(loop_length, ".  vₑ: ", vₑ, "   vₓ: ", vₓ)

            # Go to the next vertex
            vₑ = vals.vertex_list[vₓ]
            sₑ = sₓ

            # If the loop closes, terminate.
            if (vₑ == v0) && (sₑ == s0)
                break
            end

            # Otherwise, determine parameters of the new vertex and entrance leg
            m = (vₑ-1) ÷ (2*inter.n_legs_max) + 1                                                    # "time" layer
            label_dof = deepcopy(inter.hams[vals.op_ham_idx[m]].labels_dof[vals.op_lbl[m]])          # label_dof of the operator (i.e., in the form ([1,2,...], [2,1,...]))
            legₑ = (vₑ-1) % (2*inter.n_legs_max) + 1                                                 # leg of the chosen vₑ, i.e., entrance leg of the vertex
        end
    end
    
    # Randomly change dofs that don't take part in any operator (i.e., their worldlines just wind around the p.b.c.).
    for (site, fl) in enumerate(vals.first_leg)
        if fl == 0
            vals.dofs[site] = rand(1:inter.dof_max)
        end
    end
end

function adjust_cut_off!(; vals::Vals, inter::Interaction)
    M_new = vals.n + vals.n ÷ 3
    if vals.M <= M_new
        resize!(vals.op_type, M_new)
        resize!(vals.op_bond, M_new)
        resize!(vals.op_ham_idx, M_new)   
        resize!(vals.op_lbl, M_new)
        resize!(vals.vertex_list, M_new*2*inter.n_legs_max)         
        vals.op_type[vals.M+1:end] .= Int8(0)
        vals.op_bond[vals.M+1:end] .= 0
        vals.op_ham_idx[vals.M+1:end] .= Int8(0)
        vals.op_lbl[vals.M+1:end] .= 0
        vals.vertex_list[vals.M*2*inter.n_legs_max+1:end] .= 0
        vals.M = M_new
    end
end


# # TODO: implement for interactions with several types of bonds
# function off_diagonal_update!(; vals::Vals, inter::Interaction, mode::String, connections::Vector{<:Integer})
#     if mode ≠ "directed-loops-deterministic"
#         raise(error("Wrong mode of the off-diagonal update."))
#     end
#     for v0 in 1:inter.n_legs_max
#         if vals.vertex_list[v0] == 0
#             continue
#         end
#         v = v0
#         if rand() < 0.5
#             # The loop doesn't flip, but we have to exclude it from the vertex list.
#             while true
#                 vals.vertex_list[v] = 0
#                 v1 = connections[(v-1) % (2*inter.n_legs_max) + 1] + ((v-1)÷(2*inter.n_legs_max))*2*inter.n_legs_max
#                 v = vals.vertex_list[v1]
#                 vals.vertex_list[v1] = 0
#                 if v == v0
#                     break
#                 end
#             end
#         else
#             # The loop flips and excluded from the vertex list.
#             while true
# 
#             end
#         end
#     end
# end


