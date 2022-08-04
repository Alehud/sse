export is_compatible_with_dir_loops, solve_set, solve_directed_loop_equations


"""
    is_compatible_with_dir_loops(inter)

Check if the directed loop update scheme can be used for the given interactions. 
    
An interaction is compatible with the directed loop update scheme if every allowed vertex configuration can be connected to a diagonal vertex configuration
through a number of 2-site updates on this vertex. The interaction might become compatible with the directed loop update scheme upon adding a constant 
to the Hamiltonian. In this case, a respective notification will be printed.

# Arguments
- `inter::Interaction`: interaction

"""
function is_compatible_with_dir_loops(inter::Interaction)
    compatible = true
    compatible_with_const = true

    # Iterate over Hamiltonians in the interaction.
    for H in inter.hams
        # We will create a graph `g` with vertices representing allowed vertex configurations and edges connecting vertices that can be brought to each other by
        # changing no more than 2 legs. We also include all diagonal configurations in the graph (even if they are not allowed). Ids of disallowed diagonal
        # vertices are stored in `excluded_vertices`. In `g_items` we store vertices of the graph in the format (ii, jj), where ii and jj are tuples denoting
        # the bottom and the top configurations of the vertex.
        g_items = Tuple[]
        excluded_vertices = Int[]
        
        # Construct `g_items` and `excluded_vertices`.
        range = [0:inter.dof_max-1;]
        powers = inter.dof_max .^ [0:H.n_legs-1;]
        iterator = fill(range, H.n_legs)
        for ii in IterTools.product(iterator...), jj in IterTools.product(iterator...)
            i = sum(ii .* powers) + 1
            j = sum(jj .* powers) + 1
            if H.ham[i, j] ≠ 0
                tup = (ii, jj)
                push!(g_items, tup)
            elseif i == j
                tup = (ii, jj)
                push!(g_items, tup)
                push!(excluded_vertices, length(g_items))
            end
        end

        # Construct the graph `g`.
        g = SimpleGraph(length(g_items))
        for (k, (ii, jj)) in enumerate(g_items)
            for (l, (mm, nn)) in enumerate(g_items[k+1:end])
                if 1 ≤ count((mm .- ii) .≠ 0) + count((nn .- jj) .≠ 0) ≤ 2
                    add_edge!(g, k, k+l)
                end
            end
        end

        # For each off-diagonal vertex, we find if there is a path in the graph that connects it to any of the diagonal vertices.
        # We do it both excluding and including disallowed diagonal vertices.

        # Iterate through off-diagonal elements of the graph
        for (k, (ii, jj)) in enumerate(g_items)
            if ii ≠ jj
                is_path = false
                is_path_with_const = false
                # Iterate through diagonal elements of the graph
                for (l, (mm, nn)) in enumerate(g_items)
                    if mm == nn
                        if has_path(g, k, l, exclude_vertices=excluded_vertices)
                            is_path = true
                            is_path_with_const = true
                            break
                        elseif has_path(g, k, l)
                            is_path_with_const = true
                        end
                    end
                end
                compatible &= is_path
                compatible_with_const &= is_path_with_const
            end
        end
    end

    if !compatible & compatible_with_const
        println("The interactions are not compatible with the directed loop update scheme, but they will be upon adding a constant to the Hamiltonian.")
    end

    return compatible
end


"""
    solve_set(weights, mode)

Solve one set of directed loop equations. For explanation, see Siljuåsen and Sandvik, 2002 (section III). 

For example, directed loop equations can look like
W1 = a11 + a12 + a13
W2 = a21 + a22 + a23
W3 = a31 + a32 + a33
where W1, W2, W3 are the vertex weights (i.e. matrix elements of a local Hamiltonian) determined by the legs' configuration s (s=1,2,3). 
a_ij are the probabilities that given vertex state i and some chosen entrance leg, the directed loop will transition into the vertex state j going through some exit leg. 
Because the detailed balance has to be satisfied, a_ij = a_ji.
We want to solve for probabilities a_ij. Since the set of equations is underdetermined, there are many solutions.

# Arguments
- `weights::Vector{<:Real}`: vertex weights W_i
- `mode`: how to solve the directed loop eqations. "heat-bath" mode uses the heat bath solution. "linear-programming" uses a linear programming solver to minimize the probabilities of bounces (a_ii).
"""
function solve_set(; weights::Vector{<:Real}, mode="heat-bath")
    n = length(weights)
    w_sum = sum(weights)
    A = Array{Float64}(undef, n, n)

    # Sum of variables in row i has to be equal to weights[i], since this is the directed loops equations. Eq. (28) from Siljuåsen, Sandvik 2002.
    if mode == "heat-bath"
        for i in 1:n, j in 1:n
            A[i, j] = weights[i] * weights[j] / w_sum
        end
    elseif mode == "linear-programming"
        # Initialize linear programming solver
        model = Model()
        set_optimizer(model, GLPK.Optimizer)
        
        # Define variables
        @variable(model, a[i = 1:n, j = i:n] >= 0)

        # Define Constraints
        @constraint(model, c[i = 1:n], sum([i ≤ j ? a[i,j] : a[j,i] for j in 1:n]) == weights[i])

        # Define Objective
        @objective(model, Min, sum([a[i,i] for i in 1:n]))

        # Run the opimization
        optimize!(model)

        # Check if the solver found a solution
        status = termination_status(model)
        if status == OPTIMAL
            # println("Linear programming solver found (and proved) a globally optimal solution of the directed loop equations.")
        elseif status == LOCALLY_SOLVED
            println("Linear programming solver found a locally optimal solution of the directed loop equations (which may also be globally optimal, but it could not prove so).")
        else
            raise(error("Linear programming solver could not find a solution. Termination status: ", status))
        end

        for i in 1:n, j in 1:n
            if i ≤ j
                A[i,j] = value(a[i,j])
            else
                A[i,j] = value(a[j,i])
            end
        end
    else
        raise(error("Invalid mode!"))
    end

    # Normalize each row to get transition probabilities
    for i in 1:n
        A[i,:] ./= weights[i]
    end
    
    return A
end


"""
    solve_directed_loop_equations(inter)

Solve directed loop equations and construct a tensor with probabilities that can be directly used to determine where the loop goes next once it enters a certain vertex.
The final tensor looks like P[type of hamiltonian][label, entrance leg, new state of entrance leg, exit leg, new state of exit leg].

# Arguments
- `inter::Interaction`: interaction

"""
function solve_directed_loop_equations(inter::Interaction)
    P = Vector{Array}(undef, length(inter.hams))

    for (ham_idx, H) in enumerate(inter.hams)
        prob = fill(Float32(-1.0), (length(H.labels), 2H.n_legs, H.dof_max, 2H.n_legs, H.dof_max))
        # prob[label, entrance leg, new state of entrance leg, exit leg, new state of exit leg]
        # If entrance and exit legs are the same, then the "new state of exit leg" decides the new state of this leg (not the "new state of entrance leg")!

        for (idx, (label, label_dof)) in enumerate(zip(H.labels, H.labels_dof)), legₑ in 1:2H.n_legs, sₑ in 1:H.dof_max
            
            old_s_at_legₑ = label_dof[(legₑ-1) ÷ H.n_legs + 1][(legₑ-1) % H.n_legs + 1]

            # If already assigned, skip.
            if prob[idx, legₑ, sₑ, legₑ, old_s_at_legₑ] ≥ 0
                for legₓ in 1:2H.n_legs, sₓ in 1:H.dof_max
                    if prob[idx, legₑ, sₑ, legₓ, sₓ] < 0
                        prob[idx, legₑ, sₑ, legₓ, sₓ] = 0
                    end
                end
                continue
            end

            # Construct a set of directed loop equations by iterating over all exit legs and all new states of the exit leg.
            labels_set = CartesianIndex[]  # vector of Cartesian indices corresponding to labels in rows of the set of directed loop equations
            idx_set = Integer[]
            legₓ_set = Integer[] # vector of exit legs in columns of the set
            sₓ_set = Integer[]  # vector of new exit leg states in columns of the set
            legₑ_set = Integer[] # vector of entrance legs in rows of the set
            sₑ_set = Integer[]  # vector of new entrance leg states in rows of the set
            push!(labels_set, label)
            push!(idx_set, idx)
            push!(legₓ_set, legₑ)
            push!(sₓ_set, old_s_at_legₑ)
            push!(legₑ_set, legₑ)
            push!(sₑ_set, sₑ)

            for legₓ in 1:2H.n_legs, sₓ in 1:H.dof_max
                # Construct the label of the new configuration.
                new_label_dof = (copy(label_dof[1]), copy(label_dof[2]))
                new_label_dof[(legₑ-1) ÷ H.n_legs + 1][legₑ - H.n_legs*((legₑ-1) ÷ H.n_legs)] = sₑ
                new_label_dof[(legₓ-1) ÷ H.n_legs + 1][legₓ - H.n_legs*((legₓ-1) ÷ H.n_legs)] = sₓ
                new_label = CartesianIndex(base2int(new_label_dof[1], H.dof_max, offset_zero=true) + 1, base2int(new_label_dof[2], H.dof_max, offset_zero=true) + 1)
                
                # Check if the matrix element of the new configuration is zero. If yes, skip.
                if H.ham[new_label] == 0
                    prob[idx, legₑ, sₑ, legₓ, sₓ] = 0
                    continue
                elseif (new_label == label) && (legₑ == legₓ) && (sₓ == old_s_at_legₑ)
                else
                    new_idx = findall(x -> x == new_label, H.labels)[1]  # linear index of the new label (i.e. H.labels[new_idx] == new_label)
                    new_legₑ = legₓ
                    new_sₑ = label_dof[(legₓ-1) ÷ H.n_legs + 1][(legₓ-1) % H.n_legs + 1]
                    push!(labels_set, new_label)
                    push!(idx_set, new_idx)
                    push!(legₑ_set, new_legₑ)
                    push!(sₑ_set, new_sₑ)
                    push!(legₓ_set, legₓ)
                    push!(sₓ_set, sₓ)
                end
            
            end
            weights = [H.ham[c] for c in labels_set]
            A = solve_set(weights=weights, mode="linear-programming")

            # Assign probabilities to the corresponding tensor elements of `prob`
            for (i, (lbl_idx, leg_ent, s_ent)) in enumerate(zip(idx_set, legₑ_set, sₑ_set)), (j, (leg_exit, s_exit)) in enumerate(zip(legₓ_set, sₓ_set))
                prob[lbl_idx, leg_ent, s_ent, leg_exit, s_exit] = A[i,j]
            end
        end

        if any(sum(prob, dims=(4,5)) .≠ 1.0)
            raise(error("Probabilities don't sum up to 1."))
        end
        P[ham_idx] = prob
    end
    println("Directed loop equations solved.")
    return P
end
