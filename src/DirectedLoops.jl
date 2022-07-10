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

        # println(model)
        # println(solution_summary(model, verbose=true))

        # Check if the solver found a solution
        status = termination_status(model)
        if status == OPTIMAL
            println("Linear programming solver found (and proved) a globally optimal solution of the directed loop equations.")
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


function solve_directed_loop_equations(inter::Interaction)
    probabilities = Vector{Array}(undef, length(inter.hams))
    println("inter.bond_map:"); display(inter.bond_map); println("\n")

    for (i, (ham, ham_bonds)) in enumerate(zip(inter.hams, inter.int_bonds))
        println(i)
        display(ham); println("\n")
        println(ham_bonds)
        num_legs = 2*length(inter.bonds[:, ham_bonds[1]])   # number of operator's legs
        num_states = length(ham)  # number of possible vertex states (i.e., the number of matrix elements in the Hamiltonian)
        P = Array{Any}(undef, num_legs, num_states) # first dimension is the number of legs, the second dimension is the number of possible vertex states
        
        println(size(P))
        for s_down in 1:size(ham)[1], s_up in 1:size(ham)[2]
            s = (s_down, s_up)  # states of bottom legs (s_down) and top legs (s_up) of the vertex
            s_ind = (s_down - 1)*size(ham)[1] + (s_up - 1) + 1  # index of the state (from 1 to num_states) -> to write down the probability in P
            for leg in 1:num_legs
                if !isassigned(P, leg, s_ind)
                    println("")
                end
            end
        end
    end
end
