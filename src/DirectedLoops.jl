module DirectedLoops

using ..Interactions
using LinearAlgebra
using IterTools
using JuMP
using GLPK
using Graphs

export solve_set, is_compatible_with_dir_loops


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
    for ham in inter.hams
        # We will create a graph `g` with vertices representing allowed vertex configurations and edges connecting vertices that can be brought to each other by
        # changing no more than 2 legs. We also include all diagonal configurations in the graph (even if they are not allowed). Ids of disallowed diagonal
        # vertices are stored in `excluded_vertices`. In `g_items` we store vertices of the graph in the format (ii, jj), where ii and jj are tuples denoting
        # the bottom and the top configurations of the vertex.
        g_items = Tuple[]
        excluded_vertices = Int[]
        
        # Construct `g_items` and `excluded_vertices`.
        ss = size(ham)
        bond_size = round(Int, log(inter.dof_max, ss[1]))
        range = [0:inter.dof_max-1;]
        powers = inter.dof_max .^ [0:bond_size-1;]
        iterator = fill(range, bond_size)
        for ii in IterTools.product(iterator...), jj in IterTools.product(iterator...)
            i = sum(ii .* powers) + 1
            j = sum(jj .* powers) + 1
            if ham[i, j] ≠ 0
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



end