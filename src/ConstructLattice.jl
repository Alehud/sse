module ConstructLattice

using IterTools

export LatticePoint, construct_lattice, chain_NN, chain_NN_NNN, chain_NN_3_site, square_lattice_NN, square_lattice_plaquette, cubic_lattice_NN, honeycomb_lattice_NN,
       square_lattice_links_plaquette, square_lattice_links_star


"""
An arbitrary point of a lattice.
# Fields
- `unit_cell::Tuple{Vararg{<:Integer}}`: a tuple with coordinates of a unit cells
- `site::Integer`: a site within the unit cell

"""
struct LatticePoint
    unit_cell::Tuple{Vararg{<:Integer}}
    site::Integer
end


"""
    bare_ind(p, linear_size, sites_in_uc)

Calculate the bare index (i.e., an integer from 1 to N) of the lattice point `p`. 
"""
function bare_ind(p::LatticePoint, linear_size::Tuple{Vararg{<:Integer}}, sites_in_uc::Integer)
    if length(linear_size) ≠ length(p.unit_cell)
        raise(error("Dimensions don't match."))
    end
    ind = p.site
    for (i, x) in enumerate(p.unit_cell)
        ind += sites_in_uc * prod(linear_size[1:i-1]) * (x-1)
    end
    return ind
end


"""
    construct_lattice(linear_size, sites_in_uc, interaction_types)

Construct interactions (or "bonds") on a lattice. 

Each site of the lattice has a d.o.f. and is labeled with an integer (from 1 to N).
Create a 2D array `bonds`. i'th column of `bonds` contains indices of d.o.f.s participating in i'th interaction term (bond).
I.e., if `bonds[1,i] == a`, `bonds[2,i] == b`, ..., it means that i'th bond connects d.o.f.s a, b, ...
The number of rows in `bonds` is the maximal number of d.o.f.s in a single interaction. For interactions with smaller number of d.o.f.s, the remaining
rows are filled with 0.

# Arguments
- `linear_size::Tuple{Vararg{<:Integer}}`: tuple with the number of unit cells in each dimension
- `sites_in_uc::Integer`: the number of sites in a unit cell
- `interaction_types::Vector{Vector{LatticePoint}}`: vector of elementary interactions that will be translated along all dimensions of the lattice 
in order to get all interaction terms of the model. Each element of the vector is an interaction term, i.e., it is a vector of LatticePoints. Here, the 
actual operators of the interactions don't matter, only their geometry. E.g., if you have a 1D chain with two alternating types of interaction, you still 
only need one term in `interaction_types`.

"""
function construct_lattice(; linear_size::Tuple{Vararg{<:Integer}}, sites_in_uc::Integer, interaction_types::Vector{Vector{LatticePoint}})
    Nuc = prod(linear_size)     # number of unit cells
    # N = Nuc*sites_in_uc    # total number of sites
    Nb = Nuc*length(interaction_types)  # total number of interaction terms (or "bonds")
    # WARNING: provide the minimal possible set of interactions. There is no check whether any of the interactions in `interaction_types` are related
    # by translation.
    max_int_size = maximum(map(length, interaction_types))    # maximal number of degrees of freedom in a single interaction
    bonds = zeros(Int32, max_int_size, Nb)
    
    b = 1
    for interaction in interaction_types
        # Iterate over 1:linear_size[1], 1:linear_size[2], 1:linear_size[3], etc.
        for i in Iterators.product(Base.OneTo.(linear_size)...)
            # Iterate over lattice points in each interaction type
            for (p_ind, p) in enumerate(interaction)
                # Shift the unit cell by the vector `i` (with p.b.c.), calculate the bare index of the site, and fill in `bonds`.
                pp = LatticePoint((p.unit_cell .+ i .- 1 .- 1) .% linear_size .+ 1 , p.site)
                bonds[p_ind, b] = bare_ind(pp, linear_size, sites_in_uc)
            end
            b += 1
        end
    end
    
    return bonds
end


"""
    chain_NN(L)

Construct a 1D chain of length `L` with nearest neighbor two-site interactions (and p.b.c).
"""
function chain_NN(; L::Int)
    interaction_types = [[LatticePoint(Tuple(1), 1), LatticePoint(Tuple(2), 1)]]
    return construct_lattice(linear_size=Tuple(L), sites_in_uc=1, interaction_types=interaction_types)
end


"""
    chain_NN_NNN(L)

Construct a 1D chain of length `L` with nearest neighbor two-site interactions and next-nearest neighbor two-site interactions (and p.b.c).
"""
function chain_NN_NNN(; L::Int)
    interaction_types = [[LatticePoint(Tuple(1), 1), LatticePoint(Tuple(2), 1)], [LatticePoint(Tuple(1), 1), LatticePoint(Tuple(3), 1)]]
    return construct_lattice(linear_size=Tuple(L), sites_in_uc=1, interaction_types=interaction_types)
end


"""
    chain_NN_3_site(L)

Construct a 1D chain of length `L` with nearest neighbor three-site interactions (and p.b.c).
"""
function chain_NN_3_site(; L::Int)
    interaction_types = [[LatticePoint(Tuple(1), 1), LatticePoint(Tuple(2), 1), LatticePoint(Tuple(3), 1)]]
    return construct_lattice(linear_size=Tuple(L), sites_in_uc=1, interaction_types=interaction_types)
end


"""
    square_lattice_NN(Lx, Ly)

Construct a 2D square lattice with nearest neighbor two-site interactions (and p.b.c).
"""
function square_lattice_NN(; Lx::Int, Ly::Int)
    interaction_types = [[LatticePoint((1,1), 1), LatticePoint((2,1), 1)], [LatticePoint((1,1), 1), LatticePoint((1,2), 1)]]
    return construct_lattice(linear_size=(Lx, Ly), sites_in_uc=1, interaction_types=interaction_types)
end


"""
    square_lattice_plaquette(Lx, Ly)

Construct a 2D square lattice with plaquette four-site interactions (and p.b.c).
"""
function square_lattice_plaquette(; Lx::Int, Ly::Int)
    interaction_types = [[LatticePoint((1,1), 1), LatticePoint((2,1), 1), LatticePoint((1,2), 1), LatticePoint((2,2), 1)]]
    return construct_lattice(linear_size=(Lx, Ly), sites_in_uc=1, interaction_types=interaction_types)
end


"""
    cubic_lattice_NN(Lx, Ly, Lz)

Construct a 3D cubic lattice with nearest neighbor two-site interactions (and p.b.c).
"""
function cubic_lattice_NN(; Lx::Int, Ly::Int, Lz::Int)
    interaction_types = [[LatticePoint((1,1,1), 1), LatticePoint((2,1,1), 1)], [LatticePoint((1,1,1), 1), LatticePoint((1,2,1), 1)], 
    [LatticePoint((1,1,1), 1), LatticePoint((1,1,2), 1)]]
    return construct_lattice(linear_size=(Lx, Ly, Lz), sites_in_uc=1, interaction_types=interaction_types)
end


"""
    honeycomb_lattice_NN(L_zigzag, L_armchair)

Construct a 2D honeycomb lattice with nearest neighbor two-site interactions (and p.b.c). Along one direction, the lattice is glued in a zigzag way.
Along the other direction, it is glued in an armchair way. It is glued without a twist, meaning that if we choose a unit cell and translate it along the 
lattice vector L_zigzag (L_armchair) times, we will come to the initial cell, rather than wrap around the torus. WARNING: note that when we draw a honeycomb 
lattice on paper, an enticing way to glue the edges (connecting the opposite vertical/horizontal sides) results in a twist in the boundary.
"""
function honeycomb_lattice_NN(; L_zigzag::Int, L_armchair::Int)
    interaction_types = [[LatticePoint((1,1), 1), LatticePoint((1,1), 2)], [LatticePoint((1,1), 2), LatticePoint((2,1), 1)], 
    [LatticePoint((1,1), 2), LatticePoint((2,2), 1)]]
    return construct_lattice(linear_size=(L_zigzag, L_armchair), sites_in_uc=2, interaction_types=interaction_types)
end


"""
    square_lattice_links_plaquette(Lx, Ly)

Construct a 2D square lattice with d.o.f. on the links (i.e., two sites in a unit cell), with plaquette four-site interactions (and p.b.c).
"""
function square_lattice_links_plaquette(; Lx::Int, Ly::Int)
    interaction_types = [[LatticePoint((1,1), 1), LatticePoint((1,1), 2), LatticePoint((2,1), 2), LatticePoint((1,2), 1)]]
    return construct_lattice(linear_size=(Lx, Ly), sites_in_uc=2, interaction_types=interaction_types)
end


"""
    square_lattice_links_star(Lx, Ly)

Construct a 2D square lattice with d.o.f. on the links (i.e., two sites in a unit cell), with star four-site interactions (and p.b.c).
"""
function square_lattice_links_star(; Lx::Int, Ly::Int)
    interaction_types = [[LatticePoint((1,1), 1), LatticePoint((1,1), 2), LatticePoint((Lx,1), 1), LatticePoint((1,Ly), 2)]]
    return construct_lattice(linear_size=(Lx, Ly), sites_in_uc=2, interaction_types=interaction_types)
end


"""
    square_lattice_links_plaquette_star(Lx, Ly)

Construct a 2D square lattice with d.o.f. on the links (i.e., two sites in a unit cell), with plaquette and star four-site interactions (and p.b.c).
"""
function square_lattice_links_plaquette_star(; Lx::Int, Ly::Int)
    interaction_types = [[LatticePoint((1,1), 1), LatticePoint((1,1), 2), LatticePoint((2,1), 2), LatticePoint((1,2), 1)]]
    push!(interaction_types, [LatticePoint((1,1), 1), LatticePoint((1,1), 2), LatticePoint((Lx,1), 1), LatticePoint((1,Ly), 2)])
    return construct_lattice(linear_size=(Lx, Ly), sites_in_uc=2, interaction_types=interaction_types)
end







end