export Observables, ObservablesC, accumulate_bin_values!, measure!, measureC!

@with_kw mutable struct Observables
    n_accum::Integer
    n²_accum::Integer
    mag_stag²_accum::Real
    mag²_accum::Real
    n_av::Real
    n²_av::Real
    E::Real
    C::Real
    mag_stag²::Real
    mag²::Real
    χ::Real
    function Observables(n_accum::Integer, n²_accum::Integer, mag_stag²_accum::Real, mag²_accum::Real, n_av::Real, n²_av::Real, E::Real, C::Real, mag_stag²::Real, mag²::Real, χ::Real)
        # TODO: add checks
        new(n_accum, n²_accum, mag_stag²_accum, mag²_accum, n_av, n²_av, E, C, mag_stag², mag², χ)
    end
end

function Observables()
    return Observables(0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

@with_kw mutable struct ObservablesC
    enrg1::Real
    enrg2::Real
    amag2::Real
    ususc::Real
    function ObservablesC(enrg1::Real, enrg2::Real, amag2::Real, ususc::Real)
        # TODO: add checks
        new(enrg1, enrg2, amag2, ususc)
    end
end

function ObservablesC()
    return ObservablesC(0.0, 0.0, 0.0, 0.0)
end


function accumulate_bin_values!(; vals::Vals, inter::Interaction, bipartition_mask::Vector{<:Integer}, obs::Observables)
    obs.n_accum += vals.n
    obs.n²_accum += (vals.n)^2
    obs.mag²_accum += sum(vals.dofs)^2

    # TODO: implement more general static structure factor S(q) instead of the square of the staggered magnetization, which is S(π). (page 156 of Sandvik 2010)
    if vals.n ≠ 0
        mag_stag²_accum_bin = 0
        for m in 1:vals.M
            # Only collect values from layers with an operator (diagonal or off-diagonal).
            if vals.op_type[m] ≠ 0
                mag_stag = sum(bipartition_mask .* (vals.dofs .- 1.5))     # vals.dofs are = 1 or 2, but we want them to be = -0.5 or +0.5
                mag_stag²_accum_bin += mag_stag^2
            end
            # If there is an off-diagonal operator, advance the state
            if vals.op_type[m] == 2
                vals.dofs[inter.bond_map[:,vals.op_bond[m]]] = inter.hams[vals.op_ham_idx[m]].labels_dof[vals.op_lbl[m]][2]
            end
        end
        obs.mag_stag²_accum += mag_stag²_accum_bin / (vals.n * inter.n_dofs)
    else
        obs.mag_stag²_accum += sum(bipartition_mask .* (vals.dofs .- 1.5)) / inter.n_dofs
    end
end

function measure!(; vals::Vals, inter::Interaction, steps_bin::Integer, beta::Real, obs::Observables)
    obs.n_av = obs.n_accum/steps_bin
    obs.n²_av = obs.n²_accum/steps_bin
    obs.E = -obs.n_av/beta                              # internal energy
    obs.C = obs.n²_av - obs.n_av^2 - obs.n_av                   # specific heat
    obs.mag_stag² = obs.mag_stag²_accum / steps_bin     # square of the staggered magnetization
    obs.mag² = obs.mag²_accum / steps_bin               # square of the magnetization
    obs.χ = beta * obs.mag² / inter.n_dofs              # uniform magnetic susceptibility 
end

function measureC!(; vals::Vals, inter::Interaction, steps_bin::Integer, beta::Real, obs::ObservablesC)
    Lx = 4
    am = 0
    for s in 0:(inter.n_dofs - 1)
        am += (-1)^vals.dofs[s+1] * (-1)^(s%Lx + s÷Lx)
    end
    am = am÷2
    am2 = 0.0
    for m in 1:vals.M
        op = vals.op_type[m]
        bond = vals.op_bond[m]
        if op == 2
            s1 = inter.bond_map[1, bond]
            s2 = inter.bond_map[2, bond]
            vals.dofs[s1] += (-1)^(vals.dofs[s1] - 1)
            vals.dofs[s2] += (-1)^(vals.dofs[s2] - 1)
            am += 2 * (-1)^vals.dofs[s1] * (-1)^((s1-1)%Lx + (s1-1)÷Lx)
        end
        am2 += am^2
    end
    am2 = am2 / vals.M

    obs.enrg1 += vals.n
    obs.enrg2 += vals.n^2
    obs.amag2 += am2
    spin_sum = sum(vals.dofs .- 1.5)
    obs.ususc += spin_sum^2
end