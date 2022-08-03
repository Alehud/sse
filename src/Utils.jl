export int2base, base2int

"""
    int2base(n, b, num_digits, offset_zero=false)

Convert an integer `n` to a vector containing the representation of `n` in a given base `b`. The total number of vector elements is `num_digits`.

# Arguments
- `n::Integer`: a (non-negative) integer to be converted
- `b::Integer`: the base (positive integer > 1)
- `num_digits::Integer`: total number of digits in the representation (length of the resulting vector). Has to fit all non-zero digits of `n` in base `b`.
- `offset_zero::Bool=false`: if true, all the digits in the final representation are incremented by 1.

"""
function int2base(n::Integer, b::Integer, num_digits::Integer; offset_zero::Bool=false)
    if (n < 0)
        raise(error("`n` has to be non-negative."))
    elseif b ≤ 1
        raise(error("`b` has to be positive > 1."))
    elseif num_digits ≤ 0
        raise(error("`num_digits` has to be positive."))
    end

    digits = fill(0, num_digits)
    for i in 1:num_digits
        digits[i] = n % b
        n = n ÷ b
    end
    if n ≠ 0
        raise(error("The number ", n ," does not fit in ", num_digits, " digits in base ", b, "."))
    end
    if offset_zero
        digits .+= 1
    end
    return digits[end:-1:1]
end


"""
    base2int(digits, b, offset_zero=false)

Convert a vector `digits` containing the digits in a given base `b` to the corresponding integer.

# Arguments
- `digits::Vector{<:Integer}`: a vector of non-negative integers (each less than `b`)
- `b::Integer`: the base (positive integer)
- `offset_zero::Bool=false`: if true, all the digits in the final representation are incremented by 1.

"""
function base2int(digits::Vector{<:Integer}, b::Integer; offset_zero::Bool=false)
    digits_temp = copy(digits)
    if offset_zero
        digits_temp .-= 1
    end

    if b ≤ 1
        raise(error("`b` has to be positive > 1."))
    elseif any(digits_temp .< 0)
        raise(error("`digits` cannot contain negative integers."))
    elseif any(digits_temp .≥ b)
        raise(error("Integers in `digits` have to be less than base `b`."))
    end
    
    n = 0
    for (i, d) in enumerate(digits_temp[end:-1:1])
        n += d * b^(i-1)
    end
    return n
end
