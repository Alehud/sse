export int2base, base2int

function int2base(n::Integer, b::Integer, num_digits::Integer; offset_zero::Bool=false)
    if (n < 0)
        raise(error("`n` has to be non-negative."))
    elseif b ≤ 0
        raise(error("`b` has to be positive."))
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


function base2int(digits::Vector{<:Integer}, b::Integer; offset_zero::Bool=false)
    digits_temp = copy(digits)
    if offset_zero
        digits_temp .-= 1
    end
    n = 0
    for (i, d) in enumerate(digits_temp[end:-1:1])
        n += d * b^(i-1)
    end
    return n
end
