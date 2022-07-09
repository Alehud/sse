export int2base

function int2base(n::Integer, b::Integer, num_digits::Integer)
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
    return digits[end:-1:1]
end
