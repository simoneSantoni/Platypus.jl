# Correctly rounded square root of a nonnegative exact rational. Float64 supplies
# only a seed; exact squared rounding-midpoint comparisons decide the result.
# This avoids changing MPFR's process-global precision during threaded fits.
function _rational_sqrt(q::Rational{BigInt})
    iszero(q) && return 0.0
    q > 0 || throw(DomainError(q, "square-root argument must be nonnegative"))
    e = fld(ndigits(numerator(q); base=2) - ndigits(denominator(q); base=2), 2)
    power = BigInt(1) << abs(2e)
    scaled = e >= 0 ? q / power : q * power
    value = ldexp(sqrt(_rational_float(scaled)), e)
    if isinf(value)
        value = floatmax(Float64)
    end
    while true
        exact = Rational{BigInt}(value)
        if value > 0
            lower = (exact + Rational{BigInt}(prevfloat(value))) / 2
            if q < lower^2 || (q == lower^2 && isodd(reinterpret(UInt64, value)))
                value = prevfloat(value)
                continue
            end
        end
        # The upper overflow midpoint is floatmax + half its final ulp.
        upper = value == floatmax(Float64) ? exact + (BigInt(1) << 970) :
            (exact + Rational{BigInt}(nextfloat(value))) / 2
        if q > upper^2 || (q == upper^2 && isodd(reinterpret(UInt64, value)))
            value = nextfloat(value)
            isinf(value) && return value
            continue
        end
        return value
    end
end

# Round an exact rational directly to binary64, independently of MPFR's ambient
# precision. Julia's Float64(::Rational{BigInt}) can otherwise round via BigFloat.
function _rational_float(q::Rational{BigInt})
    iszero(q) && return 0.0
    negative = q < 0
    n, d = abs(numerator(q)), denominator(q)
    e = ndigits(n; base=2) - ndigits(d; base=2)
    if e >= 0 ? n < (d << e) : (n << -e) < d
        e -= 1
    end
    e > 1023 && return negative ? -Inf : Inf
    quantum = max(e - 52, -1074)
    a, b = quantum < 0 ? (n << -quantum, d) : (n, d << quantum)
    m, remainder = divrem(a, b)
    twice = 2remainder
    if twice > b || (twice == b && isodd(m))
        m += 1
    end
    value = ldexp(Float64(UInt64(m)), quantum)
    return negative ? -value : value
end

