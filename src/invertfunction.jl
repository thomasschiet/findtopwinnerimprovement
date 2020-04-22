function getC(K, δ)
    if δ == 1
        return K/2
    end

    @assert 0 < δ < 1

    x_left = BigFloat(K/2)
    x_right = BigFloat(K)

    f(x) = exp(-x)*(2*exp(1)*x/K)^(K/2) - δ

    y_left = f(x_left)
    y_mid = f((x_left + x_right) / 2)
    y_right = f(x_right)

    @assert y_left >= 0

    while y_right > 0
        x_right *= 2
        y_right = f(x_right)
    end

    y_mid = f((x_left + x_right) / 2)

    while abs(y_mid) > eps(Float64)
        if y_mid < 0
            x_left = x_left
            x_right = (x_left + x_right)/2
        elseif y_mid > 0
            x_left = (x_left + x_right)/2
            x_right = x_right
        else
            return (x_left + x_right)/2
        end
        y_mid = f((x_left + x_right) / 2)
    end

    Float64((x_left + x_right)/2)
end
