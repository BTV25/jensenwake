using Test

function f(theta)
    if theta > 20
        return 0
    else
        return .5*(1+cosd(9*theta))
    end
end

function v(u,theta)
    return u*(1-2/3*(f(theta))^2)
end

@test f(21) == 0
@test f(0) == 1
@test v(6,0) == 2
@test v(6,20) == 6

