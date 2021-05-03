using Test
using Plots
using LaTeXStrings

# calculate f(theta) for Jensen Wake Model
function f(theta)
    if abs(theta) > 20
        return 0
    else
        return .5*(1+cosd(9*theta))
    end
end

# calculate velocity of wind behind a turbine given ambient velocity u and the off axis angle theta
function vcosine(u,theta,r0,x,alpha)
    return u*(1-2/3*(r0/(r0+alpha*x))^2*f(theta))
end

# calculate velocity for the tophat model
function vtophat(u,r0,alpha,x)
    return u*(1-2/3*(r0/(r0+alpha*x))^2)
end

# generate plots of tophat vs cosine model
function generatePlots(ratio,legend,xlabel)
    theta = -30:0.1:30
    u = 5
    alpha = 0.1
    x = 40
    r0 = x/ratio
    v2 = zeros(size(theta))
    v1 = zeros(size(theta))
    vTop = vtophat(u,r0,alpha,x)
    thetaLim = 8
    for i = 1:length(theta)
        if abs(theta[i]) < thetaLim
            v1[i] = vTop
        else
            v1[i] = u
        end
        v2[i] = vcosine(u,theta[i],r0,x,alpha) 
    end

    # Check to add legend
    if legend
        p = plot(theta,v1/u,ylim=(0.7,1),title = latexstring("x/r_{0} = ",ratio),label="Top-hat Model")
        plot!(theta,v2/u,ylim=(0.7,1),label="Cosine Model",legend=:bottomright)
    else
        p = plot(theta,v1/u,ylim=(0.7,1),title = latexstring("x/r_{0} = ",ratio),label="Top-hat Model")
        plot!(theta,v2/u,ylim=(0.7,1),label="Cosine Model",legend=false)
    end

    # Check to add xlabel
    if xlabel
        xlabel!(L"$\theta\;(\degree)$")
    end
    ylabel!(L"$v/u$")
    return p
end

# testing functions
@test f(21) == 0
@test f(0) == 1
@test vcosine(6,20,8.1,40,0.1) == 6
@test vtophat(8.1,20,0.1,40) == 4.35
@test vtophat(8.1,20,0.1,100) == 5.7

p1 = generatePlots(16,true,false)
p2 = generatePlots(10,false,false)
p3 = generatePlots(6,false,true)

plot(p1,p2,p3,layout = (3,1))
plot!(size = (400,500))

savefig("jensen1983.pdf")
