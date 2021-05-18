using Plots

# This file is to redo the code from katicBody.jl to allow for a better code organization

############################################################################
# Structs

# struct to hold wake model parameters
struct model
    alpha::Float64
end

# struct to hold turbine properties (diameter and power)
struct turbine
    diameter::Float64
end

# struct to hold wind properties
struct wind
    speeds::Vector{Float64}
    directions::Vector{Float64}
    probabilities::Vector{Float64}
end

############################################################################
# Functions 

# function to rotate the frame to deal with the wind
function rotateFarmIntoWind(turbine_x,turbine_y,angle)
    cosWind = cosd(-angle)
    sinWind = sind(-angle)
    xrot = turbine_x*cosWind - turbine_y*sinWind
    yrot = turbine_x*sinWind + turbine_y*cosWind
    return xrot,yrot
end

# function to calculate the wake deficit from wake model object and distance (Jensen)
function wakeDeficit(x,wake::model,turbine::turbine)
    r = turbine.diameter/2
    alpha = wake.alpha
    return 2/3*(r/(r+alpha*x))^2
end

# function to combine wake deficits (Katic)
function combineDeficits(deficit)
    total = 0
    for def in deficit
        total += def^2
    end
    return sqrt(total)
end

# function to iterate through the turbines for a given wind

# function to determine which wakes affect the current turbine

# function to determine what percent of turbine is affected by which wakes

# calculate wind farm total power

# set up farm beyond x and y coordinates (return wind, wake, and turbine parameters)
function buildFarm()
    alpha = 0.1 #jensen wake constant
    d = 10 #turbine diameter
    angle = [0] #wind direction measured in meteorolgoical coordinates
    speed = [10] #m/s
    prob = [1] #probabilities that the wind will be that direction and at that speed

    # construct object of turbine properties (diameter and power properties)
    turbine_properties = turbine(d)

    # construct object of wake properties and constants
    model_properties = model(alpha)

    # construct object of wind properties (speeds, directions, percent)
    wind_properties = wind(speed,angle,prob)

    return turbine_properties,model_properties,wind_properties
end

############################################################################
# CIRCULAR FORMATION
n = 10 #number of turbines
D = 200 #diameter of circle in meters

# construct x and y coordinates of turbines in two arrays
turbine_x = zeros(n)
turbine_y = zeros(n)
delta_theta = 360/n
for i = 1:n
    theta = delta_theta*(i-1)
    turbine_x[i] = D/2*cosd(theta)
    turbine_y[i] = D/2*sind(theta)
end