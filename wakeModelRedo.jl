using Plots
using Test

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
    cutIn::Float64
    cutOut::Float64
    ratedSpeed::Float64
    ratedPower::Float64
    efficiency::Float64
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
# take turbine_x,turbine_y,angle and return xrot,yrot
function rotateFarmIntoWind(turbine_x,turbine_y,angle)
    cosWind = cosd(-angle)
    sinWind = sind(-angle)
    xrot = turbine_x*cosWind - turbine_y*sinWind
    yrot = turbine_x*sinWind + turbine_y*cosWind
    return xrot,yrot
end

# function to calculate the wake deficit from wake model object and distance (Jensen)
function wakeDeficit(x,wake::model,turbines::turbine)
    r = turbines.diameter/2
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

# function to determine the parameters of the wakes effect
function wakeHits(leftWakeEdge,rightWakeEdge,leftTurbineEdge,rightTurbineEdge,diameter)
    right = false
    left = false

    if rightWakeEdge <= rightTurbineEdge # right edge of wake hits turbine
        right = true
    end
    if leftWakeEdge >= leftTurbineEdge # left edge of wake hits turbine
        left = true
    end

    # define wake effect area with center and radius relative to the diameter
    if right && left
        # find where both edges of the wake hit the turbine relative to turbine's left edge
        # in terms of the diameter
        leftHit = (leftWakeEdge - leftTurbineEdge)/diameter
        rightHit = (rightWakeEdge - leftTurbineEdge)/diameter
    elseif right
        # find where right edge of wake hits the turbine relative to left edge
        leftHit = 0
        rightHit = (rightWakeEdge - leftTurbineEdge)/diameter
    elseif left
        # find where left edge of wake hits the turbine relative to the right edge
        leftHit = (leftWakeEdge - leftTurbineEdge)/diameter
        rightHit = 1
    elseif !right && !left
        leftHit = 0
        rightHit = 1
    end

    return leftHit,rightHit
end

# d is in terms of the diameter (0.5 D)
# calculates the area of a slice that goes to the edges of the turbine
function calculateAreaOfTurbineSegment(d)
    theta = 2*acos(1-2*d)
    return (0.5*(theta-sin(theta)))/pi
end

# function to calculate slice size from d1 and d2
function calculatePercentAreaOfTurbineSlice(d1,d2)
    if d1 == d2
        return 0
    end

    if d2 < d1
        temp = d1
        d1 = d2
        d2 = temp
    end

    if d1 < 0.5
        if d2 <= 0.5
            return calculateAreaOfTurbineSegment(d2) - calculateAreaOfTurbineSegment(d1)
        elseif d2 > 0.5
            return 1 - calculateAreaOfTurbineSegment(d1) - calculateAreaOfTurbineSegment(1-d2)
        end

    elseif d1 >= 0.5
        if d2 > 0.5
            return calculateAreaOfTurbineSegment(1-d1) - calculateAreaOfTurbineSegment(1-d2)
        end
    end
    error("Error in calculatePercentAreaOfTurbineSlice")
end

# function to calculate power from single turbine
function deficit_for_single_turbine(xrot,yrot,turbine_properties,model_properties,wind_properties,currentTurbine)
    leftTurbineEdge = xrot[currentTurbine] - .5*turbine_properties.diameter
    rightTurbineEdge = xrot[currentTurbine] + .5*turbine_properties.diameter

    n = length(xrot)
    # arrays to hold the deficit of the wake and where it hits the turbine
    leftHits = zeros(n)
    rightHits = zeros(n)
    deficits = zeros(n)

    # iterate through each other turbine and calculate the effect of each turbine
    for i = 1:n
        # if upwind turbine is current turbine then ignore
        if i == currentTurbine
            continue

        # if upwind turbine is behind then ignore
        elseif yrot[i] >= yrot[currentTurbine]
            continue
        end
        
        # calculate wake range at the currentTurbine
        range = 0.5*turbine_properties.diameter+model_properties.alpha*(abs(yrot[i]-yrot[currentTurbine]))
        rightWakeEdge = xrot[i] + range
        leftWakeEdge = xrot[i] - range

        # if the wake does not overlap the turbine ignore it
        if (leftTurbineEdge >= rightWakeEdge) || (rightTurbineEdge <= leftWakeEdge)
            continue
        end

        deficits[i] = wakeDeficit(abs(yrot[i]-yrot[currentTurbine]),model_properties,turbine_properties)
        
        leftHits[i],rightHits[i] = wakeHits(leftWakeEdge,rightWakeEdge,leftTurbineEdge,rightTurbineEdge,turbine_properties.diameter)
    end
    # remove wakes with no effect
    remove = []
    for i = 1:n
        if deficits[i] == 0
            remove = vcat(remove,i)
        end
    end
    deleteat!(deficits,remove)
    deleteat!(leftHits,remove)
    deleteat!(rightHits,remove)
    if isempty(leftHits)
        # no deficit exists
        return 0
    end

    # determine left to right order of wake hits
    wakes = [leftHits rightHits deficits]
    order = zeros(2*length(leftHits),3)

    for i = 1:length(leftHits)
        order[i,1] = i
        order[i,2] = 1
        order[i,3] = wakes[i,1]
        order[i+length(leftHits),1] = i
        order[i+length(leftHits),2] = 2
        order[i+length(leftHits),3] = wakes[i,2]
    end
    order = sortslices(order,dims=1,by=x->x[3],rev=false)

    # order is a list ordered by the third column with columns 1 and 2 being
    # the row and column in wake at which the edge is located
    active = convert(Array{Int8},zeros(length(leftHits),1))

    active[convert(Int32,order[1,1])] = 1
    d1 = order[1,3]

    deficit = 0
    for i = 2:size(order,1)
        col = convert(Int32,order[i,2])
        d2 = order[i,3]

        # use d1 and d2 to get slice size
        slice = calculatePercentAreaOfTurbineSlice(d1,d2)

        # use slice size and multiply by deficit for the area
        def = combineDeficits(wakes[:,3] .* active)
        deficit += def * slice

        if col == 1
            active[convert(Int32,order[i,1])] = 1
        else
            active[convert(Int32,order[i,1])] = 0
        end
        d1 = d2
    end
    return deficit
end

# calculate power for a state
function power_for_state(turbine_x,turbine_y,turbine_properties,model_properties,wind_properties,num)
    # rotate the frame
    xrot,yrot = rotateFarmIntoWind(turbine_x,turbine_y,wind_properties.directions[num])

    deficits = zeros(length(xrot))
    # calculate deficit for each turbine
    for i = 1:length(xrot)
        deficits[i] = deficit_for_single_turbine(xrot,yrot,turbine_properties,model_properties,wind_properties,i)
    end
    windSpeeds = 1 .- deficits

    energy = sum(windSpeeds.^3)/length(xrot)

    windSpeeds = windSpeeds .* wind_properties.speeds[num]

    # calculate power from deficits
    power = 0
    for j = 1:length(deficits)
        power += calculatePower(windSpeeds[j],wind_properties,turbine_properties)
    end

    return power,energy
end

function calculatePower(wind_speed,wind_properties,turbine_properties)
    cutIn = turbine_properties.cutIn
    cutOut = turbine_properties.cutOut
    ratedSpeed = turbine_properties.ratedSpeed
    ratedPower = turbine_properties.ratedPower

    if wind_speed < cutIn
        power = 0.0
    elseif wind_speed < ratedSpeed
        power = ratedPower*((wind_speed-cutIn)/(ratedSpeed-cutIn))^3
    elseif wind_speed < cutOut
        power = ratedPower
    elseif wind_speed > cutOut
        power = 0.0
    end

    return power
end

# calculate wind farm total power
function totalPower(turbine_x,turbine_y)
    turbine_properties,model_properties,wind_properties = buildFarm()
    power = 0
    energy = zeros(length(wind_properties.speeds))
    for i = 1:length(wind_properties.speeds)
        singlePower,energy[i] = power_for_state(turbine_x,turbine_y,turbine_properties,model_properties,wind_properties,i)
        power += wind_properties.probabilities[i]*singlePower
    end

    if length(wind_properties.speeds) == 1
        plots = buildWakePlots(turbine_properties,model_properties,wind_properties,turbine_x,turbine_y)
    else
        plots = []
    end

    return power,energy,plots
end

# set up farm beyond x and y coordinates (return wind, wake, and turbine parameters)
function buildFarm()
    alpha = 0.1 #jensen wake constant
    d = 10 #turbine diameter
    # global allAngles = collect(270:306) #wind direction measured in meteorolgoical coordinates
    # angle = allAngles
    angle = [288]
    speed = zeros(length(angle)) .+ 5 #m/s
    prob = zeros(length(angle)) .+ 1.0/length(angle) #probabilities that the wind will be that direction and at that speed

    in = 4.0
    out = 25.0
    ratedSpeed = 16.0
    ratedPower = 2.0E6
    efficiency = 0.944

    # construct object of turbine properties (diameter and power properties)
    turbine_properties = turbine(d,in,out,ratedSpeed,ratedPower,efficiency)

    # construct object of wake properties and constants
    model_properties = model(alpha)

    # construct object of wind properties (speeds, directions, percent)
    wind_properties = wind(speed,angle,prob)

    return turbine_properties,model_properties,wind_properties
end

function buildWakePlots(turbine_properties,model_properties,wind_properties,turbine_x,turbine_y)
    n = length(turbine_x)
    angle = wind_properties.directions[1]
    cosWind = cosd(-angle)
    sinWind = sind(-angle)
    plotPoints = zeros(n,4,2)
    wakeLength = 300
    r = turbine_properties.diameter/2
    alpha = model_properties.alpha
    for i = 1:n
        center_x = turbine_x[i]
        center_y = turbine_y[i]

        plotPoints[i,1,1] = center_x - (r+alpha*wakeLength) #left wake x
        plotPoints[i,1,2] = center_y - wakeLength #left wake y
        plotPoints[i,2,1] = center_x - r #left turbine edge x
        plotPoints[i,2,2] = center_y #left turbine edge y
        plotPoints[i,3,1] = center_x + r #right turbine edge x
        plotPoints[i,3,2] = center_y #right turbine edge y
        plotPoints[i,4,1] = center_x + (r+alpha*wakeLength) #right wake x
        plotPoints[i,4,2] = center_y - wakeLength #right wake y

        # rotate
        for j = 1:4
            x = plotPoints[i,j,1] - center_x
            y = plotPoints[i,j,2] - center_y

            plotPoints[i,j,1] = x*cosWind - y*sinWind + center_x
            plotPoints[i,j,2] = x*sinWind + y*cosWind + center_y
        end
    end

    return plotPoints
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

# turbine_x = [-50,0]
# turbine_y = [0,0]

power,energy,plots = totalPower(turbine_x,turbine_y)

# Plot
if !isempty(plots)
    p = plot(plots[1,:,1],plots[1,:,2],xlim=(-110,110),ylim=(-110,110),legend=false)
    # p = plot(plots[1,:,1],plots[1,:,2],xlim=(-55,10),ylim=(-25,25),legend=false)

    for i = 2:size(plots,1)
        plot!(p,plots[i,:,1],plots[i,:,2])
    end
    display(p)
    println(energy)
else
    plot(allAngles,energy,legend=false)
end