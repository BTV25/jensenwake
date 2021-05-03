using Test

# function to calculate the wake deficit (1-v/u) from a single turbine
function wakeDeficit(Ct,k,x,D)
    return (1-sqrt(1-Ct))/((1+2*k*x/D)^2)
end

# fuction to calculate wake deficit of combined wakes
function wakeInteractionDeficit(Ct,k,x,D)
    deficit = 0
    for i = 1:length(x)
        deficit += wakeDeficit(Ct,k,x[i],D)
    end
    return sqrt(deficit)
end

function determineIfTurbinesAreInWakes(turbines,k,D,n)
    # example inWake 1,2 means turbine 1 is in the wake of turbine 2
    inWake = zeros(Int8,n,n)
    # determine which turbines are in wake of a single turbine
    for i = 1:n # turbine causing the wake
        for j = 1:n # is this turbine in the wake?
            if i == j # turbine cant be in its own wake
                continue
            elseif turbines[1,i] >= turbines[1,j] #if turbine is in front in cant be in wake of those behind it
                continue
            else # turbine j is behind turbine i
                range = D+k*(turbines[1,j] - turbines[1,i])
                maxYWake = turbines[2,i] + range
                minYWake = turbines[2,i] - range

                if (turbines[2,j] < minYWake) || (turbines[2,j] > maxYWake)
                    continue
                end
            end
            inWake[i,j] = 1
        end
    end
    return inWake
end

function calculateTotalDeficitAtEachTurbine(turbines,inWake,Ct,k,D,n)
    deficit = zeros(Float64,n)
    for i = 1:n # do total calculation for each turbine
        x = zeros(Float64,n)
        for j = 1:n
            if inWake[j,i] == 1
                x[j] = abs(turbines[1,j]-turbines[1,i])
            else
                x[j] = 0;
            end
        end
        x = filter(a->a!=0,x)
        if isempty(x)
            deficit[i] = 0
        else
            deficit[i] = wakeInteractionDeficit(Ct,k,x,D)
        end
    end
    return deficit
end

# Deficit is calculated as if wind comes from the east, this function is help account for wind from 
# other directions by rotating the cordinate frame so the negative x is pointing into the wind
function rotateWind(theta,turbines,n)
    R = [cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; 0 0 1]
    temp = zeros(3,n)
    for i = 1:n
        temp[:,i] = R*[turbines[:,i];0]
    end
    return temp[1:2,:]
end

# test
@test round(wakeDeficit(0.4,0.1,40,20),digits=3) == 0.115
@test round(wakeDeficit(1,0.1,40,20),digits=4) == 0.5102

# Row 1 is x, Row 2 is y, each column is a turbine location
# treat each turbine as a point and wind will only come from one direction 
# n is number of turbines, xMax and yMax define the area
n = 4
xMax = 20
yMax = 20
turbines = zeros(Float64,2,n)

# set turbines in a line parallel to the x axis
for i = 1:n
    turbines[1,i] = xMax/n*i
    turbines[2,i] = yMax/2
end
originalTurbines = turbines

k = 0.1
D = 1
Ct = 0.2
theta = 0;

# rotate into wind to allow for deficit calculation (default wind comes from the west)
turbines = rotateWind(theta,turbines,n)

inWake = determineIfTurbinesAreInWakes(turbines,k,D,n)
deficit = calculateTotalDeficitAtEachTurbine(turbines,inWake,Ct,k,D,n)

# rotate away from wind to get positions 
turbines = rotateWind(-theta,turbines,n)
@test turbines == originalTurbines

println(deficit)
println(sum(deficit))