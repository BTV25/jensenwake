# Wind direction from north is 0 degrees
# Wind direction from east is 90 degrees

# Each turbine has a diameter, x coordinate, and y coordinate
# Each column in the array turbines is a turbine
# Row 1 is the x coordinate
# Row 2 is the y coordinate
# Row 3 is the diameter of the turbine
# k and Ct are treated as constants across all turbines
# All heights are assumed to be the same

# For this model, derived from Katic 1986 "A simple model for cluster efficiency"
# each turbine is treated as a disk and therefore can be affected by partial wakes
# meaning a turbine can experience Wake A on part of it's surface and by Wake B on
# a different part of its surface and by both the wakes on the rest of its surface

# Turbines are assumed to spin to face the wind

###################################################################################
# Functions
# To account for different wind directions this function rotates the coordinate
# frame of the turbines to line them with the wind
function rotateFrame(theta,turbines,n)
    R = [cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; 0 0 1]
    temp = zeros(3,n)
    for i = 1:n
        temp[:,i] = R*[turbines[1:2,i];0]
        temp[3,i] = turbines[3,i]
    end
    return temp
end

# Resets the coordinate frame to give positions of turbines
function resetTurbineFrame(theta,turbines,n)
    theta = -theta
    R = [cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; 0 0 1]
    temp = zeros(3,n)
    for i = 1:n
        temp[:,i] = R*[turbines[1:2,i];0]
    end
    return temp[1:2,:]
end

# function to calculate the wake deficit (1-v/u) from a single turbine
function wakeDeficit(Ct,k,x,D)
    return (1-sqrt(1-Ct))/((1+2*k*x/D)^2)
end

# fuction to calculate wake deficit of multiple combined wakes (pass in array of deficits)
function wakeInteractionDeficit(deficit)
    total = 0
    for def in deficit
        total += def^2
    end
    return sqrt(total)
end

# d is in terms of the diameter (0.5 D)
# calculates the area of a slice that goes to the edges of the turbine
function calculateAreaOfTurbineSegment(d)
    theta = 2*acos(1-2*d)
    return (0.5*(theta-sin(theta)))/pi
end

# d1 and d2 are in terms of the diameter d2 > d1
# uses multiple percent areas of the turbine given by
# calculateAreaOfTurbineSegment to calculate percent areas of slices
# that dont extend to the edges
function calculatePercentAreaOfTurbineSlice(d1, d2)
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

function calculateDeficitOnSingleTurbine(turbines,Ct,k,numTurbine,n)
    rightTurbineEdge = turbines[1,numTurbine] + .5*turbines[3,numTurbine]
    leftTurbineEdge = turbines[1,numTurbine] - .5*turbines[3,numTurbine]
    
    # wakes: Row 1 is where the left edge of the wake hits the turbine (.3 D)
    # row 2 is when the right edge of the wake hits the turbine
    # row 3 is the deficit of the single wake
    wakes = zeros(3,n)
    for i = 1:n
        if i == numTurbine # a turbine is not affected by its own wake
            continue
        end
        if turbines[2,i] <= turbines[2,numTurbine] # ignore turbines behind current turbine
            continue
        end
        # calculate wake range at the distance from the turbine i
        ranges = turbines[3,i]+k*(abs(turbines[2,i] - turbines[2,numTurbine]))
        rightWakeEdge = turbines[1,i] + ranges
        leftWakeEdge = turbines[1,i] - ranges
        
        # if numTurbine falls in range determine which sections are affected 
        if (leftTurbineEdge >= rightWakeEdge) || (rightTurbineEdge <= leftWakeEdge)
            # wake does not affect the turbine
            continue 
        end

        # The wake must affect the turbine, determine the center point of that effect and the radius
        right = false
        left = false
        if rightWakeEdge <= rightTurbineEdge # right edge of wake hits turbine
            right = true
        end
        if leftWakeEdge >= leftTurbineEdge # left edge of wake hits turbine
            left = true
        end

        # calculate the deficit of this wake
        wakes[3,i] = wakeDeficit(Ct,k,abs(turbines[2,i] - turbines[2,numTurbine]),turbines[3,i])

        # define wake effect area with center and radius relative to the diameter
        if right && left
            # find where both edges of the wake hit the turbine relative to turbine's left edge
            # in terms of the diameter
            leftHit = (leftWakeEdge - leftTurbineEdge)/turbines[3,numTurbine]
            rightHit = (rightWakeEdge - leftTurbineEdge)/turbines[3,numTurbine]
        elseif right
            # find where right edge of wake hits the turbine relative to left edge
            leftHit = 0
            rightHit = (rightWakeEdge - leftTurbineEdge)/turbines[3,numTurbine]
        elseif left
            # find where left edge of wake hits the turbine relative to the right edge
            leftHit = (leftWakeEdge - leftTurbineEdge)/turbines[3,numTurbine]
            rightHit = 1
        elseif !right && !left
            leftHit = 0
            rightHit = 1
        end

        # store wake edges in wakes
        wakes[1,i] = leftHit
        wakes[2,i] = rightHit
    end
    # remove non-affecting wakes
    remove = []
    for i = 1:n
        if wakes[2,i] == 0
            remove = vcat(remove,i)
        end
    end
    wakes = wakes[:,setdiff(1:end, tuple(remove...))]
    if isempty(wakes)
        return 0
    end
    
    # determine where wakes overlap and determine how much area is affected by
    # different single wakes and wake combinations
    order = zeros(size(wakes,1)*size(wakes,2),3)
    for i = 1:size(wakes,1)
        f = (i-1)*size(wakes,2)+1
        l = f + size(wakes,2)-1
        order[f:l,1] .= i
        order[f:l,2] .= 1:size(wakes,2)
    end

    for i = 1:size(order,1)
        row = convert(Int32,order[i,1])
        col = convert(Int32,order[i,2])
        order[i,3] = wakes[row,col]
    end

    order = sortslices(order,dims=1,by=x->x[3],rev=false)
    # order is a list ordered by the third column with columns 1 and 2 being
    # the row and column in wake at which the edge is located
    active = convert(Array{Int8},zeros(length(wakes),1))

    active[convert(Int32,order[1,2])] = 1
    d1 = order[1,3]

    totalDeficit = 0
    def = transpose(wakes[3,:])
    for i = 2:size(order,1)
        d2 = order[i,3] # what is the distance of this point
        currentDeficit = wakeInteractionDeficit(def .* active)
        totalDeficit += calculatePercentAreaOfTurbineSlice(d1,d2) * currentDeficit

        if order[i,1] == 1 # start of stop the wake
            active[convert(Int32,order[i,2])] = 1
        else
            active[convert(Int32,order[i,2])] = 0
        end
        
        d1 = d2
    end

    return totalDeficit
end

function calculateFarmDeficit(turbines,Ct,k,n)
    deficit = 0
    for i = 1:n
        deficit += calculateDeficitOnSingleTurbine(turbines,Ct,k,i,n)
    end
    return deficit
end
###################################################################################
# Coordinate area BOX (mins are zero)
xMax = 10
yMax = 10

# Create turbines
n = 6 #number of turbines
turbines = zeros(3,n)
for i = 1:n
    # Set the turbines in a straight line along the center of the space provided with 
    # equal diameter
    turbines[1,i] = xMax/(n-1)*(i-1)
    turbines[2,i] = yMax/2
    turbines[3,i] = 5
end

# Set constants
k = 0.1 #determines the angle at which the wake expands at
Ct = 0.9 #thrust coefficient of the turbine

# Declare wind direction from North is 0, from East is 90
windDirection = 90

# Rotate the frame
turbines = rotateFrame(windDirection,turbines,n)

# Calculate the deficit
farmDeficit = calculateFarmDeficit(turbines,Ct,k,n)
println(farmDeficit)

# Plot wakes?

# Reset the frames
turbines = resetTurbineFrame(windDirection,turbines,n)
# reset wake frame