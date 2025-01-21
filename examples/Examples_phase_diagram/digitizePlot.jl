using GLMakie, FileIO, DelimitedFiles


"""
    calc_X_Y(point, X_BC, Y_BC, pixel)

Calculate the X and Y coordinates based on the given point, boundary conditions, and pixel information.

# Arguments
- `point::Tuple{Float64, Float64}`: The coordinates of the point to be transformed.
- `X_BC::Tuple{Float64, Float64}`: The boundary conditions for the X coordinate.
- `Y_BC::Tuple{Float64, Float64}`: The boundary conditions for the T coordinate.
- `pixel::Tuple{Int, Int}`: The pixel dimensions of the image.

# Returns
- `Tuple{Float64, Float64}`: The transformed X and Y coordinates.
"""
function calc_X_Y(point, X_BC, Y_BC, pixel)
    m_X = (X_BC[2]- X_BC[1])/pixel[1]
    m_Y = (Y_BC[2]- Y_BC[1])/pixel[2]

    X = linear_func(X_BC, point[1], m_X)
    Y = linear_func(Y_BC, point[2], m_Y) 
    
    return (X, Y)
end

"""
    linear_func(X, x, m)

Compute the linear function `Y = m * x + b`.

# Arguments
- `X::Vector{T}`: A vector where `X[1]` is used as the y-intercept.
- `x::T`: The independent variable.
- `m::T`: The slope of the linear function.

# Returns
- `Y::T`: The computed value of the linear function.
"""
function linear_func(X, x, m)
    Y = m * x + X[1]
    return Y
end

"""
    digitizePlot(X_BC, Y_BC, file_name, export_name)

Digitize a plot by clicking on the points of interest and exporting the coordinates.

# Arguments
- `X_BC::Tuple{Float64, Float64}`: The boundary conditions for the X coordinate.
- `Y_BC::Tuple{Float64, Float64}`: The boundary conditions for the Y coordinate.
- `file_name::String`: The name of the file containing the plot.
- `export_name::String`: The name of the file to which the coordinates will be exported.

# Returns
- Vector{Vector} where each vector contains the coordinates of the points in the corresponding line.

"""

function digitizePlot(X_BC::Tuple, Y_BC::Tuple, file_name::String, export_name::String = "digitized_data")
    #Load image---------------------------------------------------------
    img   = rotr90(load(file_name))
    pixel = size(img)
    #Create scene-------------------------------------------------------
    scene       = Scene(camera = campixel!, size=pixel)
    image!(scene,img) 
    #Initialize Observables---------------------------------------------
    lines       = Observable([[]])
    lines_conv  = Observable([[]])
    points      = Observable(Point2f[])
    points_conv = Observable(Point2f[])
    colors      = [RGBf(rand(3)...)]
    active      = 1
    #Flag to stop the interaction---------------------------------------
    interaction = true                      
    #Plot the points and lines-------------------------------------------
    p1 = lines!(scene, points, color = colors[active], linewidth = 3) 
    p2 = scatter!(scene, points, color = :red, marker = '+',markersize = 25)
    #Define interactions-------------------------------------------------
    on(events(scene).mousebutton) do event
        while interaction == true
            if event.button == Mouse.left
                if event.action == Mouse.press
                    #Delete points
                    if Keyboard.d in events(scene).keyboardstate
                        deleteat!(points[], length(points[]))
                        deleteat!(points_conv[], length(points_conv[]))
                        notify(points)
                        notify(points_conv)
                        println("------------------------")
                        println("The last point has been deleted.")
                        println("------------------------")
                    #Add new line
                    elseif Keyboard.n in events(scene).keyboardstate
                        new_points = []
                        new_points_conv = []
                        new_color = RGBf(rand(3)...)
                        push!(lines[], new_points)
                        push!(lines_conv[], new_points_conv)
                        push!(colors, new_color)
                        println("------------------------")
                        println("A new line has been added.")
                        println("Number of lines: $(length(lines[]))")
                        println("------------------------")
                    #Switch line
                    elseif Keyboard.s in events(scene).keyboardstate
                        #Check if there is only one line--------------------------------
                        if length(lines[]) == 1
                            println("------------------------")
                            println("There is only one line. To add a new line, press 'n' and click the left mouse button.")
                            println("------------------------")
                            return
                        end
                        #Switch line-----------------------------------------------------
                        lines[][active] = copy(points[])
                        lines_conv[][active] = copy(points_conv[])
                        active = active + 1
                        #Start from the first line if the last line is reached-----------
                        if active > length(lines[])
                            active = 1
                        end
                        #Update the points in active line---------------------------------
                        points[] = []
                        points_conv[] = []
                        
                        push!(points[] , lines[][active]...)
                        push!(points_conv[] , lines_conv[][active]...)

                        p1.color = colors[active]
                        #Notify the changes------------------------------------------------
                        println("------------------------")
                        println("Switched to line $active.")
                        println("------------------------")
                        notify(points)
                        notify(lines)
                    #Export line
                    elseif Keyboard.e in events(scene).keyboardstate
                        println("------------------------")
                        println("Exporting line_$active.")
                        println("------------------------")
                        writedlm(export_name * "_$active.csv", points_conv[])
                    #Set new points
                    else
                        mp = events(scene).mouseposition[]
                        pos = to_world(scene, mp)
                        pos_conv = calc_X_Y(pos, X_BC, Y_BC, pixel)
                        println("------------------------")
                        println("Your new coordinate [x,y] within the picture is: $pos.")
                        println("\n X, Y = $pos_conv \n")
                        println("To exit the function, please press the right mouse button.")
                        println("------------------------")
                        push!(points[], mp)
                        push!(points_conv[], pos_conv)
                        notify(points)
                        notify(points_conv)
                    end
                end
            end
            #Condition to stop digitalization---------------------------------
            if event.button == Mouse.right
                if event.action == Mouse.press
                    println("You finished the digitalization.")
                    println("Please close the GLMakie window you used.")
                    break
                end
            end
            return points_conv
        end
    end
    
    display(scene)

    return lines_conv[]
end

file_name = "Examples_phase_diagram/Ol_Phase_diagram_without_framework.png"
# file_name = "Examples_phase_diagram/test.png"
X_BC       = (0.0,   1.0)   #min max of X in the phase diagram
Y_BC       = (1273.0, 1873.0)   #min max of Y in the phase diagram
lines      = digitizePlot(X_BC, Y_BC, file_name)