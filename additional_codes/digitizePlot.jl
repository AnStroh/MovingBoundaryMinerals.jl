#= 
This is a short summery of details related to our Olivine example and the digitalization of the
phase diagram.We are using the repository PlotDigitizer.jl in GitHub from Stroh & Frasunkiewicz
(07.01.2025,https://github.com/AnStroh/PlotDigitizer)
Their example is the same as in our case. The phase diagram was created using Perple_X PERPLE_X software package (Connolly, 2005,
https://doi.org/10.1016/j.epsl.2005.04.033; Connolly, 2009, https://doi.org/10.1029/2009GC002540; Xiang & Connolly, 2022,
https://doi.org/10.1111/jmg.12626). 
thermodynamic data set: 
hp02ver.dat (Holland & Powell, 1998, https://doi.org/10.1111/j.1525-1314.1998.00140.x)
solution models: 
O(HP)       (Holland & Powell, 1998, https://doi.org/10.1111/j.1525-1314.1998.00140.x)   
melt(HP)    (Holland & Powell, 2001, https://doi.org/10.1093/petrology/42.4.673; White, Powell & Holland, 2001, https://doi.org/10.1046/j.0263-4929.2000.00303.x)  

Reminder: PlotDigitizer.jl is not included into this Pkg. Please download PlotDigitizer.jl to 
use it or use your prefered digitizer. The code shown here is only a previous version of the package. 
More options are possible with the help of the package. To proceed, two output files containing X-T values for the 
reaction lines are needed.
CAUTION: Be careful, that limits are captured correctly. Sometimes fewer points give better results.
=#


using GLMakie, FileIO, DelimitedFiles


"""
    calc_X_Y(point, X_BC, Y_BC, pixel)

Calculate the X and Y coordinates based on the given point, boundary conditions, and pixel information.

# Arguments
- `point::Tuple{Float64, Float64}`: The coordinates of the point to be transformed.
- `X_BC::Tuple{Float64, Float64}`: The boundary conditions for the X coordinate.
- `Y_BC::Tuple{Float64, Float64}`: The boundary conditions for the T coordinate.
- `pixel::Tuple{IntThe code shown here is only a previous version of the package. More options are possible with the help of the package Int}`: The pixel dimensions of the image.

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
- `X::Vector{T}`: A vector where `X[1]` is used as the y-intercept (b in the equation).
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
    digitizePlot(X_BC::Tuple, Y_BC::Tuple, file_name::String, export_name::String = "digitized_data")

Digitizes points from an image file and exports the coordinates.

# Arguments
- `X_BC::Tuple`: Tuple containing the X-axis boundary conditions.
- `Y_BC::Tuple`: Tuple containing the Y-axis boundary conditions.
- `file_name::String`: Path to the image file to be digitized.
- `export_name::String`: (Optional) Base name for the exported CSV files. Default is "digitized_data".

# Description
This function loads an image and allows the user to interactively digitize points on the image. The digitized points 
are transformed based on the provided boundary conditions and can be exported to CSV files.

# Interactions
- Left Mouse Button: Add a new point.
- Right Mouse Button: Finish the digitization process.
- Keyboard `d`: Delete the last point.
- Keyboard `n`: Add a new line.
- Keyboard `s`: Switch to the next line.
- Keyboard `e`: Export the current line to a CSV file.

# Returns
- `lines_conv[]`: A list of digitized lines with their coordinates transformed based on the boundary conditions.
"""

function digitizePlot(X_BC::Tuple, Y_BC::Tuple, file_name::String, export_name::String = "digitized_data")
    #Load image---------------------------------------------------------
    img         = rotr90(load(file_name))
    pixel       = size(img)
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
    #Plot the points and lines------------------------------------------
    p1          = lines!(scene, points, color = colors[active], linewidth = 3) 
    p2          = scatter!(scene, points, color = :red, marker = :cross,markersize = 25)
    #Define interactions------------------------------------------------
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
                        new_points      = []
                        new_points_conv = []
                        new_color       = RGBf(rand(3)...)
                        push!(lines[], new_points)
                        push!(lines_conv[], new_points_conv)
                        push!(colors, new_color)
                        println("------------------------")
                        println("A new line has been added.")
                        println("Number of lines: $(length(lines[]))")
                        println("------------------------")
                    #Switch line
                    elseif Keyboard.s in events(scene).keyboardstate
                        #Check if there is only one line----------------
                        if length(lines[]) == 1
                            println("------------------------")
                            println("There is only one line. To add a new line, press 'n' and click the left mouse button.")
                            println("------------------------")
                            return
                        end
                        #Switch line------------------------------------
                        lines[][active]         = copy(points[])
                        lines_conv[][active]    = copy(points_conv[])
                        active                  = active + 1
                        #Start from the first line if the last line is reached
                        if active > length(lines[])
                            active = 1
                        end
                        #Update the points in active line---------------
                        points[]        = []
                        points_conv[]   = []
                        push!(points[] , lines[][active]...)
                        push!(points_conv[] , lines_conv[][active]...)
                        p1.color        = colors[active]
                        #Notify the changes-----------------------------
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
                        mp          = events(scene).mouseposition[]
                        pos         = to_world(scene, mp)
                        pos_conv    = calc_X_Y(pos, X_BC, Y_BC, pixel)
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
            #Condition to stop digitalization---------------------------
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

file_name = "examples/Examples_phase_diagram/Ol_Phase_diagram_without_framework.png" 
T_BC      = (1273.0, 1873.0)                                                            #min max of T in the phase diagram
X_BC      = (0.0, 1.0)                                                                  #min max of X (composition) in the phase diagram
println("Please export your data before finishing.")
lines = digitizePlot(X_BC, T_BC, file_name)


