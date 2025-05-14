#= This code uses the digitized data from digitizePlot.jl (Stroh and Frasunkiewicz, 2025; https://github.com/AnStroh/PlotDigitizer.jl) (within this Pkg) to calculate the
equation of the reaction line of two phases from a X-T phase diagram. The results are the
coefficients a, b, c of X(T) = c + b*T + a*T², which are used in Chemical_Stefan_problem_XT.jl.
=#
using FileIO, DelimitedFiles, Plots

"""
    CalculateReactionLine(data1, data2)

Calculates coefficients for two transition lines in a phase diagram using Linear Least Squares.
The composition X is a function of temperature T, where X(T) = a + b*T + c*T².

##Units
- `T` is in [K].
- `X` is dimensionless (normalized).

# Arguments
- `data1::Array{Float64,2}`: The data points for the lower transition line. Each row represents a data point, where the first column is the x-coordinate and the second column is the y-coordinate.
- `data2::Array{Float64,2}`: The data points for the upper transition line. Each row represents a data point, where the first column is the x-coordinate and the second column is the y-coordinate.

# Returns
- `coeff_up::Array{Float64,1}`: The coefficients for the upper transition line equation.
- `coeff_do::Array{Float64,1}`: The coefficients for the lower transition line equation.
"""
function CalculateReactionLine(data1,data2)
    #Data preparation-------------------------------------
    T_data_up = data2[:,2]	                                            #Extract coordinates for linear least squares
    T_data_do = data1[:,2]                                              #Extract coordinates for linear least squares
    X_data_up = data2[:,1]                                              #Extract coordinates for linear least squares
    X_data_do = data1[:,1]                                              #Extract coordinates for linear least squares
    A1        = ones(length(X_data_up),3)
    A2        = ones(length(X_data_do),3)
    for (ni,_) in enumerate(1:length(X_data_up))
        A1[ni,2] = T_data_up[ni]
        A1[ni,3] = T_data_up[ni] ^ 2.0
    end
    for (ni,_) in enumerate(1:length(X_data_do))
        A2[ni,2] = T_data_do[ni]
        A2[ni,3] = T_data_do[ni] ^ 2.0
    end
    #Linear Least Squares---------------------------------
    coeff_up = A1 \ X_data_up	                                        #Coefficients for X(T) upper transition line
    coeff_do = A2 \ X_data_do                                           #Coefficients for X(T) lower transition line
    println("Coefficients for X1(T) = a1 + b1*T + c1*T², where T is in K")
    println("These are the coefficients for the lower line.")
    println("a1 = $(coeff_up[1])")
    println("b1 = $(coeff_up[2])")
    println("c1 = $(coeff_up[3])"); println()
    println("Coefficients for X2(T) = a2 + b2*T + c2*T², where T is in K.")
    println("These are the coefficients for the upper line.")
    println("a2 = $(coeff_do[1])")
    println("b2 = $(coeff_do[2])")
    println("c2 = $(coeff_do[3])");println()
    #Exporting data---------------------------------------
    folder_path = "examples/Examples_phase_diagram/"
    writedlm(folder_path * "Coefficients_Reaction_lines_test.csv", [coeff_up coeff_do ])

    return coeff_up, coeff_do
end

#Call function--------------------------------------------
#CAUTION: CSV.file with more points give inaccurate results at the boundaries (Check ~ X = 0 and ~ X = 1).
#file_name_1 = "examples/Examples_phase_diagram/Lower_Line.csv"
#file_name_2 = "examples/Examples_phase_diagram/Upper_Line.csv"

#CAUTION: CSV.file with fewer points give more accurate results at the boundaries (Check ~ X = 0 and ~ X = 1).
file_name_1 = "examples/Examples_phase_diagram/digitized_data_2.csv"    #Lower Line
file_name_2 = "examples/Examples_phase_diagram/digitized_data_1.csv"    #Upper Line

#Read data from CSV file
data1       = readdlm(file_name_1)
data2       = readdlm(file_name_2)
#Calculate coefficients for the reaction lines
coeff_up, coeff_do  = CalculateReactionLine(data1,data2)

#Check for correctness of the result within the plot
println("CAUTION: Be careful, that limits are captured correctly. Sometimes fewer points give better results.")
T  = LinRange(1273.0,1873.0,100)                                        #Choose the temperature range of the phase diagram
y  = coeff_do[1] .+ coeff_do[2] .* T .+ coeff_do[3] * T .^ 2.0
y2 = coeff_up[1] .+ coeff_up[2] .* T .+ coeff_up[3] * T .^ 2.0
plot(T,y)
plot!(T,y2)

T_data_up = data2[:,2]	                                            #Extract coordinates for linear least squares
T_data_do = data1[:,2]                                              #Extract coordinates for linear least squares
X_data_up = data2[:,1]                                              #Extract coordinates for linear least squares
X_data_do = data1[:,1]                                              #Extract coordinates for linear least squares
plot!(T_data_do,X_data_do,seriestype = :scatter, label = "Lower Line", color = :blue)
plot!(T_data_up,X_data_up,seriestype = :scatter, label = "Upper Line", color = :red)