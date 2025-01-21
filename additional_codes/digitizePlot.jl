#= 
This is a short summery of details related to our Olivine example and the digitalization of the
phase diagram.We are using the repository PlotDigitizer.jl in GitHub from Stroh & Frasunkiewicz
(07.01.2025,https://github.com/AnStroh/PlotDigitizer)
Their example is the same as in our case. The phase diagram was created using Perple_X (Connolly,
2009, DOI: 10:Q10014 DOI:10.1029/2009GC002540) and 

solution models/
thermodynamic data set as

=#

#Our input-----------------------------------------------------------------
#=
Reminder: PlotDigitizer.jl is not included into this Pkg. Please download PlotDigitizer.jl to 
use it or use your prefered Digitizer. To proceed 2 output files containing X-T values for the 
reaction lines are needed.
CAUTION: Be careful, that limits are captured correctly. Sometimes fewer points give better results.
=#
file_name = "Examples_phase_diagram/Ol_Phase_diagram_without_framework.png" 
T_BC      = (1273.0, 1873.0)   #min max of T in the phase diagram
X_BC      = (0.0, 1.0)         #min max of X (composition) in the phase diagram
println("Please export your data before finishing.")
lines = digitizePlot(X_BC, T_BC, file_name)


