# Digitizing the phase diagram and extracting coefficients

## Phase diagrams

The phase diagrams, which shows the composition of the system as a function of temperature, include all important information about chemical reactions between mineral assemblages. Therefore, it is necessary to describe the equilibrium conditions at the interface. Equations for the phase transformation lines contain all thermodynamically important data. These equations can be described by a quadratic equation (Eq. 1). Within our code, two transformation lines are always required to determine the composition of the components at the edges of a stable assemblage fields. The user needs to define coefficients (a, b, c) for both transformation lines as input prior to the start of the model.

## Work flow

We have used the PlotDigitizer.jl package in Julia (Stroh and Frasunkiewicz, 2025) (https://github.com/AnStroh/PlotDigitizer.jl, version v0.1.0) to digitize the phase diagram. A shortened version of the package can be found in the additional_codes folder under digitizePlot.jl. The package can be used to digitize reaction lines and export the X-T coordinates. The coordinates can then be used to calculate the coefficients for the reaction lines. We provide a short guideline here.
First, the user defines the pathway of the image of the phase diagram (.png or .jpg file) as a variable. It is important that the image only shows the phase diagram (like here: examples\Examples_phase_diagram\Ol_Phase_diagram_without_framework.png) and has no border. The phase diagram should show the composition on the x-axis and the temperature on the y-axis. Next, the user is asked for the limits of the axis, which should be defined in X_BC and Y_BC. Afterwards, the user can run the function digitizePlot(X_BC,Y_BC,file_name) to set or delete points of lines, switch the line and extract coordinates. Using CalculateReactionLine.jl, the coefficient of:
X(T)=a*T^2+b*T+c              (1)
where X is the composition the phase transition, T is the temperature in K and a, b and c depict the coefficients defining the quadratic equation, which can be calculated with the least squares method and subsequently stored. 
However, this step can be skipped, if the coefficients are known. In this case, the user can enter the coefficients directly in the main code.
