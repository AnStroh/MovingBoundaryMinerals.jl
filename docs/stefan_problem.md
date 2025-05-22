# The (chemical) Stefan problem

The Stefan problem describes the movement of a reaction front in a thermodynamically constrained problem (e.g., the propagation of an ice front or the crystallization/resorption of minerals). We included the option to calculate the growth/dissolution of crystals based on a thermodynmically constained data set (phase diagram). Here, we want to include a short describtion/additional information to the code.

## Input parameters

Most of the inpu parameters are self-explaning. However, we want to mention, that this section we explain some more details to the single input parameters. The following input parameters are required: the initial position of the interface, which also determines the size of the left phase in the composition profile, the composition of the assemblage of interest CompInt. For the calculations, information on the start and end temperature (Tstart, Tstop) of the growth process as well as maximum and minimum temperatures (TMIN, TMAX) of the chosen section of the phase diagram are required. TMIN and TMAX are calculation parameters and define a temperature range in which the model can be used. They are not to be confused with the Tstart and Tstop in the model, whereby Tstar,Tstop ∈ [TMIN,TMAX]. All temperatures must be given in K. Coefficients for the phase transformation lines can be stored in the parameter eq_values. CompInt stores the composition of interest of the solid solution and therefore determines the start composition. 

## Determing the interface composition
By digitizing two adjacent reaction lines of the phase diagram, we create a binary phase diagram to which we can apply our code. The reaction lines are described by two polynomials of the second degree (Eq. 1). 
X(T)=a*T^2+b*T+c              (1)
We use these polynomials to determine the compositions of the two phases as a function of temperature. These concentrations are linked to the composition of the solid solution/assemblage via the Lever rule. The initial composition profile is a step function with two homogeneous parts referring to the two materials. The respective compositions are based on the data from the phase diagram for the composition of CompInt at the starting temperature Tstart. Using already given information, we calculate the length of the whole modelling domain (Eqs. 2 and 3). 

