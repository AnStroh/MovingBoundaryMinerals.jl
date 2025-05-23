module MOBILE

using LinearAlgebra, BenchmarkTools, SpecialFunctions,Interpolations
using SparseArrays, Revise


include("Benchmarks.jl")
export analytical_sol_step_function, analytical_sol_half_space, calc_sinus_sphere, crank_time_transformation1, crank_time_transformation2, lasaga, rayleigh_fractionation ,smith,crank_time_transformation3

include("functions.jl")
export advect_interface_regrid!, blkdiag, calculate_dt, calc_mass_vol, calc_mass_vol_simple_diff, calc_volume, create_grid!, find_dt, fill_matrix!, linear_interpolation_1D, linspace_interface, preallocations, regrid!, set_inner_bc_mb!, set_inner_bc_flux!,set_inner_bc_Lasaga!, set_outer_bc!, trapezoidal_integration, update_time!, update_t_dependent_param!, update_t_dependent_param_simple!, construct_matrix_fem, solve_soe,calc_mass_err, make_dx_right, newton_solver, define_new_grid, sinusoid_profile

include("functionsStefan.jl")
export calculate_density, coeff_trans_line, composition, digitise_plot, ndgrid, set_inner_bc_stefan!,values_between_known_indices!

end
