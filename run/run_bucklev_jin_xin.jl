import Roots.find_zero
using Tenkai
using Tenkai.TenkaicRK
# Submodules

Eq = Tenkai.EqBuckleyLeverett1D
EqJinXin = Tenkai.TenkaicRK.EqJinXin1D
# Set backend

#------------------------------------------------------------------------------
xmin, xmax = -1.0, 1.0
initial_value = Eq.hatbuck_iv
boundary_condition = (periodic, periodic)
exact_solution = Eq.hatbuck_exact_a025
final_time = 0.4 # Larger epsilon is needed for time > 0.15

epsilon_relaxation = 1e-4

eq_bucklev = Eq.get_equation()
A = () -> nothing

advection_jin_xin = (x, u, eq_bucklev) -> nothing

function advection_jin_xin_plus(ul, ur, F, eq_bucklev)
    nothing
end
function advection_jin_xin_minus(ul, ur, F, eq_bucklev)
    nothing
end

nx = 50
equation_jin_xin = EqJinXin.get_equation(eq_bucklev, advection_jin_xin,
                                         advection_jin_xin_plus,
                                         advection_jin_xin_minus, epsilon_relaxation, nx,
                                         thresholds = (1.5e-12, 1e-3),
                                         jin_xin_dt_scaling = 1.0)

# Is it really a struct?
initial_value_struct = EqJinXin.JinXinICBC(Eq.hatbuck_iv, equation_jin_xin)
initial_value = initial_value_struct
boundary_value_struct = EqJinXin.JinXinICBC(Eq.hatbuck_exact, equation_jin_xin)
boundary_value = boundary_value_struct
exact_solution = EqJinXin.JinXinICBC(Eq.hatbuck_exact_a025, equation_jin_xin)

degree = 3
solver = cSSP2IMEX433()
solution_points = "gl"
correction_function = "radau"
bflux = evaluate
numerical_flux = EqJinXin.rusanov
bound_limit = "no"

cfl = 0.0
bounds = ([0.0], [1.0])
tvbM = 0.0
save_iter_interval = 0
save_time_interval = 0.9
compute_error_interval = 0
cfl_safety_factor = 0.9
animate = true
#------------------------------------------------------------------------------
grid_size = nx
domain = [xmin, xmax]
problem = Problem(domain, initial_value, boundary_value, boundary_condition,
                  final_time, exact_solution, source_terms = EqJinXin.jin_xin_source)
limiter = Tenkai.setup_limiter_blend(blend_type = Tenkai.fo_blend(equation_jin_xin),
                                     indicating_variables = Tenkai.conservative_indicator!,
                                     reconstruction_variables = Tenkai.conservative_reconstruction,
                                     indicator_model = "gassner",
                                     debug_blend = false)
# limiter = setup_limiter_none()
scheme = Scheme(solver, degree, solution_points, correction_function,
                numerical_flux, bound_limit, limiter, bflux)
param = Parameters(grid_size, cfl, bounds, save_iter_interval,
                   save_time_interval, compute_error_interval,
                   cfl_safety_factor = cfl_safety_factor,
                   animate = animate, saveto = joinpath(@__DIR__, "bucklev_jin_xin_nx$nx"))
#------------------------------------------------------------------------------
sol = Tenkai.solve(equation_jin_xin, problem, scheme, param)

println(sol["errors"])

sol["plot_data"].p_ua
