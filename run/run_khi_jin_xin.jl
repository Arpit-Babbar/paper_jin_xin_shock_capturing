using Tenkai
Eq = Tenkai.EqEuler2D
using Tenkai.TenkaicRK
using StaticArrays
EqJinXin = Tenkai.TenkaicRK.EqJinXin2D
#------------------------------------------------------------------------------
xmin, xmax = -1.0, 1.0
ymin, ymax = -1.0, 1.0

boundary_condition = (periodic, periodic, periodic, periodic)
γ = 1.4
equation_euler = Eq.get_equation(γ)

nx, ny = 32, 32

epsilon_relaxation = 1e-12

equation_jin_xin = EqJinXin.get_equation(equation_euler, epsilon_relaxation, nx, ny,
                                         thresholds = (1e-12, 1.0e-6),
                                        #  thresholds = (1e-14, 0.5e-6), Doesn't work
                                         jin_xin_dt_scaling = 0.5)

function inital_data_khi_chan(x, y)
    # change discontinuity to tanh
    # typical resolution 128^2, 256^2
    # domain size is [-1,+1]^2
    RealT = typeof(x)
    slope = 15
    B = tanh(slope * y + 7.5f0) - tanh(slope * y - 7.5f0)
    rho = 0.5f0 + 0.75f0 * B
    v1 = 0.5f0 * (B - 1)
    v2 = convert(RealT, 0.1) * sinpi(2 * x)
    p = 1.0
    gamma = 1.4
    return SVector(rho, rho * v1, rho * v2,
                   p / (gamma - 1.0) + 0.5 * (rho * v1 * v1 + rho * v2 * v2))
end

initial_value = EqJinXin.JinXinICBC(inital_data_khi_chan, equation_jin_xin)
exact_solution = EqJinXin.JinXinICBC((x,y,t) -> inital_data_khi_chan(x,y), equation_jin_xin)
boundary_value = exact_solution

degree = 3
solver = cBPR343()
solution_points = "gl"
correction_function = "radau"
numerical_flux = EqJinXin.rusanov

bound_limit = "no"
bflux = extrapolate
final_time = 10.0

cfl = 0.0
bounds = ([-Inf], [Inf]) # Not used in Euler
save_iter_interval = 0
save_time_interval = final_time / 20.0
animate = true # Factor on save_iter_interval or save_time_interval
compute_error_interval = 0

cfl_safety_factor = 0.5

#------------------------------------------------------------------------------
grid_size = [nx, ny]
domain = [xmin, xmax, ymin, ymax]
problem = Problem(domain, initial_value, boundary_value, boundary_condition,
                  final_time, exact_solution, source_terms = EqJinXin.jin_xin_source)
# limiter = setup_limiter_tvb(equation; tvbM = tvbM)
limiter = setup_limiter_blend(blend_type = fo_blend(equation_jin_xin),
                              indicating_variables = EqJinXin.rho_p_indicator!,
                              reconstruction_variables = conservative_reconstruction,
                              indicator_model = "gassner",
                              amax = 0.002,
                              debug_blend = false)
# limiter = setup_limiter_none()
scheme = Scheme(solver, degree, solution_points, correction_function,
                numerical_flux, bound_limit, limiter, bflux)
param = Parameters(grid_size, cfl, bounds, save_iter_interval,
                   save_time_interval, compute_error_interval,
                   animate = animate,
                   cfl_safety_factor = cfl_safety_factor,
                   saveto = joinpath(@__DIR__, "output_khi_chan2022_jin_xin_nx$nx"))
#------------------------------------------------------------------------------
sol = Tenkai.solve(equation_jin_xin, problem, scheme, param);

println(sol["errors"])

return sol;
