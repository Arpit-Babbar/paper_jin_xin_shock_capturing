include(joinpath(@__DIR__, "base_plotting.jl"))

rcParams = plt.PyDict(plt.matplotlib["rcParams"]) # Otherwise, the dictionary will be copied
rcParams["text.usetex"] = true # This will use LaTeX fonts (slow)
rcParams["font.size"] = 14.65
rcParams["axes.titlesize"] = 16
rcParams["axes.labelsize"] = 16

function plot_sol(ax, u, degree; lw = 2, label = "Numerical", ls = "--", color = "red", y_index = 2)
    nu = max(2, degree + 1)
    nx = Int(size(u, 1) / nu)
    ax.plot(u[1:nu, 1], u[1:nu, y_index], fillstyle = "none",
            color = color, linewidth = lw, label = label, ls = ls)
    for ix in 2:nx
        i1 = (ix - 1) * nu + 1
        i2 = ix * nu
        ax.plot(u[i1:i2, 1], u[i1:i2, y_index], fillstyle = "none",
                color = color, linewidth = lw, ls = ls)
    end
end

function plot_avg(ax, u; lw = 2, label = "Numerical", ls = "--", color = "red", y_index = 2,
                  type = :line)
    if type == :line
        ax.plot(u[:, 1], u[:, y_index], fillstyle = "none",
                color = color, linewidth = lw, label = label, ls = ls)
    elseif type == :scatter
        ax.scatter(u[:, 1], u[:, y_index],  color = color, label = label)
    else
        @assert false "Unknown plot type: $type"
    end
end

# Burgers' equation solution
fig, ax = plt.subplots()
ax.set_xlabel("\$x\$")
ax.set_ylabel("\$ u \$")
ax.grid(true, linestyle = "--")
jin_xin = readdlm(joinpath(@__DIR__, "..", "run", "burg1d_jinxin_nx20", "sol.txt"))
blend = readdlm(joinpath(@__DIR__, "..", "run", "burg1d_blend_nx20", "sol.txt"))
reference = readdlm(joinpath(@__DIR__, "..", "run", "burg1d_blend_nx2000", "sol.txt"), skipstart = 1)
# plot_sol(ax, reference, 2, ls = "-", label = "Reference", color = "k")
ax.plot(reference[:, 1], reference[:, 2], ls = "-", label = "Reference", color = "k")
plot_sol(ax, jin_xin, 3, ls = "--", label = "Jin-Xin", color = "red")
plot_sol(ax, blend, 3, ls = "-.", label = "Blend", color = "blue")
ax.legend()
filename = joinpath("burg1d.pdf")
fig
fig.savefig(filename)
run(`bash pdfbb $filename`)

# Bucklev's equation solution
fig, ax = plt.subplots()
exact_solution = x -> Tenkai.EqBuckleyLeverett1D.hatbuck_exact_a025(x, 0.15)
ax.set_xlabel("\$x\$")
ax.set_ylabel("\$ u \$")
ax.grid(true, linestyle = "--")
jin_xin = readdlm(joinpath(@__DIR__, "..", "run", "bucklev_jin_xin_nx50", "sol.txt"))
blend = readdlm(joinpath(@__DIR__, "..", "run", "bucklev_blend_nx50", "sol.txt"))
# reference = readdlm(joinpath(@__DIR__, "..", "run", "bucklev_blend_nx2000", "sol.txt"), skipstart = 1)
# plot_sol(ax, reference, 2, ls = "-", label = "Reference", color = "k")

domain = LinRange(-1.0, 1.0, 1000)
ax.plot(domain, exact_solution.(domain), ls = "-", label = "Exact", color = "k")
plot_sol(ax, jin_xin, 3, label = "Jin-Xin", color = "red", ls = "--")
plot_sol(ax, blend, 3, label = "Blend", color = "blue", ls = "-.")
ax.legend()
filename = joinpath("bucklev.pdf")
fig
fig.savefig(filename)
run(`bash pdfbb $filename`)

# Jin Xin Burgers' equation solution

epsilon_relaxation_array = (1e-1, 1e-2, 1e-3, 1e-4, 1e-12)
eps2string = Dict(1e-1 => "1", 1e-2 => "2", 1e-3 => "3", 1e-4 => "4", 1e-12 => "12")
burger_sols_array = [readdlm(joinpath(@__DIR__, "..", "run", "jin_xin", "jin_xin_nx20_eps$(eps2string[epsilon_relaxation])", "sol.txt")) for epsilon_relaxation in epsilon_relaxation_array]

using Tenkai.EqBurg1D: find_zero
function initial_value_burg_marco(x)
    k = 1.0
    u = 2.0 + sinpi(k * (x[1] - 0.7))
    return u
end

function exact_solution_burger_marco(x, t_)
    t = min(0.3, t_)
    implicit_eqn(u) = u - initial_value_burg_marco(x - t * u)
    seed = initial_value_burg_marco(x)
    value = find_zero(implicit_eqn, seed)
    return value
end

exact(x) = exact_solution_burger_marco(x, 0.25)
for i in 1:5
    burgers_sol = burger_sols_array[i]
    global fig, ax
    fig, ax = plt.subplots()
    ax.set_xlabel("\$x\$", fontsize = 22)
    ax.set_ylabel("\$ u \$", fontsize = 22)
    ax.tick_params(axis="both", which="major", labelsize=23)
    ax.tick_params(axis="both", which="minor", labelsize=23)
    ax.grid(true, linestyle = "--")
    ax.plot(burgers_sol[:, 1], exact.(burgers_sol[:, 1]), label = "Exact", c = "k", lw = 2,
            ls = "-")
    epsilon_string = eps2string[epsilon_relaxation_array[i]]
    ax.set_title("\$ \\varepsilon = 10^{-$epsilon_string} \$", fontsize = 25)
    # plot_sol(ax, burgers_sol, 3)
    ax.plot(burgers_sol[:, 1], burgers_sol[:, 2], label = "Numerical", color = "red",
            lw = 2, ls = "--")
    ax.legend()
    epsilon = epsilon_relaxation_array[i]
    global filename = joinpath("jin_xin_burgers_eps$(epsilon_string).pdf")
    fig.savefig(filename)
    run(`bash pdfbb $filename`)
end

# Blast

# Density
fig, ax = plt.subplots()
ax.set_xlabel("\$x\$", fontsize = 16)
ax.set_ylabel("\$ \\rho \$", fontsize = 16)
ax.tick_params(axis="both", which="major", labelsize=18)
ax.tick_params(axis="both", which="minor", labelsize=18)
ax.grid(true, linestyle = "--")
blast_jin_xin = readdlm(joinpath(@__DIR__, "..", "run", "blast_jin_xin_nx400/avg.txt"))
blast_blend = readdlm(joinpath(@__DIR__, "..", "run", "blast_blend_nx400/avg.txt"))
exact_data = Tenkai.EqEuler1D.exact_solution_data("blast")
# approx_2000 = readdlm(joinpath(@__DIR__, "..", "run", "tenmom_super_stiff_1d/tmp_super_stiff_nx2000/sol.txt"))
# ref = readdlm(joinpath(@__DIR__, "..", "run", "tenmom_super_stiff_1d/ref_nx100000/avg.txt"))
# ax.plot(ref[:, 1], ref[:, 2], label = "Reference", lw = 2, c = "k")
plot_sol(ax, exact_data, 0, ls = "-", label = "Reference", color = "k")
plot_sol(ax, blast_jin_xin, 3, ls = "--", label = "Jin-Xin", color = "red", y_index = 2)
plot_sol(ax, blast_blend, 3, ls = "-.", color = "blue", label = "Blend", y_index = 2)
ax.set_xlim(0.55, 0.9)
ax.legend()
filename = joinpath("blast_density.pdf")
fig
fig.savefig(filename)
run(`bash pdfbb $filename`)

# Pressure
fig, ax = plt.subplots()
ax.set_xlabel("\$x\$", fontsize = 16)
ax.set_ylabel("\$ p \$", fontsize = 16)
ax.tick_params(axis="both", which="major", labelsize=18)
ax.tick_params(axis="both", which="minor", labelsize=18)
ax.grid(true, linestyle = "--")
blast_jin_xin = readdlm(joinpath(@__DIR__, "..", "run", "blast_jin_xin_nx400/avg.txt"))
blast_blend = readdlm(joinpath(@__DIR__, "..", "run", "blast_blend_nx400/avg.txt"))
exact_data = Tenkai.EqEuler1D.exact_solution_data("blast")
# approx_2000 = readdlm(joinpath(@__DIR__, "..", "run", "tenmom_super_stiff_1d/tmp_super_stiff_nx2000/sol.txt"))
# ref = readdlm(joinpath(@__DIR__, "..", "run", "tenmom_super_stiff_1d/ref_nx100000/avg.txt"))
# ax.plot(ref[:, 1], ref[:, 2], label = "Reference", lw = 2, c = "k")
plot_sol(ax, exact_data, 0, ls = "-", label = "Reference", color = "k", y_index = 4)
plot_sol(ax, blast_jin_xin, 3, ls = "--", label = "Jin-Xin", color = "red", y_index = 4)
plot_sol(ax, blast_blend, 3, ls = "-.", color = "blue", label = "Blend", y_index = 4)
plt.xlim(0.55, 0.9)
ax.legend()
filename = joinpath("blast_pressure.pdf")
fig
fig.savefig(filename)
run(`bash pdfbb $filename`)
