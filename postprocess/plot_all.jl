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


fig, ax = plt.subplots()
ax.set_xlabel("\$x\$")
ax.set_ylabel("\$ h \$")
ax.grid(true, linestyle = "--")
approx = readdlm(joinpath(@__DIR__, "..", "run", "ssw_roll_1d", "ssw_roll_wave2", "sol.txt"))
exact_brock2 = readdlm(joinpath(@__DIR__, "exact_brock2.txt"), skipstart = 1)
h0 = 5.33e-3
ax.scatter(exact_brock2[:, 1] * 1.8, exact_brock2[:, 2] * h0, label = "Experimental", lw = 2,
           c = "k")
plot_sol(ax, approx, 1, ls = "--")
ax.legend()
filename = joinpath("roll_wave2.pdf")
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

fig, ax = plt.subplots()
ax.set_xlabel("\$x\$", fontsize = 16)
ax.set_ylabel("\$ v_1 \$", fontsize = 16)
ax.tick_params(axis="both", which="major", labelsize=18)
ax.tick_params(axis="both", which="minor", labelsize=18)
ax.grid(true, linestyle = "--")
approx_500 = readdlm(joinpath(@__DIR__, "..", "run", "tenmom_super_stiff_1d/tmp_super_stiff_nx500/sol.txt"))
approx_1000 = readdlm(joinpath(@__DIR__, "..", "run", "tenmom_super_stiff_1d/tmp_super_stiff_nx1000/sol.txt"))
# approx_2000 = readdlm(joinpath(@__DIR__, "..", "run", "tenmom_super_stiff_1d/tmp_super_stiff_nx2000/sol.txt"))
ref = readdlm(joinpath(@__DIR__, "..", "run", "tenmom_super_stiff_1d/ref_nx100000/avg.txt"))
ax.plot(ref[:, 1], ref[:, 3], label = "Reference", lw = 2, c = "k")
plot_sol(ax, approx_500, 3, ls = "dashed", label = "\$M=500\$", y_index = 3)
plot_sol(ax, approx_1000, 3, ls = "dotted", color = "blue", label = "\$M=1000\$", y_index = 3)
ax.legend()
filename = joinpath("tmp_super_stiff_v1.pdf")
fig
fig.savefig(filename)
run(`bash pdfbb $filename`)

# Burgers and super Burgers nonlinear source

## Burgers' test (shock)
fig, ax = plt.subplots()
ax.set_xlabel("\$ x \$", fontsize = 22)
ax.set_ylabel("\$ u \$", fontsize = 22)
ax.tick_params(axis="both", which="major", labelsize=23)
ax.tick_params(axis="both", which="minor", labelsize=23)
ax.grid(true, linestyle = "--")
burgers_coarse = readdlm(joinpath(@__DIR__, "..", "run", "burg_stiff", "burg_stiff_shock_coarse_crk_imex433", "avg.txt"))
burgers_coarse_homogeneous = readdlm(joinpath(@__DIR__, "..", "run", "burg_stiff", "burg_stiff_shock_coarse_crk_imex433", "avg_homogeneous.txt"))
burgers_fine = readdlm(joinpath(@__DIR__, "..", "run", "burg_stiff", "burg_stiff_shock_fine_crk_imex433", "avg.txt"))
burgers_fine_homogeneous = readdlm(joinpath(@__DIR__, "..", "run", "burg_stiff", "burg_stiff_shock_fine_crk_imex433", "avg_homogeneous.txt"))
function exact_solution_burg1d(x, t)
    if x - 0.5 * t >= 0.0
        return 0.0
    else
        return 1.0
    end
end
exact(x) = exact_solution_burg1d(x, 4.0)
domain = LinRange(-5.0, 5.0, 1000)
# plot_sol(ax, burgers_sol, 1)
ax.scatter(burgers_coarse_homogeneous[:, 1], burgers_coarse_homogeneous[:, 2], marker = "d", label = "\$ \\Delta x = 1\$, Homg.", facecolors = "none", edgecolors = "orange")
ax.scatter(burgers_coarse[:, 1], burgers_coarse[:, 2], marker = "d", label = "\$ \\Delta x = 1\$", facecolors = "none", edgecolors = "red")
ax.scatter(burgers_fine_homogeneous[:, 1], burgers_fine_homogeneous[:, 2], marker = "o", label = "\$ \\Delta x = 0.1 \$, Homg.", facecolors = "none", edgecolors = "cyan" )
ax.scatter(burgers_fine[:, 1], burgers_fine[:, 2], marker = "o", label = "\$ \\Delta x = 0.1 \$", facecolors = "none", edgecolors = "blue" )
ax.plot(domain, exact.(domain), label = "Exact", c = "k", lw = 2, ls = "-")
ax.legend()
fig
filename = joinpath("burgers_non_linear.pdf")
fig.savefig(filename)
run(`bash pdfbb $filename`)

## Burgers' test (rare faction)
fig, ax = plt.subplots()
ax.set_xlabel("\$ x \$")
ax.set_ylabel("\$ u \$")
ax.set_xlabel("\$ x \$", fontsize = 22)
ax.set_ylabel("\$ u \$", fontsize = 22)
ax.tick_params(axis="both", which="major", labelsize=23)
ax.tick_params(axis="both", which="minor", labelsize=23)
ax.grid(true, linestyle = "--")
burgers_coarse = readdlm(joinpath(@__DIR__, "..", "run", "burg_stiff", "burg_stiff_rare_coarse_crk_imex433", "avg.txt"))
burgers_coarse_homogeneous = readdlm(joinpath(@__DIR__, "..", "run", "burg_stiff", "burg_stiff_rare_coarse_crk_imex433", "avg_homogeneous.txt"))
burgers_fine = readdlm(joinpath(@__DIR__, "..", "run", "burg_stiff", "burg_stiff_rare_fine_crk_imex433", "avg.txt"))
burgers_fine_homogeneous = readdlm(joinpath(@__DIR__, "..", "run", "burg_stiff", "burg_stiff_rare_fine_crk_imex433", "avg_homogeneous.txt"))
function exact_solution_burg1d(x, t)
    if x - 0.9 * t >= 0.0
        return 1.0
    else
        return 0.0
    end
end
# plot_sol(ax, burgers_sol, 1)
ax.scatter(burgers_coarse_homogeneous[:, 1], burgers_coarse_homogeneous[:, 2], marker = "d", label = "\$ \\Delta x = 1\$, Homg.", facecolors = "none", edgecolors = "orange")
ax.scatter(burgers_coarse[:, 1], burgers_coarse[:, 2], marker = "d", label = "\$ \\Delta x = 1\$", facecolors = "none", edgecolors = "red")
ax.scatter(burgers_fine_homogeneous[:, 1], burgers_fine_homogeneous[:, 2], marker = "o", label = "\$ \\Delta x = 0.05 \$, Homg.", facecolors = "none", edgecolors = "cyan" )
ax.scatter(burgers_fine[:, 1], burgers_fine[:, 2], marker = "o", label = "\$ \\Delta x = 0.05 \$", facecolors = "none", edgecolors = "blue" )
exact(x) = exact_solution_burg1d(x, 4.0)
domain = LinRange(-5.0, 5.0, 1000)
ax.plot(domain, exact.(domain), label = "Exact", c = "k", lw = 2, ls = "-")
ax.legend()
fig
filename = joinpath("burgers_rare_non_linear.pdf")
fig.savefig(filename)
run(`bash pdfbb $filename`)

## Super Burgers' test
fig, ax = plt.subplots()
ax.set_xlabel("\$ x \$")
ax.set_ylabel("\$ u \$")
ax.grid(true, linestyle = "--")
ax.set_xlabel("\$ x \$", fontsize = 22)
ax.set_ylabel("\$ u \$", fontsize = 22)
ax.tick_params(axis="both", which="major", labelsize=23)
ax.tick_params(axis="both", which="minor", labelsize=23)
burgers_coarse = readdlm(joinpath(@__DIR__, "..", "run", "sup_burg_stiff", "burg_stiff_shock_1d_coarse_crk_imex433", "avg.txt"))
burgers_coarse_homogeneous = readdlm(joinpath(@__DIR__, "..", "run", "sup_burg_stiff", "burg_stiff_shock_1d_coarse_crk_imex433", "avg_homogeneous.txt"))
burgers_fine = readdlm(joinpath(@__DIR__, "..", "run", "sup_burg_stiff", "burg_stiff_shock_1d_fine_crk_imex433", "avg.txt"))
burgers_fine_homogeneous = readdlm(joinpath(@__DIR__, "..", "run", "sup_burg_stiff", "burg_stiff_shock_1d_fine_crk_imex433", "avg_homogeneous.txt"))
function exact_solution_burg1d(x, t)
    if x - 0.25 * t >= 0.0
        return 0.0
    else
        return 1.0
    end
end
# plot_sol(ax, burgers_sol, 1)
ax.scatter(burgers_coarse_homogeneous[:, 1], burgers_coarse_homogeneous[:, 2], marker = "d", label = "\$ \\Delta x = 2\$, Homg.",  facecolors = "none", edgecolors = "orange")
ax.scatter(burgers_coarse[:, 1], burgers_coarse[:, 2], marker = "d", label = "\$ \\Delta x = 2\$", facecolors = "none", edgecolors = "red")
ax.scatter(burgers_fine_homogeneous[:, 1], burgers_fine_homogeneous[:, 2], marker = "o", label = "\$ \\Delta x = 0.25 \$, Homg.", facecolors = "none", edgecolors = "cyan" )
ax.scatter(burgers_fine[:, 1], burgers_fine[:, 2], marker = "o", label = "\$ \\Delta x = 0.25 \$",
           facecolors="none", edgecolors="green")
exact(x) = exact_solution_burg1d(x, 4.0)
domain = LinRange(burgers_fine[1, 1], burgers_fine[end, 1], 1000)
ax.plot(domain, exact.(domain), label = "Exact", c = "k", lw = 2, ls = "-")
ax.legend()
filename = joinpath("sup_burgers_non_linear.pdf")
fig
fig.savefig(filename)
run(`bash pdfbb $filename`)


# Reactive euler RP solution

# density
fig, ax = plt.subplots()
ax.set_xlabel("\$ x \$", fontsize = 16)
ax.set_ylabel("\$ \\rho \$", fontsize = 16)
ax.tick_params(axis="both", which="major", labelsize=18)
ax.tick_params(axis="both", which="minor", labelsize=18)
ax.grid(true, linestyle = "--")
# blend = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_nx800/avg.txt"))
nx100 = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_nx100/avg.txt"))
nx100_cfl2 = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_nx100_cfl2/avg.txt"))
nx1000 = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_nx1000/avg.txt"))
deg0 = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_deg0/avg.txt"))
exact_rho(x) = x <= 7.1247 ? 1.6812 : 1.0
xdomain = LinRange(-5.0, 25.0, 1000)
ax.plot(xdomain, exact_rho.(xdomain), label = "Exact", lw = 2, c = "k")
# ax.plot(approx[:,1], approx[:,2], label = "Numerical", lw = 2, c = "red")
# ax.plot(xc, ua[1, 1:end-1], label = "Numerical", lw = 2, c = "red")
# plot_sol(ax, tvb, 3, ls = "-", label = "TVB", color = "blue", y_index = 5)
ax.scatter(nx100[:,1], nx100[:,2], marker = "d", label = "\$ \\Delta x = 0.3\$", facecolors = "none", edgecolors = "orange")
ax.scatter(nx100_cfl2[:,1], nx100_cfl2[:,2], marker = "d", label = "\$ \\Delta x = 0.3, C_{\\mathrm{CFL}} = 2\$", facecolors = "none", edgecolors = "red")
ax.scatter(nx1000[:,1], nx1000[:,2], marker = "o", label = "\$ \\Delta x = 0.03\$", facecolors = "none", edgecolors = "blue")
ax.scatter(deg0[:,1], deg0[:,2], marker = "o", label = "\$ \\Delta x = 0.03, \\mathrm{FVM} \$", facecolors = "none", edgecolors = "green")
ax.legend()
ax.set_xlim(-2.0, 9.0)
filename = joinpath("reactive_rp1_density.pdf")
fig
fig.savefig(filename)
run(`bash pdfbb $filename`)

# Pressure
fig, ax = plt.subplots()
ax.set_xlabel("\$ x \$", fontsize = 16)
ax.set_ylabel("\$ p\$", fontsize = 16)
ax.tick_params(axis="both", which="major", labelsize=18)
ax.tick_params(axis="both", which="minor", labelsize=18)
ax.grid(true, linestyle = "--")
# blend = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_nx800/avg.txt"))
nx100 = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_nx100/avg.txt"))
nx100_cfl2 = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_nx100_cfl2/avg.txt"))
nx1000 = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_nx1000/avg.txt"))
deg0 = readdlm(joinpath(@__DIR__, "..", "run", "reactive_rp/reactive_rp1_deg0/avg.txt"))
exact_pressure(x) = x <= 7.1247 ? 21.562 : 1.0
xdomain = LinRange(-5.0, 25.0, 1000)
ax.plot(xdomain, exact_pressure.(xdomain), label = "Exact", lw = 2, c = "k")
ax.scatter(nx100[:,1], nx100[:,4], marker = "d", label = "\$ \\Delta x = 0.3\$", facecolors = "none", edgecolor = "orange")
ax.scatter(nx100_cfl2[:,1], nx100_cfl2[:,4], marker = "d", label = "\$ \\Delta x = 0.3, C_{\\mathrm{CFL}} = 2\$", facecolors = "none", edgecolor = "red")
ax.scatter(nx1000[:,1], nx1000[:,4], marker = "o", label = "\$ \\Delta x = 0.03\$", facecolors = "none", edgecolor = "blue")
ax.scatter(deg0[:,1], deg0[:,4], marker = "o", label = "\$ \\Delta x = 0.03, \\mathrm{FVM} \$", facecolors = "none",
           edgecolor = "green")
ax.legend()
ax.set_xlim(-2.0, 9.0)
filename = joinpath("reactive_rp1_pressure.pdf")
fig
fig.savefig(filename)
run(`bash pdfbb $filename`)

# Variable advection
exact(x) = sinpi(x^3 - 3.0 * 1.0) # final_time = 1.0
fig, ax = plt.subplots()
ax.set_xlabel("\$ x \$", fontsize = 16)
ax.set_ylabel("\$ u \$", fontsize = 16)
ax.tick_params(axis="both", which="major", labelsize=18)
ax.tick_params(axis="both", which="minor", labelsize=18)
ax.grid(true, linestyle = "--")
vardv_sol = readdlm(joinpath(@__DIR__, "..", "run", "varadv", "output_varadv", "sol.txt"))
x = LinRange(0.1, 1.0, 1000)
ax.plot(x, exact.(x), label = "Exact", c = "k", lw = 2, ls = "-")
plot_sol(ax, vardv_sol, 3)
ax.legend()
filename = joinpath("var_adv.pdf")
fig.savefig(filename)
fig
run(`bash pdfbb $filename`)
