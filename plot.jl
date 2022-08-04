using Pkg

Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("LaTeXStrings")
Pkg.add("Plots")

using CSV
using DataFrames
using LaTeXStrings
using Plots

#CHOOSE CORRECT LIBRARY!
#The numerical library used to compute values (either eigen or armadillo)
#########################################################################
#########################################################################
#########################################################################
num_lib = "eigen"
#########################################################################
#########################################################################
#########################################################################

#Colour theme of plots
theme(:solarized)

# #Time steps when computing the hamiltonians
hamilton_h = [0.1, 0.01]
hamiltonians_files = ["./$(num_lib)/output/hamiltonians_$(h).csv" for h in hamilton_h]

#File names of poincare map values
poincare_rk(h) = "./$(num_lib)/output/poincare_rk4_$(h).csv"
poincare_sb(h) = "./$(num_lib)/output/poincare_sb_$(h).csv"
poincare_kahans(h) = "./$(num_lib)/output/poincare_kahans_$(h).csv"
poincare_sv(h) = "./$(num_lib)/output/poincare_sv_$(h).csv"

# Initial energy of system
H_0 = 1.0/12.0 

function plot_commands(x_label, y_label, title = "")
    #Function containing the usual commands when plotting with Plots
    xlabel!(x_label)
    ylabel!(y_label)
    if !isempty(title)
        title!(title)
    end
end

function read_hamiltonians(ham_file)

    #Read all of the hamiltonians
    df = Matrix(CSV.read(ham_file, DataFrame; header = 0))
    t, H = df[:,1], df[:,2:end]

    return t, H
end

function plot_all_hamiltonians(h, ham_file, fig_name)
    # Size of plot
    fig_size = (800,600)

    #Read the hamiltonians from file
    t, H = read_hamiltonians(ham_file)

    #Plot the lines
    plot(
        t,
        H,
        label = ["Kutta's method" "Shampine-Bogacki method" "Kahan's method" "Störmer-Verlet method"],
        size = fig_size,
        legend = :right
    )
    
    #Attach labels and save figure
    plot_commands("Time [s] (dt = $(h))", "H(p,q)", "Plot of the H(p,q) for various numerical methods")
    savefig(fig_name)
end

function plot_all_hamiltonians_individual(h, ham_file, fig_name)
    #Plot every method as subplots in the same figure
    #The plots show the deviation of a particular numerical method,
    #from the constant energy the system actually has

    # Size of plot
    fig_size = (800,600)

    #Read the hamiltonians from file
    t, H = read_hamiltonians(ham_file)

    const_energy = fill(H_0, size(t))

    #Kutta's method
    p1 = plot(
        t,
        [H[:,1], const_energy],
        title = "Kutta's method",
    )

    #Shampine-Bogacki
    p2 = plot(
        t,
        [H[:,2], const_energy],
        title = "Shampine-Bogacki",
    )

    #Kahan's method
    p3 = plot(
        t,
        [H[:,3], const_energy],
        title = "Kahan's method",
    )

    #Störmer-Verlet
    p4 = plot(
        t,
        [H[:,4], const_energy],
        title = "Störmer-Verlet",
    )

    p = plot(
        p1, p2, p3, p4,
        layout = (2,2),
        title = ["Kutta's method" "Shampine-Bogacki method" "Kahan's method" "Störmer-Verlet method"],
        plot_title = "Error of Hamiltonian for Hénon-Heiles system for different numerical methods",
        label = ["Numerical" "Exact"],
        size = fig_size,
        lw = 2
    )

    #Attach labels and save figure
    plot_commands("Time [s] (dt = $(h))", "H(p,q)")
    savefig(fig_name)
end

function read_poincare(h)
    #Read the poincare files and store them in arrays
    P_rk = Matrix(CSV.read(poincare_rk(h), DataFrame; header = 0))
    P_sb = Matrix(CSV.read(poincare_sb(h), DataFrame; header = 0))
    P_kahans = Matrix(CSV.read(poincare_kahans(h), DataFrame; header = 0))
    P_sv = Matrix(CSV.read(poincare_sv(h), DataFrame; header = 0))

    return P_rk, P_sb, P_kahans, P_sv
end

function poincare_all_methods(h_arr, labels, fig_names)
    #Plot the Poincaré map for all methods for the different time steps used

    # Size of plot
    fig_size = (800,600)

    for (h, label, fig_name) in zip(h_arr, labels, fig_names)
        P_rk, P_sb, P_kahans, P_sv = read_poincare(h)

        #Kutta's method
        p1 = plot(
            P_rk[1,:],
            P_rk[2,:],
            legend = false,
            seriestype = :scatter
        )

        #Shampine-Bogacki
        p2 = plot(
            P_sb[1,:],
            P_sb[2,:],
            legend = false,
            seriestype = :scatter
        )

        #Kahan's method
        p3 = plot(
            P_kahans[1,:],
            P_kahans[2,:],
            legend = false,
            seriestype = :scatter
        )

        #Störmer-Verlet
        p4 = plot(
            P_sv[1,:],
            P_sv[2,:],
            legend = false,
            seriestype = :scatter
        )

        plot(
            p1, p2, p3, p4,
            layout = (2,2),
            title = ["Kutta's method" "Shampine-Bogacki method" "Kahan's method" "Störmer-Verlet method"],
            plot_title = label,
            size = fig_size
        )
        
        plot_commands(L"q_2", L"p_2")
        savefig(fig_name)
    end
end

#Folder for storing plots
output_folder = "./$(num_lib)/plots"
mkpath(output_folder)
for (h, ham_file) in zip(hamilton_h, hamiltonians_files)
    plot_all_hamiltonians(h, ham_file, "$(output_folder)/hamiltonians_$(h).png")
    plot_all_hamiltonians_individual(h, ham_file, "$(output_folder)/hamiltonians_individual_$(h).png")
end

#All the time steps for poincare maps
h_arr = [0.0001, 0.001, 0.01, 0.1]
t_arr = [3e3, 3e4, 3e5, 3e6]

fig_names = ["$(output_folder)/poincare_$(h).png" for h in h_arr]

labels = []
for (t, h) in zip(t_arr, h_arr)
    append!(labels, "T = $(t)s & h = $(h)s")
end

poincare_all_methods(h_arr, labels, fig_names)