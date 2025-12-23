##? Exercise 2 - MuLoe
##* load packages 
using Revise        # for development purposes
using LinearAlgebra # for solving linear systems
using Random        # for random number generation
using Plots         # for plotting
plotly()            # use the plotly backend

#* load functions in the files
include("fun.jl")

##* Burn of natural gas
# Elemental matrix E
# 1 - CH4, 2 - C2H6, 3 - C3H8, 4 - O2, 5 - N2, 6 - CO2, 7 - H2O
E = [1 2 3 0 0 1 0 ;  # C
     4 6 8 0 0 0 2 ;  # H
     0 0 0 2 0 2 1 ;  # O
     0 0 0 0 2 0 0 ]; # N

# composition of the natural gas
x_f = [0.93 0.03 0.02 0.01 0.01 0]  

# composition of the air
x_air_in = [0.21; 0.79; 0; 0]          # mol/mol   

# investigated λ range
λ = 1.0:0.1:2.0

# calculate the stoechiometric matrix
nu,index_iskey,isfirst = stoechiometrie(E)

# preallocate variables with **7 columns** to match exhaust gas size
X_air_min = zeros(length(λ))
X_air_in  = zeros(length(λ))
X_O2_min  = zeros(length(λ))
X_O2_in   = zeros(length(λ))
X_in      = zeros(length(λ), 7)
X_ex_out  = zeros(length(λ), 7)
x_ex_out  = zeros(length(λ), 7)

# calculate the air demand for the given lambda values
for i = 1:length(λ)
    X_air_min[i], X_air_in[i], X_O2_min[i], X_O2_in[i] = air_demand_molar(x_f, λ[i], x_air_in, nu)
    X_in[i, 1:length(x_air_in)] = X_air_in[i] .* x_air_in   # only fill the first 4
    X_ex_out[i, :] = exhaust_gas(x_f, X_in[i, :], nu)
    x_ex_out[i, :] = X_ex_out[i, :] ./ sum(X_ex_out[i, :], dims=1)
end

# plot the results
plot(λ, X_air_in, xlabel="λ", ylabel="Air demand", label="X_air_in")
plot!(λ, X_air_min, label="X_air_min")
plot(λ, X_ex_out[:, 1:4], xlabel="λ", ylabel="X_ex_out", label=["O2" "N2" "CO2" "H2O"])
plot(λ, x_ex_out[:, 1:4], xlabel="λ", ylabel="x_ex_out", label=["O2" "N2" "CO2" "H2O"])

##* Burn of black coal
# 1 = C, 2 = H, 3 = S, 4 = O2, 5 = N2, 6 = H2O, 7 = solid inerts
w_f = [0.72 0.05 0.01 0.1 0.01 0.06 0.05] # kg/kg

# molar mass of coal components 1=C,2=H,3=S,4=O2,5=N2,6=H2O
MW_f = [12.011 1.008 32.065 32 28 18]./1000  

# composition of the air
x_air_in = [0.21 0.79 0 0];         # mol/mol   
MW_air = [32 28 44 18]./1000         # kg/mol

# Stoichiometric coefficients of exhaust gas
# 1 = O2, 2 = N2, 3 = CO2, 4 = H2O, 5 = SO2
nu       =[-1    0  1 0   0;  # C + O2 -> CO2
           -0.25 0  0 0.5 0;  # H + 1/4 O2 -> 1/2 H2O
           -1    0  0 0   1]   # S + O2 -> SO2

# molar mass exhaust gas components
MW_ex = [15.999*2 2*14.007 44.01 18.0153 64.066]./1000 # kg/mol

# convert to mass fractions
w_air_in, _ = mol2mass(reshape(x_air_in,1,4), MW_air)

# investigated λ range
λ = 1.0:0.1:2.0

# preallocate variables
W_air_min = zeros(length(λ))
W_air_in  = zeros(length(λ))
W_O2_min  = zeros(length(λ))
W_O2_in   = zeros(length(λ))
W_in      = zeros(length(λ), 5)
W_ex_out  = zeros(length(λ), 5)
w_ex_out  = zeros(length(λ), 5)
x_ex_out_mass = zeros(length(λ), 5)

# calculate the air demand for the given lambda values
for i = 1:length(λ)
    W_air_min[i], W_air_in[i], W_O2_min[i], W_O2_in[i] = air_demand_mass(w_f, λ[i], w_air_in, MW_f, nu)
    W_in[i, 1:4] = W_air_in[i] .* w_air_in
    W_ex_out[i, :] = exhaust_gas_mass(w_f, MW_f, W_in[i, :], MW_ex, nu)
    w_ex_out[i, :] = W_ex_out[i, :] ./ sum(W_ex_out[i, :], dims=1)
    x_ex_out_mass[i, :], _ = mass2mol(reshape(w_ex_out[i, :], 1, 5), MW_ex)
end

# plot the results
plot(λ, W_air_in, xlabel="λ", ylabel="Air demand", label="W_air_in")
plot!(λ, W_air_min, label="W_air_min")
plot(λ, W_ex_out, xlabel="λ", ylabel="X_ex_out", label=["O2" "N2" "CO2" "H2O" "SO2"])
plot(λ, w_ex_out, xlabel="λ", ylabel="w_ex_out / kg/kg", label=["O2" "N2" "CO2" "H2O" "SO2"])
plot(λ, x_ex_out_mass, xlabel="λ", ylabel="x_ex_out / mol/mol", label=["O2" "N2" "CO2" "H2O" "SO2"])
