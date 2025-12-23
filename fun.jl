# fun.jl
using LinearAlgebra

"""
mol2mass converts mole fractions to mass fractions
x: mole fractions size KxN
MW: molecular weights size 1xN
"""
function mol2mass(x::Matrix, MW::Matrix)
    if size(x, 2) != size(MW, 2)
        error("The number of species in x and MW must be the same")
    end
    MWges = sum(x .* MW, dims=2)
    xm = x .* MW ./ MWges
    return xm, MWges
end

"""
mass2mol converts mass fractions to mole fractions
xm: mass fractions size KxN
MW: molecular weights size 1xN
"""
function mass2mol(xm::Matrix, MW::Matrix)
    if size(xm, 2) != size(MW, 2)
        error("The number of species in xm and MW must be the same")
    end
    x = xm ./ MW ./ sum(xm ./ MW, dims=2)
    MWges = sum(MW .* x, dims=2)
    return x, MWges
end

"""
stoichiometric matrix calculation
E: element matrix size LxN
"""
function stoechiometrie(E)
    L, N = size(E)
    rankE = rank(E)
    if rankE < L
        error("Rank of E smaller than number of elements")
    end
    R = N - rankE
    N_notkey = N - R
    N_key = N - N_notkey

    index_notkey = N-N_notkey+1:N
    index_iskey = setdiff(1:N, index_notkey)
    E_stern = E[:, index_notkey]

    nu_notkey = zeros(R, N_notkey)
    nu = zeros(R, N)

    if rank(E_stern) >= rankE
        isfirst = true
        nu_key = -Diagonal(ones(R))
        for i = 1:R
            nu_notkey[i,:] = E_stern \ (-E[:, index_iskey[i]] * nu_key[i,i])
        end
        nu[:, index_notkey] = nu_notkey
        nu[:, index_iskey] = nu_key
    else
        @warn("First R species are not key species, auto guess")
        isfirst = false
        fail = true
        counter = 0
        while fail
            try
                index_notkey = randperm(N)[1:N_notkey]
                index_iskey = setdiff(1:N, index_notkey)
                E_stern = E[:, index_notkey]
                nu_key = -Diagonal(ones(R))
                for i = 1:R
                    nu_notkey[i,:] = E_stern \ (-E[:, index_iskey[i]] * nu_key[i,i])
                end
                nu[:, index_notkey] = nu_notkey
                nu[:, index_iskey] = nu_key
                if norm(E*nu') == 0
                    fail = false
                else
                    fail = true
                end
            catch
                fail = true
            end
            counter += 1
            println("counter: ", counter)
        end
    end
    return nu, index_iskey, isfirst
end

"""
Air demand (molar mode)
"""
function air_demand_molar(x_f, λ, x_air_in, nu)
    Nf = length(x_f)
    R, N = size(nu)
    x_f_ext = zeros(N)
    x_f_ext[1:Nf] = x_f
    O2_idx = 4
    O2_air_idx = 1
    O2_min = -sum(nu[:,O2_idx] .* x_f_ext[1:R])
    X_air_min = O2_min / x_air_in[O2_air_idx]
    X_air_in = λ * X_air_min
    X_O2_min = O2_min
    X_O2_in = X_air_in * x_air_in[O2_air_idx]
    return X_air_min, X_air_in, X_O2_min, X_O2_in
end

"""
Air demand (mass mode)
"""
function air_demand_mass(w_f, λ, w_air_in, MW_f, nu)
    Nf = length(w_f)
    w_f_ext = zeros(7)
    w_f_ext[1:Nf] = w_f
    w_air_ext = zeros(5)
    w_air_ext[1:length(w_air_in)] = w_air_in
    MW_C, MW_H, MW_S, MW_O2 = MW_f[1:4]
    nC = w_f_ext[1]/MW_C
    nH = w_f_ext[2]/MW_H
    nS = w_f_ext[3]/MW_S
    nO2_fuel = w_f_ext[4]/MW_O2
    O2_min = nC*1 + nH/4 + nS*1 - nO2_fuel
    O2_in = λ * O2_min
    W_air_min = O2_min * MW_O2 / w_air_ext[1]
    W_air_in = λ * W_air_min
    W_O2_min = O2_min * MW_O2
    W_O2_in = O2_in * MW_O2
    return W_air_min, W_air_in, W_O2_min, W_O2_in
end

"""
Exhaust gas (molar mode)
"""
function exhaust_gas(x_f, X_in, nu)
    Nf = length(x_f)
    R, N = size(nu)
    x_f_ext = zeros(N)
    x_f_ext[1:Nf] = x_f
    X_in_ext = zeros(N)
    X_in_ext[1:length(X_in)] = X_in
    X_ex_out = x_f_ext + X_in_ext + nu' * x_f_ext[1:R]
    return X_ex_out
end

"""
Exhaust gas (mass mode)
"""
function exhaust_gas_mass(w_f, MW_f, W_in, MW_ex, nu)
    Nf = length(w_f)
    Na = length(W_in)
    w_f_ext = zeros(7)
    w_f_ext[1:Nf] = w_f
    W_in_ext = zeros(5)
    W_in_ext[1:Na] = W_in
    n_O2 = W_in_ext[1]/MW_ex[1] + w_f_ext[4]/MW_f[4]
    n_N2 = W_in_ext[2]/MW_ex[2] + w_f_ext[5]/MW_f[5]
    n_CO2 = w_f_ext[1]/MW_f[1]
    n_H2O = w_f_ext[2]/MW_f[2]/2 + W_in_ext[4]/MW_ex[4]
    n_SO2 = w_f_ext[3]/MW_f[3]
    mass = [n_O2*MW_ex[1], n_N2*MW_ex[2], n_CO2*MW_ex[3], n_H2O*MW_ex[4], n_SO2*MW_ex[5]]
    W_ex_out = mass ./ sum(mass)
    return W_ex_out
end



