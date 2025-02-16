using Random, DelimitedFiles
using Plots

# Global variables
global MagMean 
global EnMean 
global Mag 
global En 
global En2Mean 
global Mag2Mean
global Check

function initialize_lattice(L, IC)
    if IC == 1
        spin = ones(Int, L, L)
    elseif IC == 2
        spin = -ones(Int, L, L)
    else
        spin = sign.(0.5 .- rand(L, L))
    end
    
    En = 0
    for x in 1:L, y in 1:L
        En -= spin[x, y] * (spin[mod1(x+1, L), y] + spin[x, mod1(y+1, L)])
    end
    
    Mag = sum(spin)
    return spin, En, Mag
end

function monte_carlo_step!(spin, L, T)
    global Mag, En 
    row, col = rand(1:L), rand(1:L)
    dE = 2 * spin[row, col] * (
        spin[mod1(row-1, L), col] + spin[mod1(row+1, L), col] + 
        spin[row, mod1(col-1, L)] + spin[row, mod1(col+1, L)]
    )
    if dE >= 0 
        if rand() < exp(-dE / T)
            spin[row, col] *= -1
        end
    else
        spin[row, col] *= -1
    end
    for x in 1:L, y in 1:L
        En -= spin[x, y] * (spin[mod1(x+1, L), y] + spin[x, mod1(y+1, L)])
    end
    
    Mag = sum(spin)
    return En, Mag
end

# Parameters
L = 30 # Lattice size
IC = 3 # 1 for all spin up, 2 for all spin down, 3 for random
T = 2.5 # Temperature
steps = 10^6 # Number of simulation steps

# Initialize lattice
spin, En, Mag = initialize_lattice(L, IC)

# Run simulation


i = 1
Temp = 1.5:0.1:3
Energy = Float64[]
Magnetization = Float64[]
MagSus = Float64[]
SpesHeat = Float64[]

for T in Temp
    for t in 1:steps
        En, Mag = monte_carlo_step!(spin, L, T)
        global EnMean = En / L^2
        global MagMean = Mag / L^2
    end
    push!(Energy, EnMean)
    push!(Magnetization, MagMean)
    push!(MagSus, ((L^2)/T) * (Mag2Mean - MagMean^2))
    push!(SpesHeat, ((L/T)^2) * (En2Mean - EnMean^2))
    i += 1
end

plot(Temp, Energy, label="Mean Energy", seriestype=:scatter)
plot(Temp, Magnetization, label="Mean Magnetization", seriestype=:scatter)
plot(Temp, SpesHeat, label="Specific Heat", seriestype=:scatter)
plot(Temp, MagSus, label="Magnetic Susceptibility", seriestype=:scatter)
