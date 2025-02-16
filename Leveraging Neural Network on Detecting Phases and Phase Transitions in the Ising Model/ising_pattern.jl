using Random, Plots

function initialize_lattice(L, IC)
    if IC == 1
        spin = ones(Int, L, L)
    elseif IC == 2
        spin = -ones(Int, L, L)
    else
        spin = sign.(0.5 .- rand(L, L))
    end
    return spin
end

function monte_carlo_step!(spin, L, T)
    row, col = rand(1:L), rand(1:L)
    dE = 2 * spin[row, col] * (
        spin[mod1(row-1, L), col] + spin[mod1(row+1, L), col] + 
        spin[row, mod1(col-1, L)] + spin[row, mod1(col+1, L)]
    )
    
    if dE >= 0 
        if rand() < exp(-dE / T)
            spin[row, col] *= -1
            return 0  
        end
    else
        spin[row, col] *= -1
            return 0 
    end

end

# Parameters
L = 200 # Lattice size
IC = 3 # 1 for all spin up, 2 for all spin down, 3 for random
T = 1 # Temperature
steps = 10^5 # Number of simulation steps

# Initialize lattice
spin = initialize_lattice(L, IC)

# Run simulation with animation
anim = Animation()
for t in 1:steps
    monte_carlo_step!(spin, L, T)
    if t % 1000 == 1 
    heatmap(spin, color = :thermal, title = "lattice size=$(L)×$(L), temperature=$(T), time=$(t)")
    frame(anim)
    end
end

gif(anim, "ising.gif", fps = 50)


