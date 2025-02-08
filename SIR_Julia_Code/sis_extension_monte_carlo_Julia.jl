using Random, Images, FileIO, Plots

# Parameters
grid = 500 # Matrix (space) size
times = 5000  # Number of iterations
initial = 500  # Initial number of infected individuals
StoE = 80

EtoI = 0.05
ItoR = 0.04
ItoQ = 0.02
QtoR = 0.04
RtoS = 0.003
# Initialize matrices
M = zeros(Int, grid, grid)

# Randomly generate initial infected locations
for _ in 1:initial
    M[rand(1:grid), rand(1:grid)] = 2
end
N=copy(M)
heatmap(M)
heatmap(N)
# Simulation loop
anim = Animation()
S, E, I, R, Q, n = 0, 0, 0, 0, 0, 0
for t in 1:times      
    for i in 1:grid
    for j in 1:grid
        if M[i, j] == -1 # R to S transition
            if rand() < RtoS
                N[i, j] = 0
            end
        end
        if M[i, j] == 1 # E to I transition
            if rand() < EtoI
                N[i, j] = 2
            end
            if rand(1:StoE) == 1
                neighbors = [(1, -1), (1, 0), (1, 1), (0, -1), (0, 1), (-1, -1), (-1, 0), (-1, 1)]
                dx, dy = rand(neighbors)
                ni, nj = i + dx, j + dy
                if 1 ≤ ni ≤ grid && 1 ≤ nj ≤ grid && M[ni, nj] in (0, 1)
                    N[ni, nj] = 1
                    n += 1
                end
            end
        end
        if M[i, j] == 2 # I to R or Q transition
            x = rand()
            if x < ItoR
                N[i, j] = -1
            end
            if ItoR<=x && x<ItoQ+ItoR
                N[i,j]=-2;
            end
            if x>=ItoQ+ItoR
                rand(1:StoE) == 1 
                neighbors = [(1, -1), (1, 0), (1, 1), (0, -1), (0, 1), (-1, -1), (-1, 0), (-1, 1)]
                dx, dy = rand(neighbors)
                ni, nj = i + dx, j + dy
                if 1 ≤ ni ≤ grid && 1 ≤ nj ≤ grid && M[ni, nj] in (0, 1)
                    N[ni, nj] = 1
                    n += 1
                end
            end
        end
        if M[i, j] == -2 # Q to R transition
            if rand() < QtoR
                N[i, j] = -1
            end
        end
    end
    end

    if N == M
        break
    end
    M = copy(N)

    # Count population categories
    for i in 1:grid, j in 1:grid
        if M[i, j] == -1
            R += 1
        elseif M[i, j] == 0
            S += 1
        elseif M[i, j] == 1
            E += 1
        elseif M[i, j] == 2
            I += 1
        elseif M[i, j] == -2
            Q += 1
        end
    end
    # Save animation frames every 100 steps
    if t % 10 == 0 || t == 2
        heatmap(M)
        frame(anim)
    end
end
heatmap(N)
heatmap(M)
gif(anim, "sis.gif", fps=10)