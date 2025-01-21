using Plots

# Constants for the 2D infinite square well
Lx = 1.0  # Width of the well in x-direction
Ly = 1.0  # Width of the well in y-direction
ħ = 1.0  # Reduced Planck's constant (set to 1 in natural units)
m = 1.0  # Mass of the particle (set to 1 in natural units)

# Define the eigenfunction for the (nx, ny)-th quantum state
function eigenfunction_2D(x, y, nx, ny)
    return sqrt(4 / (Lx * Ly)) * sin(nx * π * x / Lx) * sin(ny * π * y / Ly)  # Normalized eigenfunction
end

# Compute the Fourier series by summing eigenfunctions for (nx, ny) combinations
function fourier_series_2D(x, y, max_nx, max_ny)
    ψ = 0.0
    for nx in 1:max_nx
        for ny in 1:max_ny
            ψ += eigenfunction_2D(x, y, nx, ny)
        end
    end
    return ψ
end

# Plot the Fourier series of summed eigenfunctions and the individual eigenfunctions
function plot_fourier_series_and_eigenfunctions_2D(max_nx, max_ny, filename="plot_2D.png")
    # Define the x and y ranges
    x = range(0, Lx, length=100)
    y = range(0, Ly, length=100)

    # Compute the Fourier series
    ψ_fourier = [fourier_series_2D(xi, yi, max_nx, max_ny) for xi in x, yi in y]

    # Create the plot for the Fourier series
    plt_fourier = heatmap(x, y, ψ_fourier, title="Fourier Series of Summed Eigenfunctions",
                          xlabel="x", ylabel="y", colorbar_title="ψ(x, y)")

    # Create plots for individual eigenfunctions
    plt_eigenfunctions = plot(layout=(max_nx, max_ny), size=(1200, 1000))
    for nx in 1:max_nx
        for ny in 1:max_ny
            ψ_eigen = [eigenfunction_2D(xi, yi, nx, ny) for xi in x, yi in y]
            heatmap!(plt_eigenfunctions[(nx - 1) * max_ny + ny], x, y, ψ_eigen,
                     title="Eigenfunction (nx=$nx, ny=$ny)", xlabel="x", ylabel="y",
                     colorbar=false)
        end
    end

    # Save the plots
    savefig(plt_fourier, "fourier_series_" * filename)
    savefig(plt_eigenfunctions, "eigenfunctions_" * filename)
    println("Plots saved as fourier_series_$filename and eigenfunctions_$filename")
end

# Simulate the 2D Fourier series of summed eigenfunctions and plot alongside individual eigenfunctions
plot_fourier_series_and_eigenfunctions_2D(3, 3, "2D_plot.png")
