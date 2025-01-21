using Plots

# Constants for the 1D infinite square well
L = 1.0  # Width of the well
ħ = 1.0  # Reduced Planck's constant (set to 1 in natural units)
m = 1.0  # Mass of the particle (set to 1 in natural units)

# Define the eigenfunction for the n-th quantum state
function eigenfunction(x, n)
    return sqrt(2 / L) * sin(n * π * x / L)  # Normalized eigenfunction
end

# Compute the Fourier series by summing eigenfunctions from n = 1 to max_n
function fourier_series(x, max_n)
    ψ = 0.0
    for n in 1:max_n
        ψ += eigenfunction(x, n)
    end
    return ψ
end

# Plot the Fourier series of summed eigenfunctions and the individual eigenfunctions
function plot_fourier_series_and_eigenfunctions(max_n, filename="plot.png")
    # Define the x range
    x = range(0, L, length=1000)

    # Compute the Fourier series
    ψ_fourier = [fourier_series(xi, max_n) for xi in x]

    # Create the plot
    plt = plot(layout=(2, 1), size=(1200, 1000))

    # Plot the individual eigenfunctions
    plot!(plt[1], title="Eigenfunctions", xlabel="x", ylabel="ψ(x)")
    for n in 1:max_n
        ψ_eigen = [eigenfunction(xi, n) for xi in x]
        plot!(plt[1], x, ψ_eigen, label="n=$n")
    end

    # Plot the Fourier series
    plot!(plt[2], x, ψ_fourier, label="Fourier Series (sum of n=1 to $max_n)",
          title="Fourier Series of Summed Eigenfunctions", xlabel="x", ylabel="ψ(x)")

    # Save the plot
    savefig(plt, filename)
    println("Plot saved as $filename")
end

# Simulate the Fourier series of summed eigenfunctions and plot alongside eigenfunctions
plot_fourier_series_and_eigenfunctions(10, "fourier_series_plot.png")
