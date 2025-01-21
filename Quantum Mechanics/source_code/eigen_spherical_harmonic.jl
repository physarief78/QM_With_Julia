using Plots
plotlyjs()  # Set PlotlyJS as the backend

# Compute the factorial of a number (custom implementation)
function my_factorial(n::Int)
    n <= 1 ? 1 : prod(1:n)
end

# Compute the associated Legendre polynomial P_l^m(x) manually
function assoc_legendre(l, m, x)
    if m < 0
        throw(ArgumentError("Order m must be non-negative"))
    end

    # Initial values for recurrence relation
    p_prev = 1.0  # P_0^0(x)
    p_curr = x    # P_1^0(x) for m = 0

    # Apply recurrence for m = 0 to get P_l^0(x)
    for n in 1:l-1
        p_next = ((2n + 1) * x * p_curr - n * p_prev) / (n + 1)
        p_prev, p_curr = p_curr, p_next
    end

    p_lm = p_curr  # P_l^0(x)

    # Adjust for m > 0 using recurrence relation
    for k in 1:m
        p_lm = -sqrt(1 - x^2) * (2k + 1) * p_lm
    end
    return p_lm
end

# Define the spherical harmonics Y_l^m(θ, φ)
function spherical_harmonic(l, m, θ, φ)
    normalization = sqrt((2l + 1) / (4 * π) * my_factorial(l - abs(m)) / my_factorial(l + abs(m)))
    legendre = assoc_legendre(l, abs(m), cos(θ))
    phase = exp(im * m * φ)
    return normalization * legendre * phase
end

# Generate the 3D perspective plot of spherical harmonics
function plot_spherical_harmonics(l, m, filename="spherical_harmonic.png")
    # Create a grid of θ (colatitude) and φ (longitude)
    θ = range(0, π, length=100)
    φ = range(0, 2π, length=100)

    # Compute |Y_l^m|^2 on the spherical surface
    r = [abs(spherical_harmonic(l, m, θi, φj))^2 for θi in θ, φj in φ]

    # Convert spherical coordinates to Cartesian coordinates for 3D plotting
    x = [r[i, j] * sin(θ[i]) * cos(φ[j]) for i in 1:length(θ), j in 1:length(φ)]
    y = [r[i, j] * sin(θ[i]) * sin(φ[j]) for i in 1:length(θ), j in 1:length(φ)]
    z = [r[i, j] * cos(θ[i]) for i in 1:length(θ), j in 1:length(φ)]

    # Create the 3D plot
    plt = surface(x, y, z, color=:viridis, title="Spherical Harmonic Y_$l^$m", xlabel="x", ylabel="y", zlabel="z")

    # Save the plot
    #savefig(plt, filename)
    #println("3D Spherical Harmonic plot saved as $filename")

    # Display the interactive plot
    display(plt)
end

# Plot the 3D spherical harmonic for given l and m
plot_spherical_harmonics(5, 3, "spherical_harmonic_3_2.png")
