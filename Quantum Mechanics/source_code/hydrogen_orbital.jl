using LinearAlgebra
using PlotlyJS

# Constants
ϵ₀ = 8.854187817e-12  # Vacuum permittivity (F/m)
e = 1.602176634e-19   # Elementary charge (C)
ħ = 1.054571817e-34   # Reduced Planck's constant (J·s)
mₑ = 9.10938356e-31   # Electron mass (kg)
a₀ = 4 * π * ϵ₀ * ħ^2 / (mₑ * e^2)  # Bohr radius (m)

# Factorial function (manual implementation)
function factorial_manual(n::Int)
    return n == 0 ? 1 : n * factorial_manual(n - 1)
end

# Generalized Laguerre polynomial
function genlaguerre(n, α, x)
    if n == 0
        return 1
    elseif n == 1
        return 1 + α - x
    else
        Lₙ₋₁ = genlaguerre(n - 1, α, x)
        Lₙ₋₂ = genlaguerre(n - 2, α, x)
        return ((2n - 1 + α - x) * Lₙ₋₁ - (n - 1 + α) * Lₙ₋₂) / n
    end
end

# Associated Legendre polynomial
function legendre_p(ℓ, m, x)
    if ℓ == 0
        return 1
    elseif ℓ == 1
        return x
    else
        Pℓ₋₁ = legendre_p(ℓ - 1, m, x)
        Pℓ₋₂ = legendre_p(ℓ - 2, m, x)
        return ((2ℓ - 1) * x * Pℓ₋₁ - (ℓ + m - 1) * Pℓ₋₂) / ℓ
    end
end

# Radial wavefunction
function radial_wavefunction(n, ℓ, r)
    ρ = 2 * r / (n * a₀)
    N = sqrt((2.0 / (n * a₀))^3 * factorial_manual(n - ℓ - 1) / (2 * n * factorial_manual(n + ℓ)))
    L = genlaguerre(n - ℓ - 1, 2ℓ + 1, ρ)
    return N * ρ^ℓ * exp(-ρ / 2) * L
end

# Spherical harmonics
function spherical_harmonic(ℓ, m, θ, φ)
    Pₗₘ = legendre_p(ℓ, abs(m), cos(θ))
    normalization = sqrt((2ℓ + 1) / (4π) * factorial_manual(ℓ - abs(m)) / factorial_manual(ℓ + abs(m)))
    return normalization * Pₗₘ * exp(1im * m * φ)
end

# Hydrogen wavefunction
function hydrogen_wavefunction(n, ℓ, m, r, θ, φ)
    R = radial_wavefunction(n, ℓ, r)
    Y = spherical_harmonic(ℓ, m, θ, φ)
    return R * Y
end

function generate_density_surface_data(n, ℓ, m, resolution=100)
    θ_vals = range(0, π, length=resolution)  # Polar angle
    φ_vals = range(0, 2π, length=resolution)  # Azimuthal angle

    x = zeros(Float64, resolution, resolution)
    y = zeros(Float64, resolution, resolution)
    z = zeros(Float64, resolution, resolution)
    prob_density = zeros(Float64, resolution, resolution)

    max_density = 0.0  # Track the maximum density for normalization

    for i in 1:resolution, j in 1:resolution
        θ = θ_vals[i]
        φ = φ_vals[j]

        # Integrate radial component with angular components
        r = 10 * a₀  # Arbitrary radial range for sampling
        ψ = hydrogen_wavefunction(n, ℓ, m, r, θ, φ)
        P = abs(ψ)^2  # Probability density

        max_density = max(max_density, P)

        # Scale radius based on probability density
        scaled_r = sqrt(P) * r  # Probability density contributes to shape

        # Convert spherical to Cartesian coordinates
        x[i, j] = scaled_r * sin(θ) * cos(φ)
        y[i, j] = scaled_r * sin(θ) * sin(φ)
        z[i, j] = scaled_r * cos(θ)
        prob_density[i, j] = P
    end

    # Normalize probability density for consistent color scaling
    prob_density ./= max_density
    return x, y, z, prob_density
end

function plot_density_surface_wavefunction(n, ℓ, m)
    x, y, z, prob_density = generate_density_surface_data(n, ℓ, m)

    surface_trace = PlotlyJS.surface(
        x=x,
        y=y,
        z=z,
        surfacecolor=prob_density,
        colorscale="Viridis",
        showscale=true
    )

    layout = PlotlyJS.Layout(
        title="Hydrogen Atom Orbital Density Shape (n=$n, ℓ=$ℓ, m=$m)",
        scene=PlotlyJS.attr(
            xaxis_title="x",
            yaxis_title="y",
            zaxis_title="z",
            aspectmode="cube"
        )
    )

    PlotlyJS.plot([surface_trace], layout)
end

# Example usage
plot_density_surface_wavefunction(4, 3, 0)



