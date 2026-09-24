using Microclimate
using Microclimate: shortwave_radiation!
using Unitful
using Test

@testset "AngstromMaxwellShortwave: clearness index invariant to zenith angle alone" begin
    # Regression test: eth (extraterrestrial irradiance) previously omitted the
    # cos(zenith) horizontal-plane projection, so diffuse_fraction spuriously
    # rose toward 1 away from solar noon even under fixed cloud cover.
    n = 4
    zenith = [0.0, 30.0, 60.0, 80.0]u"°"
    cloud = fill(0.3, n)
    doy = fill(172, n)
    # Clear-sky irradiance that itself follows the physical cos(zenith) law --
    # with the bug fixed, diffuse_fraction should come out ~equal across zenith.
    diffuse_clear_sky = 50.0u"W/m^2" .* max.(0.0, cosd.(zenith))
    direct_clear_sky = 800.0u"W/m^2" .* max.(0.0, cosd.(zenith))

    output = (; solar_radiation = (; global_horizontal = zeros(typeof(0.0u"W/m^2"), n)),
                diffuse_fraction = zeros(n))
    result = shortwave_radiation!(AngstromMaxwellShortwave(), output, cloud, diffuse_clear_sky, direct_clear_sky, zenith, doy)

    @test maximum(result.diffuse_fraction) - minimum(result.diffuse_fraction) < 0.05
end

@testset "diffuse_fraction: not pinned near 1 at solar noon under moderate cloud" begin
    out = solve(example_microclimate_problem())
    july_noon = (7 - 1) * 24 + 13  # day 7 (July), hour 12
    @test out.diffuse_fraction[july_noon] < 0.7
end
