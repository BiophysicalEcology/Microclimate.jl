module MicroclimateMakieExt

using Microclimate
using Unitful
using Makie

function Makie.plot(mr::Microclimate.MicroResult; alpha=0.7)
    fig = Figure()
    solar_ax = Axis(fig[1, 1])
    lines!(solar_ax, mr.global_solar; label="Global Solar Radiation", alpha)
    lines!(solar_ax, mr.diffuse_solar; label="Difuse Solar Radiation", alpha)
    lines!(solar_ax, mr.direct_solar; label="Direct Solar Radiation", alpha)
    axislegend(solar_ax)

    temperature_ax = Axis(fig[2, 1])
    lines!(temperature_ax, vec(mr.air_temperature); label="Air Temperature", alpha)
    lines!(temperature_ax, vec(mr.sky_temperature); label="Sky Temperature", alpha)
    labels = ["Soil depth $i" for i in 1:size(mr.soil_temperature, 1)]
    color = [Makie.RGBA(rand(RGBf), alpha) for _ in 1:size(mr.soil_temperature, 1)]
    series!(temperature_ax, ustrip.(mr.soil_temperature'); labels, color)
    axislegend(temperature_ax)

    soil_moisture_ax = Axis(fig[3, 1])
    labels = ["Soil moisture $i" for i in 1:size(mr.soil_moisture, 1)]
    color = [Makie.RGBA(rand(RGBf), alpha) for _ in 1:size(mr.soil_moisture, 1)]
    series!(soil_moisture_ax, ustrip.(mr.soil_moisture'); labels, color)
    axislegend(soil_moisture_ax)

    linkxaxes!(solar_ax, temperature_ax, soil_moisture_ax)

    return fig
end

end
