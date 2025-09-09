using Microclimate
using NCDatasets
using Unitful
using Plots

# extract monthly data and simulate middle day of each year
x = [138.7, -34.27]
folder = "c:\\globalclimate"
(; elevation, precipitation, rainy_days, wind_speed,
            air_temperature_min, air_temperature_max, humidity_min, 
            humidity_max, cloud_cover) = load_CRU_CL_v2(; x, folder)

micro_out = runmicro(; 
    latitude = x[2]u"°",
    elevation,
    daily_rainfall = precipitation,
    wind_speed_max = wind_speed,
    wind_speed_min = wind_speed * 0.1,
    air_temperature_min,
    air_temperature_max,
    cloud_cover_min = cloud_cover,
    cloud_cover_max = cloud_cover,
)

plot(u"°C".(micro_out.soil_temperature), legend = false)
