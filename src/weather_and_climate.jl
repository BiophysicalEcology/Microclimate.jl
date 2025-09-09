# build up library of functions for extracting climate/weather data from various
# sources and prepare them as forcing variables for the microclimate model
function load_CRU_CL_v2(;
    x::Vector{Float64},
    folder::String
)
    @assert length(x) == 2 "x must be a 2-element vector [lon, lat]"
    lon, lat = x
    global_climate_file = joinpath(folder, "global_climate.nc")
    ds = NCDataset(global_climate_file)
    days_per_month = [31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0]u"d"

    # --- extract data at location (nearest neighbor here) ---
    lon_vals = ds["longitude"][:]
    lat_vals = ds["latitude"][:]

    i_lon = findmin(abs.(lon_vals .- lon))[2]
    i_lat = findmin(abs.(lat_vals .- lat))[2]

    climate = vec(ds["variable"][i_lon, i_lat, :]) # <-- replace "varname" with correct variable

    # pull out variables
    elevation = climate[1]u"m"
    precipitation = climate[2:13]u"kg/m^2" ./ days_per_month
    rainy_days = climate[14:25] ./ 10
    wind_speed = (climate[26:37] ./ 10)u"m/s"
    air_temperature_min = (climate[38:49] ./ 10)u"°C"
    air_temperature_max = (climate[50:61] ./ 10)u"°C"
    humidity_min = climate[62:73] ./ 10
    humidity_max = climate[74:85] ./ 10
    cloud_cover = climate[86:97] ./ 10

    close(ds)

    return (; elevation, precipitation, rainy_days, wind_speed,
            air_temperature_min, air_temperature_max, humidity_min, 
            humidity_max, cloud_cover)
end