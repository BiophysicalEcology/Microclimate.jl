import{_ as s,o as a,c as e,ao as t}from"./chunks/framework.DniuIMox.js";const m=JSON.parse('{"title":"Snow demo","description":"","frontmatter":{},"headers":[],"relativePath":"tutorials/snow_demo.md","filePath":"tutorials/snow_demo.md","lastUpdated":null}'),l={name:"tutorials/snow_demo.md"};function p(i,n,o,r,c,d){return a(),e("div",null,[...n[0]||(n[0]=[t(`<h1 id="Snow-demo" tabindex="-1">Snow demo <a class="header-anchor" href="#Snow-demo" aria-label="Permalink to &quot;Snow demo {#Snow-demo}&quot;">​</a></h1><p>This demonstrates <a href="/Microclimate.jl/dev/api#Microclimate.SnowModel"><code>SnowModel</code></a> buffering soil temperature through a winter onset — the mechanism behind Kearney (2020)&#39;s finding that shallow soil warms at roughly a third the rate of air temperature at snow-affected sites (see <a href="/Microclimate.jl/dev/manual/snow#Background">Snow</a>).</p><h2 id="A-synthetic-cold-winter" tabindex="-1">A synthetic cold winter <a class="header-anchor" href="#A-synthetic-cold-winter" aria-label="Permalink to &quot;A synthetic cold winter {#A-synthetic-cold-winter}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>using Microclimate, Unitful, CairoMakie</span></span>
<span class="line"><span></span></span>
<span class="line"><span>site = example_site()</span></span>
<span class="line"><span>depths = Microclimate.DEFAULT_DEPTHS</span></span>
<span class="line"><span>n = length(depths)</span></span>
<span class="line"><span>days = collect(335:15:365)   # late Nov through Dec, one representative day per fortnight</span></span>
<span class="line"><span></span></span>
<span class="line"><span>reference_temperature_min = (range(-2.0, -18.0; length = length(days)) .+ 273.15)u&quot;K&quot;</span></span>
<span class="line"><span>reference_temperature_max = (range(4.0, -8.0; length = length(days)) .+ 273.15)u&quot;K&quot;</span></span>
<span class="line"><span>environment_minmax = MonthlyMinMaxEnvironment(; forcings = minmax_forcings(;</span></span>
<span class="line"><span>    reference_temperature_min, reference_temperature_max,</span></span>
<span class="line"><span>    reference_wind_speed_min = fill(0.5, length(days))u&quot;m/s&quot;,</span></span>
<span class="line"><span>    reference_wind_speed_max = fill(4.0, length(days))u&quot;m/s&quot;,</span></span>
<span class="line"><span>    reference_humidity_min = fill(0.6, length(days)),</span></span>
<span class="line"><span>    reference_humidity_max = fill(0.95, length(days)),</span></span>
<span class="line"><span>    cloud_cover_min = fill(0.3, length(days)),</span></span>
<span class="line"><span>    cloud_cover_max = fill(0.6, length(days)),</span></span>
<span class="line"><span>))</span></span>
<span class="line"><span>environment_daily = example_daily_environment(days; rainfall = fill(15.0, length(days))u&quot;kg/m^2&quot;)</span></span>
<span class="line"><span>environment_hourly = example_hourly_environment(days; elevation = site.elevation)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>snow_model = SnowModel(; snow_temperature_threshold = 1.0u&quot;°C&quot;)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>function run(; snow_model)</span></span>
<span class="line"><span>    model = MicroModel(;</span></span>
<span class="line"><span>        soil_properties_model = example_soil_properties_model(),</span></span>
<span class="line"><span>        soil_hydraulic_model = example_soil_hydraulic_model(),</span></span>
<span class="line"><span>        snow_model, depths,</span></span>
<span class="line"><span>    )</span></span>
<span class="line"><span>    inputs = MicroInputs(;</span></span>
<span class="line"><span>        site, soil_profile = example_soil_profile(depths),</span></span>
<span class="line"><span>        environment_minmax, environment_daily, environment_hourly,</span></span>
<span class="line"><span>        initial_soil_temperature = fill(u&quot;K&quot;(2.0u&quot;°C&quot;), n),</span></span>
<span class="line"><span>        initial_soil_moisture = fill(0.2, n),</span></span>
<span class="line"><span>    )</span></span>
<span class="line"><span>    solve(MicroProblem(model, inputs; days, time_mode = ConsecutiveDayMode()))</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span></span></span>
<span class="line"><span>with_snow = run(; snow_model)</span></span>
<span class="line"><span>without_snow = run(; snow_model = NoSnow())</span></span>
<span class="line"><span>nothing # hide</span></span></code></pre></div><h2 id="Snow-depth-and-its-effect-on-soil-temperature" tabindex="-1">Snow depth and its effect on soil temperature <a class="header-anchor" href="#Snow-depth-and-its-effect-on-soil-temperature" aria-label="Permalink to &quot;Snow depth and its effect on soil temperature {#Snow-depth-and-its-effect-on-soil-temperature}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>fig = Figure(size = (700, 500))</span></span>
<span class="line"><span>ax1 = Axis(fig[1, 1]; ylabel = &quot;Snow depth (cm)&quot;)</span></span>
<span class="line"><span>lines!(ax1, ustrip.(u&quot;cm&quot;, with_snow.snow_depth))</span></span>
<span class="line"><span>ax2 = Axis(fig[2, 1]; xlabel = &quot;Hour&quot;, ylabel = &quot;Surface soil temperature (°C)&quot;)</span></span>
<span class="line"><span>lines!(ax2, ustrip.(u&quot;°C&quot;, with_snow.soil_temperature[:, 1]); label = &quot;with snow&quot;)</span></span>
<span class="line"><span>lines!(ax2, ustrip.(u&quot;°C&quot;, without_snow.soil_temperature[:, 1]); label = &quot;without snow (rain instead)&quot;)</span></span>
<span class="line"><span>axislegend(ax2; position = :rb)</span></span>
<span class="line"><span>fig</span></span></code></pre></div><p>As the pack builds, surface soil temperature decouples from the falling air temperature and settles close to 0°C, buffered by the latent heat of fusion at the snow-soil interface — while the no-snow run (precipitation still falls, as rain) continues to track the cooling air directly.</p><h2 id="Next-steps" tabindex="-1">Next steps <a class="header-anchor" href="#Next-steps" aria-label="Permalink to &quot;Next steps {#Next-steps}&quot;">​</a></h2><ul><li><p><a href="./../manual/snow">Snow</a> for the density, albedo, conductivity and melt equations behind this behaviour.</p></li><li><p><a href="./../manual/soil_thermal_properties">Soil thermal properties</a> for the underlying heat equation and phase-transition treatment near 0°C.</p></li></ul>`,9)])])}const h=s(l,[["render",p]]);export{m as __pageData,h as default};
