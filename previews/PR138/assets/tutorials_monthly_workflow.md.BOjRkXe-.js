import{_ as a,o as s,c as n,ao as l}from"./chunks/framework.DCtpE0KF.js";const h=JSON.parse('{"title":"Monthly workflow","description":"","frontmatter":{},"headers":[],"relativePath":"tutorials/monthly_workflow.md","filePath":"tutorials/monthly_workflow.md","lastUpdated":null}'),t={name:"tutorials/monthly_workflow.md"};function i(o,e,p,r,c,u){return s(),n("div",null,[...e[0]||(e[0]=[l(`<h1 id="Monthly-workflow" tabindex="-1">Monthly workflow <a class="header-anchor" href="#Monthly-workflow" aria-label="Permalink to &quot;Monthly workflow {#Monthly-workflow}&quot;">​</a></h1><p>One way to run the model is with long-term monthly climate normals: <a href="/Microclimate.jl/previews/PR138/api#Microclimate.MonthlyMinMaxEnvironment"><code>MonthlyMinMaxEnvironment</code></a> supplies one representative day per month (day of year 15, 46, 74, …, <code>DEFAULT_DAYS</code>), each expanded to an hourly cycle from that month&#39;s minimum/maximum values (see <a href="./../manual/diel_time_handling">Diel curves and time handling</a>). This is driven here by <a href="/Microclimate.jl/previews/PR138/api#Microclimate.example_monthly_weather-Tuple{}"><code>example_monthly_weather</code></a>&#39;s Madison, Wisconsin climate normals (from the CRU CL v2.0 dataset).</p><h2 id="Building-the-run" tabindex="-1">Building the run <a class="header-anchor" href="#Building-the-run" aria-label="Permalink to &quot;Building the run {#Building-the-run}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>using Microclimate, Unitful, CairoMakie</span></span>
<span class="line"><span></span></span>
<span class="line"><span>site = example_site()</span></span>
<span class="line"><span>soil_profile = example_soil_profile()</span></span>
<span class="line"><span>model = MicroModel(;</span></span>
<span class="line"><span>    soil_properties_model = example_soil_properties_model(),</span></span>
<span class="line"><span>    soil_hydraulic_model = example_soil_hydraulic_model(),</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>depths = model.depths</span></span>
<span class="line"><span></span></span>
<span class="line"><span>inputs = MicroInputs(;</span></span>
<span class="line"><span>    site, soil_profile,</span></span>
<span class="line"><span>    environment_minmax = example_monthly_weather(),</span></span>
<span class="line"><span>    environment_daily = example_daily_environment(),</span></span>
<span class="line"><span>    environment_hourly = example_hourly_environment(),</span></span>
<span class="line"><span>    initial_soil_temperature = nothing,  # each representative day resets to its own mean air temperature</span></span>
<span class="line"><span>    initial_soil_moisture = fill(0.42 * 0.25, length(depths)),</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>problem = MicroProblem(model, inputs)</span></span>
<span class="line"><span>out = solve(problem)</span></span>
<span class="line"><span>size(out.soil_temperature)</span></span></code></pre></div><p>Each of the 12 representative days is iterated 3 times by default (<a href="/Microclimate.jl/previews/PR138/api#Microclimate.NonConsecutiveDayMode"><code>NonConsecutiveDayMode</code></a>&#39;s own <code>iterations_per_day</code>), starting fresh each pass from a uniform initial profile (the day&#39;s mean air temperature), to reach a steady-periodic solution before the output is kept — see <a href="./../manual/configuring_the_model">Configuring the model</a> and <a href="/Microclimate.jl/previews/PR138/manual/diel_time_handling#Time-modes">Time modes</a>.</p><h2 id="Soil-temperature-by-depth,-across-the-year" tabindex="-1">Soil temperature by depth, across the year <a class="header-anchor" href="#Soil-temperature-by-depth,-across-the-year" aria-label="Permalink to &quot;Soil temperature by depth, across the year {#Soil-temperature-by-depth,-across-the-year}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>fig = Figure(size = (700, 400))</span></span>
<span class="line"><span>ax = Axis(fig[1, 1]; xlabel = &quot;Hour of representative day (by month)&quot;, ylabel = &quot;Soil temperature (°C)&quot;)</span></span>
<span class="line"><span>depth_labels = [&quot;$(round(Int, ustrip(u&quot;cm&quot;, d))) cm&quot; for d in depths]</span></span>
<span class="line"><span>for (col, label) in zip((1, 4, 7, 10), depth_labels[[1, 4, 7, 10]])</span></span>
<span class="line"><span>    lines!(ax, ustrip.(u&quot;°C&quot;, out.soil_temperature[:, col]); label)</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span>axislegend(ax; position = :rt)</span></span>
<span class="line"><span>fig</span></span></code></pre></div><p>The surface node tracks the diurnal cycle closely; deeper nodes are progressively damped and lagged, as the <a href="./../manual/soil_thermal_properties">heat equation</a> predicts.</p><h2 id="Substrate-type-changes-the-amplitude,-not-just-the-level" tabindex="-1">Substrate type changes the amplitude, not just the level <a class="header-anchor" href="#Substrate-type-changes-the-amplitude,-not-just-the-level" aria-label="Permalink to &quot;Substrate type changes the amplitude, not just the level {#Substrate-type-changes-the-amplitude,-not-just-the-level}&quot;">​</a></h2><p>Swapping <a href="/Microclimate.jl/previews/PR138/api#Microclimate.SoilProfile"><code>SoilProfile</code></a> for a rock rather than a soil substrate lets more heat penetrate deeper, reducing the amplitude near the surface while increasing it at depth. The main lever for this is porosity, i.e. <code>bulk_density</code> relative to <code>mineral_density</code>: soil&#39;s default (<code>bulk_density = 1.3 Mg/m^3</code>, <code>mineral_density = 2.56 Mg/m^3</code>) is about half pore space, while solid rock has almost none, so <code>bulk_density</code> needs to sit close to <code>mineral_density</code>. Rock&#39;s own mineral conductivity and heat capacity are also higher than typical soil minerals, and with almost no pore space, there&#39;s almost no room for soil moisture either:</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>rock_profile = example_soil_profile(;</span></span>
<span class="line"><span>    bulk_density = 2.55u&quot;Mg/m^3&quot;, mineral_density = 2.6u&quot;Mg/m^3&quot;,</span></span>
<span class="line"><span>    mineral_conductivity = 2.5u&quot;W/m/K&quot;, mineral_heat_capacity = 750.0u&quot;J/kg/K&quot;,</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>rock_inputs = MicroInputs(;</span></span>
<span class="line"><span>    site, soil_profile = rock_profile,</span></span>
<span class="line"><span>    environment_minmax = example_monthly_weather(),</span></span>
<span class="line"><span>    environment_daily = example_daily_environment(),</span></span>
<span class="line"><span>    environment_hourly = example_hourly_environment(),</span></span>
<span class="line"><span>    initial_soil_temperature = nothing,</span></span>
<span class="line"><span>    initial_soil_moisture = fill(0.01, length(depths)),</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>rock_out = solve(MicroProblem(model, rock_inputs))</span></span>
<span class="line"><span></span></span>
<span class="line"><span>fig2 = Figure(size = (700, 400))</span></span>
<span class="line"><span>ax2 = Axis(fig2[1, 1]; xlabel = &quot;Hour (July, day 196)&quot;, ylabel = &quot;Surface temperature (°C)&quot;)</span></span>
<span class="line"><span>july = (7 - 1) * 24 + 1:7 * 24</span></span>
<span class="line"><span>lines!(ax2, ustrip.(u&quot;°C&quot;, out.soil_temperature[july, 1]); label = &quot;soil&quot;)</span></span>
<span class="line"><span>lines!(ax2, ustrip.(u&quot;°C&quot;, rock_out.soil_temperature[july, 1]); label = &quot;rock&quot;)</span></span>
<span class="line"><span>axislegend(ax2)</span></span>
<span class="line"><span>fig2</span></span></code></pre></div><h2 id="Next-steps" tabindex="-1">Next steps <a class="header-anchor" href="#Next-steps" aria-label="Permalink to &quot;Next steps {#Next-steps}&quot;">​</a></h2><ul><li><p><a href="./daily_hourly_workflow">Daily/hourly workflow</a> — driving the model from real hourly weather-station data instead of monthly normals.</p></li><li><p><a href="./../manual/configuring_the_model">Configuring the model</a> — every other slot on <code>MicroModel</code> this tutorial left at its default.</p></li></ul>`,13)])])}const m=a(t,[["render",i]]);export{h as __pageData,m as default};
