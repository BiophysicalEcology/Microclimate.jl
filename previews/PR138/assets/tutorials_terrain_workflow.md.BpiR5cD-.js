import{_ as s,o as n,c as e,ao as l}from"./chunks/framework.DCtpE0KF.js";const h=JSON.parse('{"title":"Terrain","description":"","frontmatter":{},"headers":[],"relativePath":"tutorials/terrain_workflow.md","filePath":"tutorials/terrain_workflow.md","lastUpdated":null}'),t={name:"tutorials/terrain_workflow.md"};function p(i,a,o,r,u,c){return n(),e("div",null,[...a[0]||(a[0]=[l(`<h1 id="Terrain" tabindex="-1">Terrain <a class="header-anchor" href="#Terrain" aria-label="Permalink to &quot;Terrain {#Terrain}&quot;">​</a></h1><p>Slope, aspect and hillshade change the radiation a surface receives, and therefore its soil temperature — see <a href="./../manual/radiation">Radiation</a> for the equations. This tutorial runs the same Madison, Wisconsin site used elsewhere in these docs (<a href="/Microclimate.jl/previews/PR138/api#Microclimate.example_site-Tuple{}"><code>example_site</code></a>) under a few different terrain configurations.</p><h2 id="Slope-and-aspect" tabindex="-1">Slope and aspect <a class="header-anchor" href="#Slope-and-aspect" aria-label="Permalink to &quot;Slope and aspect {#Slope-and-aspect}&quot;">​</a></h2><p>The default is flat ground (<code>slope = 0°</code>). Here&#39;s the effect of a 45° south-facing slope, compared to flat ground, on the 0%-shade simulation:</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>using Microclimate, Unitful, CairoMakie</span></span>
<span class="line"><span></span></span>
<span class="line"><span>site = example_site()</span></span>
<span class="line"><span>soil_profile = example_soil_profile()</span></span>
<span class="line"><span>model = MicroModel(;</span></span>
<span class="line"><span>    soil_properties_model = example_soil_properties_model(),</span></span>
<span class="line"><span>    soil_hydraulic_model = example_soil_hydraulic_model(),</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>depths = model.depths</span></span>
<span class="line"><span></span></span>
<span class="line"><span>function run_at(site)</span></span>
<span class="line"><span>    inputs = MicroInputs(;</span></span>
<span class="line"><span>        site, soil_profile,</span></span>
<span class="line"><span>        environment_minmax = example_monthly_weather(),</span></span>
<span class="line"><span>        environment_daily = example_daily_environment(),</span></span>
<span class="line"><span>        environment_hourly = example_hourly_environment(),</span></span>
<span class="line"><span>        initial_soil_temperature = nothing,  # each representative day resets to its own mean air temperature</span></span>
<span class="line"><span>        initial_soil_moisture = fill(0.105, length(depths)),</span></span>
<span class="line"><span>    )</span></span>
<span class="line"><span>    solve(MicroProblem(model, inputs))</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span></span></span>
<span class="line"><span>flat_out = run_at(site)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>sloped_site = example_site(; slope = 45.0u&quot;°&quot;, aspect = 180.0u&quot;°&quot;)</span></span>
<span class="line"><span>sloped_out = run_at(sloped_site)</span></span>
<span class="line"><span>nothing # hide</span></span></code></pre></div><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>july = (7 - 1) * 24 + 1:7 * 24</span></span>
<span class="line"><span>december = (12 - 1) * 24 + 1:12 * 24</span></span>
<span class="line"><span></span></span>
<span class="line"><span>fig = Figure(size = (700, 600))</span></span>
<span class="line"><span>ax1 = Axis(fig[1, 1]; ylabel = &quot;Surface temperature (°C)&quot;, title = &quot;July (day 196)&quot;)</span></span>
<span class="line"><span>lines!(ax1, ustrip.(u&quot;°C&quot;, flat_out.soil_temperature[july, 1]); label = &quot;flat&quot;)</span></span>
<span class="line"><span>lines!(ax1, ustrip.(u&quot;°C&quot;, sloped_out.soil_temperature[july, 1]); label = &quot;45° south-facing&quot;)</span></span>
<span class="line"><span>axislegend(ax1)</span></span>
<span class="line"><span>ax2 = Axis(fig[2, 1]; xlabel = &quot;Hour&quot;, ylabel = &quot;Surface temperature (°C)&quot;, title = &quot;December (day 349)&quot;)</span></span>
<span class="line"><span>lines!(ax2, ustrip.(u&quot;°C&quot;, flat_out.soil_temperature[december, 1]); label = &quot;flat&quot;)</span></span>
<span class="line"><span>lines!(ax2, ustrip.(u&quot;°C&quot;, sloped_out.soil_temperature[december, 1]); label = &quot;45° south-facing&quot;)</span></span>
<span class="line"><span>axislegend(ax2)</span></span>
<span class="line"><span>fig</span></span></code></pre></div><p>At this latitude, with the sun low in the sky in winter, the slope&#39;s effect is strongest on the December temperatures — a south-facing slope catches the low winter sun far more directly than flat ground does.</p><h2 id="Hillshade" tabindex="-1">Hillshade <a class="header-anchor" href="#Hillshade" aria-label="Permalink to &quot;Hillshade {#Hillshade}&quot;">​</a></h2><p>Gullies and gorges receive a shorter window of direct sun and have a smaller view of the sky, so they don&#39;t cool down as much at night. <code>Site.horizon_angles</code> sets this: a vector of horizon angles (the angle to the nearest obstruction, starting due north and running clockwise at equal steps — 32 directions, 11.25° apart, matching <a href="https://github.com/Deltares/Geomorphometry.jl" target="_blank" rel="noreferrer">Geomorphometry.jl</a>&#39;s convention), which also determine <code>sky_view_fraction</code> (the fraction of the sky not blocked by terrain — see <a href="./../manual/radiation">Radiation</a>). Here&#39;s a steep-sided gully running north-south — open along the north-south axis, walled in to the east and west:</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>gully_horizon_angles = vcat(</span></span>
<span class="line"><span>    fill(0.0, 2), fill(65.0, 14),</span></span>
<span class="line"><span>    fill(0.0, 2), fill(65.0, 14),</span></span>
<span class="line"><span>)u&quot;°&quot;</span></span>
<span class="line"><span>gully_sky_view_fraction = 1 - sum(sin, gully_horizon_angles) / length(gully_horizon_angles)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>gully_site = example_site(;</span></span>
<span class="line"><span>    horizon_angles = gully_horizon_angles, sky_view_fraction = gully_sky_view_fraction,</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span>gully_out = run_at(gully_site)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>fig2 = Figure()</span></span>
<span class="line"><span>ax = Axis(fig2[1, 1]; xlabel = &quot;Hour (July)&quot;, ylabel = &quot;Surface temperature (°C)&quot;)</span></span>
<span class="line"><span>lines!(ax, ustrip.(u&quot;°C&quot;, flat_out.soil_temperature[july, 1]); label = &quot;open&quot;)</span></span>
<span class="line"><span>lines!(ax, ustrip.(u&quot;°C&quot;, gully_out.soil_temperature[july, 1]); label = &quot;north-south gully&quot;)</span></span>
<span class="line"><span>axislegend(ax)</span></span>
<span class="line"><span>fig2</span></span></code></pre></div><p>Most of the day the gully floor is shaded and stays well below the open-ground temperature, but around solar noon the sun briefly aligns with the narrow north-south opening and direct beam pours in — a sharp, short-lived spike rather than a gradual warming, since the sun crosses from fully blocked to fully visible within about an hour at this resolution.</p><h2 id="Next-steps" tabindex="-1">Next steps <a class="header-anchor" href="#Next-steps" aria-label="Permalink to &quot;Next steps {#Next-steps}&quot;">​</a></h2><ul><li><p><a href="./../manual/radiation">Radiation</a> for the shortwave slope/aspect and longwave hillshade/view-factor equations behind both examples.</p></li><li><p><a href="./../manual/configuring_the_model">Configuring the model</a> — every other slot on <code>MicroModel</code> this tutorial left at its default.</p></li></ul>`,13)])])}const _=s(t,[["render",p]]);export{h as __pageData,_ as default};
