"""
    gads(latitude, longitude, relhum, season, data_root)

Function to extract data from the Global Aerosol Data Set (originally a Fortran program),
for a specific location, season and relative humidity

Arguments
- latitude, longitude::Float64 : input location
- relhum::Int64    : humidity class index 1:8 corresponding to [0,50,70,80,90,95,98,99]%
- season::Int64    : 1 = winter, 0 = summer
- data_root::String : path to gads data folder

Returns
- optdep::Matrix{Float64} of size (25,2): [wavelength_nm, optical_depth]

Reference
Koepke, P., Hess, M., Schult, I., & Shettle, E. P. (1997). Global Aerosol Data Set. Max-Planck-Institut für Meteorologie.
http://www.mpimet.mpg.de/fileadmin/publikationen/Reports/MPI-Report_243.pdf
"""
function gads(
    # TODO these could be keywords
    latitude::Float64=43.07305, 
    longitude::Float64=89.40123, 
    relhum::Int64=1, 
    season::Int64=0, 
    data_root::String=realpath(joinpath(dirname(pathof(Microclimate)), "../data")),
)
    _substr(s::AbstractString, i::Int, j::Int) = begin
        i < 1 && return ""
        j < i && return ""
        last_idx = lastindex(s)   # rename from n -> last_idx
        i > last_idx && return ""
        j = min(j, last_idx)
        try
            return s[i:j]
        catch
            # Fallback by iterating character indices (handles UTF-8 widths)
            inds = collect(eachindex(s))
            i > length(inds) && return ""
            j = min(j, length(inds))
            return s[inds[i]:inds[j]]
        end
    end

    _parsefloat(s) = try
        parse(Float64, strip(s))
    catch
        NaN
    end

    _parseint(s) = try
        parse(Int, strip(s))
    catch
        Int(typemax(Int)) # sentinel unlikely to match
    end

    _split_numbers(line::AbstractString) = begin
        parts = split(strip(line))
        v = Float64[]
        for p in parts
            try
                push!(v, parse(Float64, p))
            catch
                # ignore non-numeric tokens
            end
        end
        return v
    end

    lat5s = collect(-90:5:90)
    lon5s = collect(-180:5:175)
    lat5 = lat5s[argmin(abs.(lat5s .- latitude))]
    lon5 = lon5s[argmin(abs.(lon5s .- longitude))]

    optdep = zeros(25, 2)  # output array

    # -------- scalars / flags --------
    niw = 1; njh = 1; nih = 1
    lata = lat5; late = lat5; lati = 5
    lona = lon5; lone = lon5; loni = 5
    nlmal = 1; norm = 1; nprog = 4
    ip = 0; iwel = 1; ihum = relhum
    il = 1; ih = 1
    ilat = lata; imal = 1; ilon = lona
    nwel = 25

    jnopar = [1,1,1,1,1,1,0,0,1,0,0,0,0]
    nop = 7
    kop = nop

    optnam = ["ext.coef","sca.coef","abs.coef","sisc.alb","asym.par",
              "op.depth","        ","turb.fac","li.ratio","pha.func",
              "ext.rat ","abs.rat ","        "]
    opanam = fill("", 13)

    alamb = [0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,
             0.9,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.2,3.39,3.5,3.75,
             4.0,4.5,5.0,5.5,6.0,6.2,6.5,7.2,7.9,8.2,8.5,8.7,9.0,
             9.2,9.5,9.8,10.0,10.6,11.0,11.5,12.5,13.0,14.0,14.8,
             15.0,16.4,17.2,18.0,18.5,20.0,21.3,22.5,25.0,27.9,30.0,
             35.0,40.0]
    mlamb = 61

    ahum = [0,50,70,80,90,95,98,99]
    nhum = [0,50,70,80,90,95,98,99]
    mhum = 8
    chum = ["00","50","70","80","90","95","98","99"]
    comnam = ["inso","waso","soot","ssam","sscm","minm","miam",
              "micm","mitr","suso","stco","stma","cucc","cucp",
              "cuma","fog-","cir1","cir2","cir3","    "]

    atn = ["",""]
    pat = ["",""]
    rht = [0.0, 0.0]
    n = [0.0, 0.0]
    nh = [0, 0]
    njc = [0, 0]

    acnr = zeros(5,2)   # component IDs per layer
    acmr = zeros(5,2)   # mixing ratios per layer
    khum = zeros(Int,8)

    ext = zeros(Float64,2,5); sca = zeros(Float64,2,5); abs1 = zeros(Float64,2,5)
    sis = zeros(Float64,2,5); asy = zeros(Float64,2,5); bac = zeros(Float64,2,5)
    pha = zeros(Float64,112,2,5)
    bre = zeros(Float64,2,5); bim = zeros(Float64,2,5)

    # buffers (size 20)
    kbuf   = zeros(Int,20)
    extbuf = zeros(Float64,20); scabuf = zeros(Float64,20); absbuf = zeros(Float64,20)
    sisbuf = zeros(Float64,20); asybuf = zeros(Float64,20); bacbuf = zeros(Float64,20)
    phabuf = zeros(Float64,112,20)
    brebuf = zeros(Float64,20); bimbuf = zeros(Float64,20)

    extn = zeros(2); absn = zeros(2); scan = zeros(2)
    pf18n = zeros(2); supf = zeros(112); phafu = zeros(112,2)
    exta = zeros(2); absa = zeros(2); scaa = zeros(2)
    ssa = zeros(2); asf = zeros(2); pf18a = zeros(2)
    scar = zeros(2); absr = zeros(2); omer = zeros(2)

    jnangle = zeros(Int,112)  # not used under current jnopar flags
    ncomp = zeros(Int,10)
    oparam = zeros(Float64,10,2)
    nltyp = zeros(Int,10)
    parlay = zeros(10,2)
    boundl = zeros(10)
    boundu = zeros(10)

    # height profiles and constants
    iip = 7
    nil = [1,1,1,1,2,1,2,0,0]
    hfta = [10.0, 2.0, 10.0, 10.0, 8.5, 6.0, 8.5, 0.0, 0.0, 0.0]
    hstra = vcat(fill(23.0,7), [0.0,0.0,0.0])
    h0 = zeros(2,10)
    h1 = zeros(2,10)
    hm = zeros(2,10)
    h0[2, [5,7]] .= 2.0
    h1[1, 1:7] .= [2.0, 10.0, 2.0, 2.0, 2.0, 6.0, 2.0]
    h1[2, [5,7]] .= 3.5
    hm[1, 1:7] .= [2.0, 5.7, 1.77, 0.86, 0.86, 1.9, 1.77]
    hm[2, [5,7]] .= 1.5

    # extback.dat: columns: (lamb, extfta, extstr)
    extback_path = joinpath(data_root, "extback.dat")
    extback = readdlm(extback_path, skipstart=2)
    extfta = vec(extback[:,2])
    extstr = vec(extback[:,3])

    # determine active parameter names
    for i in 1:13
        if jnopar[i] == 1
            ip += 1
            opanam[ip] = optnam[i]
        end
    end

    opnam = fill(0, kop)  # numeric placeholder to mirror R (unused downstream)
    for i in 1:kop
        opnam[i] = 0
    end

    ws = season == 1 ? 'w' : 's'
    cseas = ws == 'w' ? "winter " : "summer "

    # seasonal file: winter.dat OR summer.dat
    datfile = ws == 'w' ? joinpath(data_root, "glodat", "winter.dat") : joinpath(data_root, "glodat", "summer.dat")
    ntape = readlines(datfile)

    # find the block for (ilat, ilon)
    k = 2
    line_read = ""
    for kk in 2:length(ntape)
        lr = ntape[kk]
        latx = _parsefloat(_substr(lr, 2, 4))
        lonx = _parsefloat(_substr(lr, 5, 8))
        if !isnan(latx)
            if Int(round(latx)) == ilat && Int(round(lonx)) == ilon
                k = kk
                line_read = lr
                break
            end
        end
    end

    # main wavelength loop (iwel = 1:nwel)
    for iwel in 1:nwel
        # ilamb maps 1:25 to the first 25 wavelengths in `alamb`
        ilamb = iwel
        # ------------- Read optical raw data from current line (and following, if needed) -------------
        # R substr ranges are inclusive; replicate them exactly.
        nl   = _parseint(_substr(line_read, 9, 11))
        prnr = _parseint(_substr(line_read, 12, 14))
        atn[1] = _substr(line_read, 15, 17)
        pat[1] = _substr(line_read, 18, 20)
        rht[1] = _parsefloat(_substr(line_read, 21, 23))
        n[1]   = _parsefloat(_substr(line_read, 24, 33))
        njc[1] = _parseint(_substr(line_read, 34, 37))
        acnr[1,1] = _parseint(_substr(line_read, 38, 40))
        acmr[1,1] = _parsefloat(_substr(line_read, 41, 50))
        acnr[2,1] = _parseint(_substr(line_read, 51, 53))
        acmr[2,1] = _parsefloat(_substr(line_read, 54, 63))
        acnr[3,1] = _parseint(_substr(line_read, 64, 66))
        acmr[3,1] = _parsefloat(_substr(line_read, 67, 76))

        # additional components (if njc[1] > 3) span additional lines
        k2 = k
        if njc[1] > 3
            for jc in 4:njc[1]
                k2 += 1
                line_read2 = ntape[k2]
                acnr[jc,1] = _parseint(_substr(line_read2, 38, 40))
                acmr[jc,1] = _parsefloat(_substr(line_read2, 41, 50))
            end
        end

        # if more than one layer (nl != 1), parse subsequent layer lines
        if nl != 1
            for l in 2:nl
                k2 += 1
                line_read2 = ntape[k2]
                atn[l]   = _substr(line_read2, 15, 17)
                pat[l]   = _substr(line_read2, 18, 20)
                rht[l]   = _parsefloat(_substr(line_read2, 21, 23))
                n[l]     = _parsefloat(_substr(line_read2, 24, 33))
                njc[l]   = _parseint(_substr(line_read2, 34, 37))
                acnr[1,l] = _parseint(_substr(line_read2, 38, 40))
                acmr[1,l] = _parsefloat(_substr(line_read2, 41, 50))
                acnr[2,l] = _parseint(_substr(line_read2, 51, 53))
                acmr[2,l] = _parsefloat(_substr(line_read2, 54, 63))
                acnr[3,l] = _parseint(_substr(line_read2, 64, 66))
                acmr[3,l] = _parsefloat(_substr(line_read2, 67, 76))
                if njc[l] > 3
                    for jc in 4:njc[l]
                        k2 += 1
                        line_read3 = ntape[k2]
                        acnr[jc,l] = _parseint(_substr(line_read3, 38, 40))
                        acmr[jc,l] = _parsefloat(_substr(line_read3, 41, 50))
                    end
                end
            end
        end

        # validate mixing ratios sum per layer (tolerate 0.01)
        for l in 1:nl
            s = 0.0
            for ic in 1:njc[l]
                s += acmr[ic,l]
            end
            if abs(s - 1.0) >= 0.01
                @warn "sum of mixing ratios is not 1. please check input files" layer=l sum=s
            end
        end

        # determine humidity class per layer -> nh[l]
        for ilayer in 1:nl
            rhv = rht[ilayer]
            if rhv < 30
                nh[ilayer] = 1
            elseif 30 < rhv <= 65
                nh[ilayer] = 2
            elseif 65 < rhv <= 75
                nh[ilayer] = 3
            elseif 75 < rhv <= 85
                nh[ilayer] = 4
            elseif 85 < rhv <= 92
                nh[ilayer] = 5
            elseif 92 < rhv <= 97
                nh[ilayer] = 6
            elseif rhv == 98
                nh[ilayer] = 7
            elseif rhv == 99
                nh[ilayer] = 8
            else
                nh[ilayer] = 1 # default safeguard
            end
        end

        # ---------- Calculation of optical parameters at the current grid point ----------
        ibuf = 0
        for ilayer in 1:nl
            # (nih == 0) branch not active here (nih==1), but keep structure
            if nih == 0
                for ihu in 1:mhum
                    if nh[ilayer] == ihu
                        khum[ihu] = ihu
                    end
                end
            end

            for ic in 1:njc[ilayer]
                jc = Int(round(acnr[ic, ilayer]))
                # Exclusion: swelling at insoluble, soot and mineral components and at clouds
                if jc == 1 || jc == 3 || (6 <= jc <= 9) || jc > 10
                    iht = 1
                    nta = 700 + (jc * 10) + 1
                else
                    iht = ihum
                    nta = 700 + (jc * 10) + iht
                end

                tap = joinpath(data_root, "optdat", string(comnam[jc], chum[iht]))

                # buffer lookup key for (ilamb, nta)
                nbuf = ilamb * 1000 + nta
                exist_chk = false
                mbuf = 0
                for ib in 1:20
                    if nbuf == kbuf[ib]
                        exist_chk = true
                        mbuf = ib
                        break
                    end
                end

                if exist_chk
                    ext[ilayer,ic] = extbuf[mbuf]
                    sca[ilayer,ic] = scabuf[mbuf]
                    abs1[ilayer,ic] = absbuf[mbuf]
                    sis[ilayer,ic] = sisbuf[mbuf]
                    asy[ilayer,ic] = asybuf[mbuf]
                    bac[ilayer,ic] = bacbuf[mbuf]
                    bre[ilayer,ic] = brebuf[mbuf]
                    bim[ilayer,ic] = bimbuf[mbuf]
                else
                    # rotate buffer index
                    if ibuf < 20
                        ibuf += 1
                    else
                        ibuf = 1
                    end
                    kbuf[ibuf] = nbuf

                    # read component file lines
                    ntap = readlines(tap)
                    # find the marker line: "# optical parameters:"
                    iline = 0
                    for (idx, ln) in enumerate(ntap)
                        if strip(ln) == "# optical parameters:"
                            iline = idx
                            break
                        end
                    end
                    if iline == 0
                        error("Marker '# optical parameters:' not found in $(tap)")
                    end

                    # wavelength table starts at iline+6 for mlamb rows
                    rlamb = NaN
                    extco = NaN
                    scaco = NaN
                    absco = NaN
                    sisca = NaN
                    asymf = NaN
                    exn_local = NaN
                    refr = NaN
                    refi = NaN

                    ila_last = iline + 5
                    for ila in (iline + 6):(iline + 5 + mlamb)
                        ntap_out = ntap[ila]
                        rlamb  = _parsefloat(_substr(ntap_out, 3, 12))
                        extco  = _parsefloat(_substr(ntap_out, 13, 22))
                        scaco  = _parsefloat(_substr(ntap_out, 23, 32))
                        absco  = _parsefloat(_substr(ntap_out, 33, 42))
                        sisca  = _parsefloat(_substr(ntap_out, 43, 52))
                        asymf  = _parsefloat(_substr(ntap_out, 53, 62))
                        exn_local = _parsefloat(_substr(ntap_out, 63, 72))
                        refr   = _parsefloat(_substr(ntap_out, 73, 83))
                        refi   = _parsefloat(_substr(ntap_out, 84, 94))

                        if rlamb == alamb[ilamb]
                            ext[ilayer,ic] = extco
                            sca[ilayer,ic] = scaco
                            abs1[ilayer,ic] = absco
                            sis[ilayer,ic] = sisca
                            asy[ilayer,ic] = asymf
                            bre[ilayer,ic] = refr
                            bim[ilayer,ic] = refi
                        end
                        ila_last = ila
                    end

                    # scattering phase function angles
                    ntheta = length(ntap) - (iline + 5 + 8 + mlamb)
                    # build temporary numeric matrix, but we only store column 2 later
                    # angles start at (ila_last + 8 + 1) .. (ila_last + 8 + ntheta)
                    for ii in 1:ntheta
                        rownums = _split_numbers(ntap[ila_last + 8 + ii])
                        if length(rownums) >= 2
                            # pha[theta, layer, component]
                            theta_idx = ii
                            if theta_idx <= 112
                                pha[theta_idx, ilayer, ic] = rownums[2]
                            end
                        end
                    end

                    # bac is the last pha value at this layer/component
                    bac[ilayer, ic] = pha[min(ntheta, 112), ilayer, ic]

                    # update buffers
                    extbuf[ibuf] = ext[ilayer,ic]
                    scabuf[ibuf] = sca[ilayer,ic]
                    absbuf[ibuf] = abs1[ilayer,ic]
                    sisbuf[ibuf] = sis[ilayer,ic]
                    asybuf[ibuf] = asy[ilayer,ic]
                    brebuf[ibuf] = bre[ilayer,ic]
                    bimbuf[ibuf] = bim[ilayer,ic]
                    bacbuf[ibuf] = bac[ilayer,ic]
                end
            end
        end

        # -------- aggregate optical parameters per layer (optpar) --------
        iop = 0
        kop_local = 0
        for l in 1:nl
            summe = 0.0
            summa = 0.0
            summs = 0.0
            sumssa = 0.0
            sumasf = 0.0
            supf18 = 0.0
            if jnopar[10] == 1
                for it in 1:112
                    supf[it] = 0.0
                end
            end

            for jc in 1:njc[l]
                summe  += acmr[jc,l] * ext[l,jc]
                summa  += acmr[jc,l] * abs1[l,jc]
                summs  += acmr[jc,l] * sca[l,jc]
                sumssa += acmr[jc,l] * sis[l,jc] * ext[l,jc]
                sumasf += acmr[jc,l] * asy[l,jc] * sca[l,jc]
                supf18 += acmr[jc,l] * bac[l,jc]
                if jnopar[10] == 1
                    for it in 1:112
                        # Note: original R used pha[l, jc, it] though array was pha[theta,layer,comp]
                        # This block is disabled under current flags (jnopar[10]==0).
                        supf[it] += acmr[jc,l] * pha[it, l, jc]
                    end
                end
            end

            extn[l] = summe
            absn[l] = summa
            scan[l] = summs
            pf18n[l] = supf18
            if jnopar[10] == 1
                for it in 1:112
                    phafu[it,l] = supf[it]
                end
            end

            ssa[l] = sumssa / summe
            asf[l] = sumasf / summs

            # absolute parameters
            exta[l] = extn[l] * n[l]
            absa[l] = absn[l] * n[l]
            scaa[l] = scan[l] * n[l]
            pf18a[l] = pf18n[l] * n[l]
            if jnopar[10] == 1 && norm == 1
                for it in 1:112
                    phafu[it,l] = phafu[it,l] * n[l]
                end
            end

            if norm == 1
                extn[l] = exta[l]
                absn[l] = absa[l]
                scan[l] = scaa[l]
                pf18n[l] = pf18a[l]
            end

            if jnopar[10] == 1
                itp = 1
                for it in 1:112
                    if jnangle[it] == 1
                        # phaf(itp,l) = phafu(it,l)  # phaf not defined in provided R snippet
                        itp += 1
                    end
                end
            end

            # Build oparam list in order of active jnopar flags (1..13)
            iop = 0
            kop_local = 0
            if jnopar[1] == 1
                iop += 1; oparam[iop, l] = extn[l]
            end
            if jnopar[2] == 1
                iop += 1; oparam[iop, l] = scan[l]
            end
            if jnopar[3] == 1
                iop += 1; oparam[iop, l] = absn[l]
            end
            if jnopar[4] == 1
                iop += 1; oparam[iop, l] = ssa[l]
            end
            if jnopar[5] == 1
                iop += 1; oparam[iop, l] = asf[l]
            end
            if jnopar[9] == 1
                iop += 1
                if jnopar[6] == 1; kop_local += 1; end
                if jnopar[7] == 1; kop_local += 1; end
                if jnopar[8] == 1; kop_local += 1; end
                kop_local += iop
                # li.ratio = exta / pf18a
                oparam[kop_local, l] = exta[l] / pf18a[l]
            end
            # jnopar[11..13] are 0 in the provided vector; skipped.
        end

        # -------- optical thickness (jnopar[6..8]) --------
        if jnopar[6] == 1
            if nprog == 2
                # Not taken in this config; included for parity.
                for ilayer in 1:nl
                    if nltyp[ilayer] == 1
                        hm[ilayer,1] = parlay[ilayer,1]
                    elseif nltyp[ilayer] == 2
                        hm[ilayer,1] = parlay[ilayer,2] + (exp(-boundl[ilayer] / parlay[ilayer,2]) + exp(-boundu[ilayer] / parlay[ilayer,2]))
                    end
                end
                # set extfta/extstr from parlay if that path is used
            end

            hu = h1[nl, prnr]
            ho = hu + hfta[prnr]
            z = 8.0
            hftae = z * (exp(-hu / z) - exp(-ho / z))

            odepth = 0.0
            for ilayer in 1:nl
                odepth += exta[ilayer] * hm[ilayer, prnr]
            end
            odepth += extfta[ilamb] * hftae + extstr[ilamb] * hstra[prnr]
            odeptha = odepth / log(10)

            turbr = 0.008569 * alamb[ilamb]^(-4) * (1 + 0.0113 * alamb[ilamb]^(-2) + 0.00013 * alamb[ilamb]^(-4))
            turbf = (odepth + turbr) / turbr

            # align indices in oparam for outputs 6..8
            # Determine current iop based on earlier adds
            iop2 = 0
            if jnopar[1] == 1; iop2 += 1; end
            if jnopar[2] == 1; iop2 += 1; end
            if jnopar[3] == 1; iop2 += 1; end
            if jnopar[4] == 1; iop2 += 1; end
            if jnopar[5] == 1; iop2 += 1; end
            if jnopar[9] == 1; iop2 += 1; end

            kop_index = iop2 + 1
            if jnopar[6] == 1
                oparam[kop_index, 1] = odepth
                oparam[kop_index, 2] = 0.0
                kop_index += 1
                iop2 += 1
            end
            if jnopar[7] == 1
                oparam[kop_index, 1] = odeptha
                oparam[kop_index, 2] = 0.0
                kop_index += 1
                iop2 += 1
            end
            if jnopar[8] == 1
                oparam[kop_index, 1] = turbf
                oparam[kop_index, 2] = 0.0
                iop2 += 1
            end
        end

        # -------- out4: collect results into optdep --------
        for l in 1:nl
            if nih != 0
                rht[l] = ahum[ihum]
            end
            if nop > iop
                if jnopar[9] == 1 && jnopar[10] == 1
                    oparam[iop, l] = oparam[nop - 1, l]
                elseif jnopar[9] == 1 && jnopar[10] == 0
                    oparam[iop, l] = oparam[nop, l]
                end
            end
            if iop <= 10
                if l == 1
                    optdep[ilamb, 1] = alamb[ilamb]
                    # optical depth is expected at oparam[6,l] when jnopar[6]==1
                    optdep[ilamb, 2] = oparam[6, l]
                end
            end
        end
    end

    # convert wavelength μm -> nm
    optdep[:, 1] .*= 1000.0
    return optdep
end