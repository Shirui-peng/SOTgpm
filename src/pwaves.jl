# TODO:
# - Write empty file when data is missing?
# - Prevent pairing with neighboring event?

# initialize model for ray prediction
global taup = pyimport("obspy.taup")
global model = taup.TauPyModel(model="prem")

"""
    downloadpwaves(eqname, stations)
Download the P waveforms for the event catalog specified by `eqname` and the stations (and
channels) specified in `stations`. Currently, one full hour after the event is downloaded.
"""
function downloadseisdata(eqname, tsname, stations; src="IRIS", paircat=false)
  if paircat

    # load pair catalog
    pairs = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

    # build catalog of all unique events in pair catalog
    events1 = pairs[:,[1;3:5]]
    events2 = pairs[:,[2;6:8]]
    rename!(events1, :event1=>:time, :latitude1=>:latitude, :longitude1=>:longitude, :depth1=>:depth)
    rename!(events2, :event2=>:time, :latitude2=>:latitude, :longitude2=>:longitude, :depth2=>:depth)
    events = sort(unique(vcat(events1, events2)))

  else
    # load event catalog
    events = DataFrame(CSV.File(@sprintf("data/catalogs/%s_%s.csv", eqname, tsname)))
  end
  
  # loop over stations
  for (i, s) in enumerate(stations)

    # select source (station-dependent if given as vector)
    if typeof(src) == Vector{String}
      source = src[i]
    else
      source = src
    end
    
    # path to P-waveform directory
    path = seisdatadir(eqname[1:5], s)

    # create directory if needed
    mkpath(path)

    # loop over events
    for e in eachrow(events)

      # formatted date and time
      fmttime = Dates.format(e.time, "yyyy-mm-ddTHH:MM:SS.ss")

      # filename
      filename = @sprintf("%s/%s.h5", path, fmttime)

      # check if file exists already
      if !isfile(filename)

        @printf("Downloading %s ... ", fmttime)

        # download data
        try

          S = get_data("FDSN", s; s=string(e.time), t=string(e.time + Hour(1)), src=source, xf=tempname(cleanup=true))

          # check if there are no gaps
          if length(S) == 1 && size(S[1].t, 1) > 0
            
            # remove instrument response
            remove_resp!(S)
            
            # write to file
            h5open(filename, "w") do fid

              # station latitude
              fid["latitude"] = S[1].loc.lat

              # station longitude
              fid["longitude"] = S[1].loc.lon

              # start time in microseconds since 1970-01-01T00:00:00
              fid["starttime"] = S[1].t[1,2]

              # sampling frequency (in Hz)
              fid["fs"] = S[1].fs
              
              # trace
              fid["trace"] = S[1].x

            end

            @printf("done\n")
        
          else
        
            # create empty file, so download is not attempted again
            touch(filename)

            @printf("%d gap(s)\n", length(S)-1)
          
          end
          
        catch y
          if occursin("MAJO.00", s) || occursin("WAKE.00", s)
            @printf("change to a channel without 00\n")
            sloc,smid = "00",""
            #if occursin("WAKE.10", s)
            #  sloc,smid = "10","01"
            #end
            ss = split(s, sloc)
            snew = @sprintf("%s%s%s", ss[1], smid, ss[2])
            try
              
              S = get_data("FDSN", snew; s=string(e.time), t=string(e.time + Hour(1)), src="IRIS", xf=tempname(cleanup=true))
            
              # check if there are no gaps
              if length(S) == 1 && size(S[1].t, 1) > 0
              
                remove_resp!(S)
                
                # write to file
                h5open(filename, "w") do fid
    
                  # station latitude
                  fid["latitude"] = S[1].loc.lat
    
                  # station longitude
                  fid["longitude"] = S[1].loc.lon
    
                  # start time in microseconds since 1970-01-01T00:00:00
                  fid["starttime"] = S[1].t[1,2]
    
                  # sampling frequency (in Hz)
                  fid["fs"] = S[1].fs
                  
                  # trace
                  fid["trace"] = S[1].x
    
                end
    
                @printf("done\n")
                
              else
        
                # create empty file, so download is not attempted again
                touch(filename)
    
                @printf("%d gap(s)\n", length(S)-1)
          
              end
            catch y
              @printf("failed\n")
              # create empty file, so download is not attempted again
              if isa(y, LightXML.XMLParseError)
                touch(filename)
              end
            end
            
          else
            # create empty file, so download is not attempted again
            if isa(y, LightXML.XMLParseError)
              touch(filename)
            end
          
          end
        
        end
    
      end
  
    end
  
  end

end

"""
    cutpwaves(eqname, stations, intervals, freqbands)
Cut out and filter P-waveform sections used for cross-correlation. The data are read from
file based on `eqname` and `stations`, and are cut around the predicted P-wave arrival
time. The `intervals` specify what intervals around that predicted arrivals should be cut,
e.g. `intervals = [[-3, 47], [-3, 47]]` cuts a 50 s interval that starts 3 s before the
predicted arrival for both P-wave stations to be processed. The traces are filtered (prior
to cutting) with a bandpass specified in `freqbands`, e.g. `freqbands = [[1, 3], [1.5, 2.5]`
will filter data from the first station to 1 to 3 Hz and that from the second station to
1.5 to 2.5 Hz.
"""
function cutpwaves(eqname, tsname, stations, intervals, freqbands; paircat=false, getsnr=false)

  if paircat

    # load pair catalog
    pairs = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

    # build catalog of all unique events in pair catalog
    events1 = pairs[:,[1;3:5]]
    events2 = pairs[:,[2;6:8]]
    rename!(events1, :event1=>:time, :latitude1=>:latitude, :longitude1=>:longitude, :depth1=>:depth)
    rename!(events2, :event2=>:time, :latitude2=>:latitude, :longitude2=>:longitude, :depth2=>:depth)
    events = sort(unique(vcat(events1, events2)))

  else
    # load event catalog
    events = DataFrame(CSV.File(@sprintf("data/catalogs/%s_%s.csv", eqname,tsname)))
  end

  for i in 1 : length(stations)

    # path to P-waveform directory
    datapath = seisdatadir(eqname[1:5], stations[i])
    wavepath = pwavedir(eqname, stations[i], intervals[i], freqbands[i])

    # create directory if needed
    mkpath(wavepath)

    if getsnr
      # initialize snr catalog
      eventsnr = DataFrame(event=DateTime[], latitude=Float64[], longitude=Float64[], magnitude=Float64[],snr=Float64[])
    end
    
    for e in eachrow(events)

      # filename
      datafile = @sprintf("%s/%s.h5", datapath, fmttime(e.time))
      wavefile = @sprintf("%s/%s.h5", wavepath, fmttime(e.time))

      @printf("%s %s\n", stations[i], fmttime(e.time))

      # check if input file is present and output file is not
      if isfile(datafile) && !isfile(wavefile) && filesize(datafile) > 0

        # read data
        try 
          fid = h5open(datafile, "r")
          latitude = read(fid, "latitude")
          longitude = read(fid, "longitude")
          starttime = read(fid, "starttime")
          fs = read(fid, "fs")
          trace = read(fid, "trace")
          close(fid)

          # time stamps
          time = (1:length(trace))/fs

          # predict P-wave travel time
          d = dist(e.longitude, e.latitude, longitude, latitude)
          
          traveltime = model.get_travel_times(source_depth_in_km=e.depth,
                                              distance_in_degree=rad2deg(d/earthradius),
                                              phase_list=["P","p"])[1].time
            
          # predict S-wave travel time
          traveltime2= model.get_travel_times(source_depth_in_km=e.depth,
                                              distance_in_degree=rad2deg(d/earthradius),
                                              phase_list=["S","s"])[1].time      
          dt_PtoS = traveltime2 - traveltime
          windowlength = min(6.0*dt_PtoS,intervals[i][2]-intervals[i][1])
          timeafter = windowlength+intervals[i][1]

          # time shift due to discrepancy between start time of time series and event time

          timeshift = (starttime - 10^3*(e.time - DateTime(1970, 1, 1)).value)

          # find times
          idx = intervals[i][1] .< time .+ timeshift/1e6 .- traveltime .< timeafter

          # check if data is present (25 data points needed for bandpass filter)
          if sum(idx) > 25
        
            # band pass filter
            responsetype = Bandpass(freqbands[i][1], freqbands[i][2]; fs=fs)
            designmethod = Butterworth(4)
            filttrace = filtfilt(digitalfilter(responsetype, designmethod), trace)
            
            # cut
            cutstarttime = starttime + Int(round(10^6*time[idx][1]))
            cuttrace = filttrace[idx]
          
            if getsnr
              try  
                tmbf_s,tmaf_s,tmbf_n = 0., 30., 5.
                idx_s = -tmbf_s .< time .+ timeshift/1e6 .- traveltime .< tmaf_s
                idx_n = -tmbf_n .< time .+ timeshift/1e6 .- traveltime .< -tmbf_s
                snr = max(abs.(filttrace[idx_s])...)/max(abs.(filttrace[idx_n])...)
                # record in catalog
                push!(eventsnr, [e.time, e.latitude, e.longitude, e.magnitude,snr])
              catch y
                @printf("missing SNR!\n")
                continue
               end
            end
 
            # save to file
            h5open(wavefile, "w") do fid
              fid["latitude"] = latitude
              fid["longitude"] = longitude
              fid["starttime"] = cutstarttime
              fid["trace"] = cuttrace
              fid["fs"] = fs
            end
          end
          #end
        catch y
          try
            close(fid)
            @printf("missing starttime!\n")
          catch y
            @printf("missing file!\n")
          end
          # create empty file, so open is not attempted again
          touch(datafile)
          continue
        end

      end

    end
    
    if getsnr
      # save catalog to file
      CSV.write(@sprintf("data/catalogs/%s_%s_%3.1f_%3.1f_snr.csv",eqname[1:5], stations[i], freqbands[i][1],freqbands[i][2]), eventsnr)
    end
  end

end

"""
    findpairs(eqname, stations, intervals, freqbands; saveplot=false)
Cross-correlate P waveforms to find repeating pairs. The events are read from catalogs
according to `eqname` and `stations`, and the waveforms are read from folders based on
`eqname`, `stations`, `intervals`, `freqbands`. Set `saveplot=true` if plots of the
cross-correlation functions should be saved.
"""
function findpairs(eqname, tsname, stations, intervals, freqbands; saveplot=false, dphy = [50,50,1.5], jstart=1, jstop=0, kstart=0)

  # load event catalog
  allevents = DataFrame(CSV.File(@sprintf("data/catalogs/%s_%s.csv", eqname, tsname)))

  # loop over stations
  for i = 1 : length(stations)
    
    # path to P-wave directory
    wavepath = pwavedir(eqname[1:5], stations[i], intervals[i], freqbands[i])

    # path to plot directory
    plotpath = pplotdir(eqname, stations[i], intervals[i], freqbands[i])
    if saveplot
      mkpath(plotpath)
    end

    # read traces, sampling rates, and start times from file
    latitudes = Float64[]
    longitudes = Float64[]
    starttimes = Int[]
    fs = Float64[]
    traces = Array{Float64,1}[]
    present = falses(size(allevents, 1))
    for i = 1 : size(allevents, 1)
      filename = @sprintf("%s/%s.h5", wavepath, fmttime(allevents[i,:time]))
      if isfile(filename)
        present[i] = true
        push!(latitudes, h5read(filename, "latitude"))
        push!(longitudes, h5read(filename, "longitude"))
        push!(starttimes, h5read(filename, "starttime"))
        push!(fs, h5read(filename, "fs"))
        push!(traces, h5read(filename, "trace"))
      end
    end

    # delete events without file
    events = allevents[present,:]
    @printf("%s available events number: %d\n", stations[i], size(events, 1))

    # initialize pair catalog
    pairs = DataFrame(event1=DateTime[], latitude1=Float64[], longitude1=Float64[], depth1=Float64[], magnitude1=Float64[],
                      event2=DateTime[], latitude2=Float64[], longitude2=Float64[], depth2=Float64[], magnitude2=Float64[],
                      Δτ=Float64[], cc=Float64[])
    
    ijstop = jstop                  
    if jstop<1
      ijstop = size(events, 1) - 1
    else
      @printf("try j stop time: %s\n",allevents[ijstop,:time])
      ijstop = findfirst(x->x>allevents[ijstop,:time], events[:,:time]) - 1
      if isnothing(ijstop) || ijstop<1
        continue
      end
      @printf("final j stop time: %s\n",events[ijstop,:time])
    end
    ikstart = kstart
    if kstart>0
      @printf("try k start time: %s\n",allevents[ikstart,:time])
      ikstart = findlast(x->x<allevents[ikstart,:time], events[:,:time]) + 1
      if isnothing(ikstart) || ikstart>size(events, 1)
        continue
      end
      @printf("final k start time: %s\n",events[ikstart,:time])
    end
    @printf("jstart = %d, jstop = %d, kstart = %d\n", jstart,ijstop, ikstart)

    # loop over event pairs
    for j = jstart : ijstop, k = max(ikstart,j + 1) : size(events, 1)
    
      dhypocenter = dist(events[j,:longitude], events[j,:latitude], events[k,:longitude], events[k,:latitude])/1e3
      ddepth = abs(events[j,:depth]-events[k,:depth])
      dM = abs(events[k,:magnitude]-events[j,:magnitude])
    
      # ensure sampling rate is identical
      if fs[j] == fs[k] && dhypocenter<dphy[1] && ddepth<dphy[2] && dM<dphy[3]

        # cross-correlation measurement
        if saveplot
          maxcc, Δτ, cc, lags = xcorr(traces[j], traces[k], fs[j])
        else
          maxcc, Δτ = xcorr(traces[j], traces[k], fs[j])
        end
        
        # record if CC is above 0.9, sampling rate is identical, and pairs is not a repeater
        if maxcc ≥ 0.9 && starttimes[k] - starttimes[j] > length(traces[j])/fs[j]*1e6
        
          # adjust for starttime recorded in waveforms
          originadjustment = starttimes[k] - starttimes[j] - 10^3*(events[k,:time] - events[j,:time]).value
          Δτ += originadjustment/1e6

          @printf("%s %s %s %5.3f %+6.3f\n", stations[i], fmttime(events[j,:time]), fmttime(events[k,:time]), maxcc, Δτ)

          # record in catalog
          push!(pairs, [events[j,:time], events[j,:latitude], events[j,:longitude], events[j,:depth], events[j,:magnitude],
                        events[k,:time], events[k,:latitude], events[k,:longitude], events[k,:depth], events[k,:magnitude],
                        Δτ, maxcc])

          # save catalog to file
          if jstart==1 && jstop == 0 && kstart == 0
            CSV.write(paircatfileold(eqname, tsname, stations[i], intervals[i], freqbands[i]), pairs)
          else
            CSV.write(paircatfile(eqname, tsname, stations[i], intervals[i], freqbands[i], jstart, jstop, kstart), pairs)
          end
          
          # plot cross-correlation function if desired
          if saveplot
            # lengths of waveforms
            n1 = length(traces[j])
            n2 = length(traces[k])
            fig = figure()
            plot(lags/fs[j] .+ originadjustment/1e6, circshift(cc, (n1+n2)÷2))
            xlim(Δτ-5, Δτ+5)
            xlabel("lag (s)")
            ylabel("cross-correlation")
            savefig(@sprintf("%s/%s_%s.pdf", plotpath, fmttime(events[j,:time]),
                             fmttime(events[k,:time])))
            close(fig)

            time1plot = (0:n1-1)/fs[j]
            time2plot = (0:n2-1)/fs[k].+ originadjustment/1e6
            fig, axs = subplots(2, 1)
            ax = axs[1]
            ax.plot(time1plot, traces[j], color="black", alpha=0.75, label=string("event1,M",events[j,:magnitude]))
            ax.plot(time2plot, traces[k], color="red", alpha=0.75, label=string("event2,M",events[k,:magnitude]))
            ax.legend(frameon=false, loc=4)
            ax.set_title(@sprintf("maxCC: %4.2f, offset:%+6.3f\n", maxcc, Δτ))
            ax.set_ylabel("absolute amplitude")
            ax = axs[2]
            ax.plot(time1plot, traces[j]/max(abs.(traces[j])...), color="black", alpha=0.75, label=string("event1,M",events[j,:magnitude]))
            ax.plot(time2plot.-Δτ, traces[k]/max(abs.(traces[k])...), color="red", alpha=0.75, label=string("shifted event2,M",events[k,:magnitude]))
            ax.legend(frameon=false, loc=4)
            ax.set_xlabel("time (s)")
            ax.set_ylabel("normalized amplitude")
            fig.savefig(@sprintf("%s/%s_%swaveform.pdf", plotpath, fmttime(events[j,:time]),
                             fmttime(events[k,:time])))
            close(fig)
          end

        end

      end
  
    end 
    @printf("%s findpair finished: %s\n", stations[i],Dates.format(now(), "HH:MM"))
    # save catalog to file
    #CSV.write(paircatfile(eqname, tsname, stations[i], intervals[i], freqbands[i]), pairs)

  end

end

function findpairsold(eqname, tsname, stations, intervals, freqbands; saveplot=false, dphy = [50,50,1.5],)

  # load event catalog
  allevents = DataFrame(CSV.File(@sprintf("data/catalogs/%s_%s.csv", eqname, tsname)))

  # loop over stations
  for i = 1 : length(stations)

    # path to P-wave directory
    wavepath = pwavedir(eqname, stations[i], intervals[i], freqbands[i])

    # path to plot directory
    plotpath = pplotdir(eqname, stations[i], intervals[i], freqbands[i])
    if saveplot
      mkpath(plotpath)
    end

    # read traces, sampling rates, and start times from file
    latitudes = Float64[]
    longitudes = Float64[]
    starttimes = Int[]
    fs = Float64[]
    traces = Array{Float64,1}[]
    present = falses(size(allevents, 1))
    for i = 1 : size(allevents, 1)
      filename = @sprintf("%s/%s.h5", wavepath, fmttime(allevents[i,:time]))
      if isfile(filename)
        present[i] = true
        push!(latitudes, h5read(filename, "latitude"))
        push!(longitudes, h5read(filename, "longitude"))
        push!(starttimes, h5read(filename, "starttime"))
        push!(fs, h5read(filename, "fs"))
        push!(traces, h5read(filename, "trace"))
      end
    end

    # delete events without file
    events = allevents[present,:]

    # initialize pair catalog
    pairs = DataFrame(event1=DateTime[], latitude1=Float64[], longitude1=Float64[],
                      event2=DateTime[], latitude2=Float64[], longitude2=Float64[],
                      Δτ=Float64[], cc=Float64[])

    # loop over event pairs
    for j = 1 : size(events, 1) - 1, k = j + 1 : size(events, 1)

      # lengths of waveforms
      n1 = length(traces[j])
      n2 = length(traces[k])

      # pad with zeros
      padtrace1 = [traces[j]; zeros(n2)]
      padtrace2 = [traces[k]; zeros(n1)]

      # calculate cross-correlation function
      norm = sqrt(sum(traces[j].^2)*sum(traces[k].^2))
      cc = irfft(conj.(rfft(padtrace1)).*rfft(padtrace2), n1+n2)/norm

      # find index of max cross-correlation
      maxidx = argmax(cc)

      # calculate max CC using quadratic interpolation
      maxcc = maxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)])

      # record if CC is above 0.9, sampling rate is identical, and pairs is not a repeater
      if maxcc ≥ 0.9 && fs[j] == fs[k] && starttimes[k] - starttimes[j] > n1/fs[j]*1e6

        # integer lags (depending on whether n1 + n2 is even or odd)
        lags = iseven(n1+n2) ? (-(n1+n2)÷2 : (n1+n2)÷2-1) : (-(n1+n2)÷2 : (n1+n2)÷2)

        # calculate lag from grid max CC
        Δτ = mod(maxidx-1, lags)/fs[j]

        # adjust using quadratic interpolation
        Δτ += argmaxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)], 1/fs[j])

        # adjust for starttime recorded in waveforms
        originadjustment = starttimes[k] - starttimes[j] - 10^3*(events[k,:time]
                                                                 - events[j,:time]).value
        Δτ += originadjustment/1e6

        @printf("%s %s %s %4.2f %+6.3f\n", stations[i], fmttime(events[j,:time]),
                fmttime(events[k,:time]), maxcc, Δτ)

        # record in catalog
        push!(pairs, [events[j,:time], events[j,:latitude], events[j,:longitude],
                      events[k,:time], events[k,:latitude], events[k,:longitude],
                      Δτ, maxcc])

        # plot cross-correlation function if desired
        if saveplot
          fig = figure()
          plot(lags/fs[j] .+ originadjustment/1e6, circshift(cc, (n1+n2)÷2))
          xlim(Δτ-5, Δτ+5)
          xlabel("lag (s)")
          ylabel("cross-correlation")
          savefig(@sprintf("%s/%s_%s.pdf", plotpath, fmttime(events[j,:time]),
                           fmttime(events[k,:time])))
          close(fig)
          
          time1plot = (0:n1-1)/fs[j]
          time2plot = (0:n2-1)/fs[k] .+ originadjustment/1e6
          fig, axs = subplots(2, 1)
          ax = axs[1]
          ax.plot(time1plot, traces[j], color="black", alpha=0.75, label=string("event1,M",events[j,:magnitude]))
          ax.plot(time2plot, traces[k], color="red", alpha=0.75, label=string("event2,M",events[k,:magnitude]))
          ax.legend(frameon=false, loc=4)
          ax.set_title(@sprintf("maxCC: %4.2f, offset:%+6.3f\n", maxcc, Δτ))
          ax.set_ylabel("absolute amplitude")
          ax = axs[2]
          ax.plot(time1plot, traces[j]/max(abs.(traces[j])...), color="black", alpha=0.75, label=string("event1,M",events[j,:magnitude]))
          ax.plot(time2plot.-Δτ, traces[k]/max(abs.(traces[k])...), color="red", alpha=0.75, label=string("shifted event2,M",events[k,:magnitude]))
          ax.legend(frameon=false, loc=4)
          ax.set_xlabel("time (s)")
          ax.set_ylabel("normalized amplitude")
          fig.savefig(@sprintf("%s/%s_%swaveform.pdf", plotpath, fmttime(events[j,:time]),fmttime(events[k,:time])))
          close(fig)
        end

      end

    end
    
    # save catalog to file
    #CSV.write(paircatfile(eqname, stations[i], intervals[i], freqbands[i]), pairs)

  end

end

"""
    measurepairs(eqname, stations, intervals, freqbands; saveplot=false)
Cross-correlate P waveforms for cataloged pairs. The event pairs are read from catalogs
according to `eqname` and `stations`, and the waveforms are read from folders based on
`eqname`, `stations`, `intervals`, `freqbands`. Set `saveplot=true` if plots of the
cross-correlation functions should be saved.
"""
function measurepairs(eqname, stations, intervals, freqbands)

  # load event pair catalog
  allpairs = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

  # loop over stations
  for i = 1 : length(stations)

    # path to P-wave directory
    wavepath = pwavedir(eqname, stations[i], intervals[i], freqbands[i])

    # initialize pair catalog
    pairs = DataFrame(event1=DateTime[], latitude1=Float64[], longitude1=Float64[], event2=DateTime[], latitude2=Float64[], longitude2=Float64[], Δτ=Float64[], cc=Float64[])

    # loop over event pairs
    for pair in eachrow(allpairs)

      # data file names
      filename1 = @sprintf("%s/%s.h5", wavepath, fmttime(pair.event1))
      filename2 = @sprintf("%s/%s.h5", wavepath, fmttime(pair.event2))

      @printf("%s %s %s ", stations[i], fmttime(pair.event1), fmttime(pair.event2))

      # ensure files are present
      if !isfile(filename1)
        println("file1 missing")
      elseif !isfile(filename2)
        println("file2 missing")
      else

        # read data from files
        latitude1 = h5read(filename1, "latitude")
        latitude2 = h5read(filename2, "latitude")
        longitude1 = h5read(filename1, "longitude")
        longitude2 = h5read(filename2, "longitude")
        starttime1 = h5read(filename1, "starttime")
        starttime2 = h5read(filename2, "starttime")
        fs1 = h5read(filename1, "fs")
        fs2 = h5read(filename2, "fs")
        trace1 = h5read(filename1, "trace")
        trace2 = h5read(filename2, "trace")

        # ensure sampling rate is identical
        if fs1 != fs2
          println("fs1 ≠ fs2")
        else

          # cross-correlation measurement
          maxcc, Δτ = xcorr(trace1, trace2, fs1)

          # record if CC is above 0.9 and pair is not a self-repeater
          if maxcc < 0.9
            @printf("%5.3f\n", maxcc)
          elseif starttime2 - starttime1 ≤ length(trace1)/fs1*1e6
            println("suspected self-repeater")
          else

            # adjust for starttime recorded in waveforms
            originadjustment = starttime2 - starttime1 - 10^3*(pair.event2 - pair.event1).value
            Δτ += originadjustment/1e6

            @printf("%5.3f %+6.3f\n", maxcc, Δτ)

            # record in catalog
            push!(pairs, [pair.event1, pair.latitude1, pair.longitude1, pair.event2, pair.latitude2, pair.longitude2, Δτ, maxcc])

          end

        end

      end

    end
    
    # save catalog to file
    CSV.write(paircatfile(eqname, stations[i], intervals[i], freqbands[i]), pairs)

  end

end

"Cross-correlation measurement"
function xcorr(trace1, trace2, fs; maxlag=12)

  # lengths of waveforms
  n1 = length(trace1)
  n2 = length(trace2)

  # pad with zeros
  padtrace1 = [trace1; zeros(n2)]
  padtrace2 = [trace2; zeros(n1)]

  # calculate covariance function
  cov = irfft(conj.(rfft(padtrace1)).*rfft(padtrace2), n1+n2)

  # indicator functions
  ind1 = [ones(n1); zeros(n2)]
  ind2 = [ones(n2); zeros(n1)]

  # calculate normalization factors
  norm1 = irfft(conj.(rfft(padtrace1.^2)).*rfft(ind2), n1+n2)
  norm2 = irfft(conj.(rfft(ind1)).*rfft(padtrace2.^2), n1+n2)

  # exclude negative and zero normalization factors (due to roundoff error and non-overlapping segments)
  norm1[norm1.≤0] .= Inf
  norm2[norm2.≤0] .= Inf

  # cross-correlation
  cc = cov./sqrt.(norm1.*norm2)

  # find index of max covariance in allowable range (±maxlag, in seconds)
  k = Int(round(maxlag*fs))
  cca = copy(cc); cca[k+2:n1+n2-k] .= 0
  maxidx = argmax(cca)
  
  pks, vals = findmaxima(cca)
  
  # check whether CC = 0
  if cc[maxidx] == 0 || length(findall(x->x>0.9,vals))>3

    return 0.0, NaN, 0.0.*cca, 0.0.*cca
 
  else

    # integer lags (depending on whether n1 + n2 is even or odd)
    lags = iseven(n1+n2) ? (-(n1+n2)÷2 : (n1+n2)÷2-1) : (-(n1+n2)÷2 : (n1+n2)÷2)
    
    # index shift from grid max covariance
    Δi = mod(maxidx-1, lags)

    # calculate max CC using quadratic interpolation
    maxcc = maxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)])

    # time shift adjusted using quadratic interpolation
    Δτ = Δi/fs + argmaxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)], 1/fs)

    return maxcc, Δτ, cca, lags

  end

end

"Time format"
fmttime(time) = Dates.format(time, "yyyy-mm-ddTHH:MM:SS.ss")

"Directory for raw P-wave data"
seisdatadir(eqname, station) = @sprintf("data/seisdata/%s_%s", eqname, station)

"Directory for processed P waveforms"
pwavedir(eqname, station, interval, freqband) = @sprintf("data/pwaves/%s_%s_%+03d_%+03d_%3.1f_%3.1f", eqname, station, interval[1], interval[2], freqband[1], freqband[2])

"Directory for P-wave cross-correlation plots"
pplotdir(eqname, station, interval, freqband) = @sprintf("data/pplots/%s_%s_%+03d_%+03d_%3.1f_%3.1f", eqname, station, interval[1], interval[2], freqband[1], freqband[2])

"File name for P-wave pair catalogs"
paircatfileold(eqname, tsname, station, interval, freqband) = @sprintf("data/catalogs/%s_%s_%s_%+03d_%+03d_%3.1f_%3.1f.csv", eqname, tsname, station, interval[1], interval[2], freqband[1], freqband[2])
paircatfile(eqname, tsname, station, interval, freqband, jstart, jstop, kstart) = @sprintf("data/catalogs/%s_%s_%s_%+03d_%+03d_%3.1f_%3.1f_j%dto%d_k%d.csv", eqname, tsname, station, interval[1], interval[2], freqband[1], freqband[2], jstart, jstop,kstart)
