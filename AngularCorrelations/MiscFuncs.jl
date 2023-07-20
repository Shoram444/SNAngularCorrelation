module MiscFuncs
using DataFrames, DataFramesMeta, FHist, StatsBase, Polynomials

export  y,                   # Returns y(x) = x 
        get_theta_esc_slice,    # Returns a slice of ϕ within given data cuts
        Gaus,                   # Returns the value of normal distribution
        get_chi2,               # Returns chi2 of the fit
        get_reco_efficiency,    # DEPRECATED! Returns the efficiency of the reconstructed slice
        get_diagonal_sums,      # Returns a vector of the sums of the diagonals of the 2d histogram.
        get_k_factors,          # Returns a vector of the k-factors, representing by how much there is over/under-estimate.
        get_diagonal_sums_abs,  # DEPRECATED! Returns a vector of the absolute sums of the diagonals of the 2d histogram.
        get_k_factors_abs,      # DEPRECATED! Returns a vector of the absolute k-factors, representing by how much there is deviation.
        get_slice_stats,        # Returns a vector of stats for the given slice of ϕ 
        get_cut_edges,          # Returns a pair of cut edges for dϕ
        get_gs,                 # DEPRECATED! Returns g(k, θ) 
        get_gs_rel,             # DEPRECATED! Returns g(k, θ) relative
        get_diagonal_sums_cdf,  # DEPRECATED! Returns a vector of the cdf sums of the diagonals of the 2d histogram.
        shift_angle,            # Returns shifted angle ϕ' = ϕ + s
        get_gs_mean,            # Returns mean value of the g(k) vector as if from histogram!
        get_rms,                # Returns rms value of the g(k) vector as if from histogram or from a vector!
        turn90_cc,              # Returns nxm matrix with elements turned by 90° counter-clockwise (original matrix would be mxn) 
        get_gs_df,              # Returns a DataFrame from tree with columns: k, dϕ₁, dϕ₂, ... 
        df_to_mat,              # Returns a Matrix from dfG for plotting in heatmaps, surfaces, contours, etc. (the matrix is turned so the x-axis corresponds to k, y-axis to dϕ)
        get_bin_center,         # Returns bin center of angle in range 0-180
        get_rms_set,            # Returns a set of RMS values for each dϕ slice (each row represents one dϕ slice, each value in row is rms with given shift)
        get_min_rms,            # Returns a vector of minized RMS values (minimum from the set of RMS values obtained by get_rms_set)
        get_min_shifts,         # Returns a vector of minimized shifts 
        get_phi_prime,          # Returns the value of ϕ' shifted by the provided shift vector
        get_hist_mean,          # Returns the mean value of the histogram
        get_hist_median         # Returns the !!binCenter!! of the bin containing the median value (50% statistic)

y(x) = x # generic linear line y = x 

function get_theta_esc_slice(_θemiMin::Real, _θemiMax::Real, _tree::DataFrame)
    return θesc = _tree[
        (_tree.thetaEmitted .> _θemiMin) .& (_tree.thetaEmitted .< _θemiMax),
        :thetaEscaped,
    ]
end
function get_theta_esc_slice(
    _θemiMin::Real,
    _θemiMax::Real,
    _EMin::Real,
    _EMax::Real,
    _tree::DataFrame,
)
    return θesc = _tree[
        (_tree.thetaEmitted .> _θemiMin) .& (_tree.thetaEmitted .< _θemiMax) .& (_tree.ESum .> _EMin) .& (_tree.ESum .< _EMax),
        :thetaEscaped,
    ] # return slice where all conditions are fulfilled
end

function Gaus(fit, x, N)
    return N / (fit.σ * sqrt(2 * π)) * exp(-0.5 * ((x - fit.μ) / fit.σ)^2)
end

function Gaus(μ, σ, x, N)
    return N / (σ * sqrt(2 * π)) * exp(-0.5 * ((x - μ) / σ)^2)
end

function get_chi2(x::Vector{<:Real}, y::Vector{<:Real}, f::Function)
    xi = 0
    ndf = 0
    for idx in eachindex(y)
        if (y[idx] <= 0)
            continue
        else
            xi += (y[idx] - f(x[idx]))^2 / y[idx]
            ndf += 1
        end
    end
    return xi, ndf - 3 # there are 3 fit parameters (μ. σ, N) 
end

function get_chi2(y::Vector{<:Real}, fx::Vector{<:Real})
    xi = 0
    ndf = 0
    for idx in eachindex(y)
        if (y[idx] <= 0)
            continue
        else
            xi += (y[idx] - fx[idx])^2 / y[idx]
            ndf += 1
        end
    end
    return xi, ndf - 3 # there are 3 fit parameters (μ. σ, N) 
end

"""
get_pmSigma_mode(_hist; _sigmaHalf=0.3414)

## Description of get_pmSigma_mode

Returns how many bins to the left and to the right from mode it is to cover ±σ (σ = 0.3414, or as stated in argument).
The returned values are [l,r] vector where l represents mode-sigmaBin, and r represents sigmaBin-mode.

Input arguments are:

  - Hist1D objed from FHist.jl
  - Sigma portion (by default 0.3414)
"""
function get_pmSigma_mode(_hist; _sigmaHalf = 0.3414)
    totalEvents = sum(bincounts(_hist))
    mode_bin = 0 # in which bin the mode is placed

    if (length(argmax(bincounts(_hist))) == 1 || sum(bincounts(_hist)) == 0) # makes sure there is a mode to begin with, if not return 0
        mode_bin = argmax(bincounts(_hist))
    else
        return [0, 0]
    end

    tempSumLeft = (bincounts(_hist)[mode_bin] / 2) / totalEvents   # to begin with the sum is initialized to 1/2 of mode's bin
    tempSumRight = (bincounts(_hist)[mode_bin] / 2) / totalEvents

    binLeftCounter = mode_bin - 1
    binRightCounter = mode_bin + 1

    while (tempSumLeft <= _sigmaHalf && binLeftCounter != 0) #68.27% / 2 half of one sigma
        tempSumLeft += bincounts(_hist)[binLeftCounter] / total # add the content of the bin to the left
        binLeftCounter -= 1
    end

    while (tempSumRight <= _sigmaHalf && binRightCounter != length(bincounts(_hist))) #68.27% / 2 half of one sigma
        tempSumRight += bincounts(_hist)[binRightCounter] / total # add the content of the bin to the right
        binRightCounter += 1
    end

    return [mode_bin - binLeftCounter, binRightCounter - mode_bin]
end

"""
get_reco_efficiency(_θemiMin::Real, _θemiMax::Real,_slice::Vector{Real})

## Description of get_reco_efficiency

Returns the efficiency of the reconstructed slice as:

`ε = (n of events escaped within (θ_{min}, θ_{max})/(n of all events in (θ_{min}, θ_{max})`

Input arguments are:

  - Hist1D objed from FHist.jl
  - Sigma portion (by default 0.3414)
"""
function get_reco_efficiency(_θemiMin::Real, _θemiMax::Real, _slice::Vector{<:Real})
    n_passed = filter(θ -> (_θemiMin <= θ <= _θemiMax), _slice) |> length
    return length(_slice) > 0 ? n_passed / length(_slice) : 0
end

function get_reco_efficiency(
    _θemiMin::Real,
    _θemiMax::Real,
    _EMin::Real,
    _EMax::Real,
    _tree::DataFrame,
)
    return get_reco_efficiency(
        _θemiMin,
        _θemiMax,
        get_theta_esc_slice(_θemiMin, _θemiMax, _EMin, _EMax, _tree),
    )
end

"""
get_diagonal_sums(h2d)

## Description of get_diagonal_sums

Returns a vector of the sums of the diagonals of the 2d histogram.
(Diagonals from top-left to bottom-right)

Input arguments are:

  - Hist2D object from FHist.jl
"""
function get_diagonal_sums(h2d)
    mat = bincounts(h2d)      # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]        # get the number of rows/cols of mxm matrix
    ndiagonals = 2 * (m) - 1       # there are 2m-1 diagonals in a mxm matrix
    sums = zeros(ndiagonals)   # initialize the array of sums of the diagonal elements

    for r in 1:m           # loop over rows r
        for c in 1:m       # loop over columns c
            sums[m - (c - r)] += mat[r, c]
        end
    end
    return sums
end

function get_diagonal_sums_abs(h2d)
    mat = bincounts(h2d)       # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]         # get the number of rows/cols of mxm matrix
    sums = zeros(m + 1)           # initialize the array of sums of the diagonal elements
    # there are m+1 diagonals in a mxm matrix

    for i in 1:m                # loop over rows i
        for j in 1:m            # loop over columns j
            sums[Int(abs(i - j) + 1)] += mat[i, j]
        end
    end
    return sums
end

function get_diagonal_sums_cdf(h2d)
    mat = bincounts(h2d)      # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]        # get the number of rows/cols of mxm matrix
    ndiagonals = 2 * (m) - 1       # there are 2m-1 diagonals in a mxm matrix
    sums = zeros(ndiagonals)   # initialize the array of sums of the diagonal elements

    for r in 1:m           # loop over rows r
        for c in 1:m       # loop over columns c
            sums[m - (c - r)] += mat[r, c]
        end
    end

    integral = sum(sums) # get the total sum of the array
    sums = cumsum(sums) ./ integral # get cumulative sum and normalize to the integral

    return sums
end

"""
get_k_factors(h2d)

## Description of get_k_factors

Returns a vector of the k-factors, representing by how much there is over/under-estimate.
(Diagonals from top-left to bottom-right)

Input arguments are:

  - Hist2D object from FHist.jl
"""
function get_k_factors(h2d)
    mat = bincounts(h2d)        # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]            # get the number of rows/cols of mxm matrix
    ks = zeros(2m - 1)     # initialize array of k-factors

    for d in 1:(2m - 1)            # loop over the diagnoal index in sums
        k = m - d 
        ks[d] = k
    end
    return ks
end

function get_k_factors_abs(h2d)
    mat = bincounts(h2d)      # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]        # get the number of rows/cols of mxm matrix

    return [k for k in 0:m]
end

function get_slice_stats(
    _θemiMin::Real,
    _θemiMax::Real,
    _EMin::Real,
    _EMax::Real,
    _theta_esc::Vector{<:Real},
    _binWidth = 1.0,
    _weights::AbstractWeights = Weights(ones(eltype(_theta_esc),length(_theta_esc))) # if no weights are specified, use 1.0
)
    meanData = 0
    modeData = 0
    medianData = 0
    varianceData = 0
    escEffData = 0

    hist = Hist1D(_theta_esc, _weights, (0:_binWidth:180)) # used to get mode from the binned histogram

    if (length(_theta_esc) > 0)
        meanData = get_hist_mean(hist)  
        medianData = get_hist_median(hist)
        varianceData = sqrt(var(_theta_esc; corrected = :true))
        escEffData = get_reco_efficiency(_θemiMin, _θemiMax, _theta_esc)
        modeData = argmax(bincounts(hist)) * _binWidth
    end

    return [
        _θemiMin,
        _θemiMax,
        _EMin,
        _EMax,
        meanData,
        modeData,
        medianData,
        varianceData,
        escEffData,
    ]
end

function get_slice_stats(
    _θemiMin::Real,
    _θemiMax::Real,
    _EMin::Real,
    _EMax::Real,
    _tree::DataFrame,
    _binWidth = 1.0,
)
    meanData = 0
    modeData = 0
    medianData = 0
    varianceData = 0
    escEffData = 0

    theta_esc = get_theta_esc_slice(_θemiMin, _θemiMax, _EMin, _EMax, _tree)
    hist = Hist1D(theta_esc, (0:_binWidth:180)) # used to get mode from the binned histogram

    if (length(theta_esc) > 0)
        meanData = mean(theta_esc)
        medianData = median(theta_esc)
        varianceData = sqrt(var(theta_esc; corrected = :true))
        escEffData = get_reco_efficiency(_θemiMin, _θemiMax, theta_esc)
        modeData = argmax(bincounts(hist)) * _binWidth
    end

    return [
        _θemiMin,
        _θemiMax,
        _EMin,
        _EMax,
        meanData,
        modeData,
        medianData,
        varianceData,
        escEffData,
    ]
end

function get_cut_edges(midPoint, multiple, dθ, sign)
    l = midPoint - multiple * dθ
    u = midPoint + multiple * dθ
    if (l < 0)
        l = 0
    end
    if (u > 180)
        u = 180
    end

    if (sign == "pm")
        return (l, u)
    elseif (sign == "p")
        return (midPoint, u)
    elseif (sign == "m")
        return (l, midPoint)
    else
        error("Must chose a valid sign from options: pm, m, p")
        return 0
    end
end

function get_gs(θ, k, gs) # return the value of g corresponding to the θ and k
    r = Int(floor(θ))+1
    c = Int(size(gs)[1] - floor(k)) 
    return gs[r][c]
end

function get_gs_rel(θ, k, gs) # return the value of g corresponding to the θ and k
    r = Int(floor(θ))+1
    c = (k+length(gs[1])) - 179
    return gs[r][c] ./ maximum(gs[r])
end

function shift_angle(angle, shifter::Vector{Float64})
    dϕ = 180/length(shifter)
    shiftIndex = Int(floor(angle/dϕ) +1 ) < 180/dϕ ? Int(floor(angle/dϕ) +1 ) : Int(180/dϕ) 
    newAngle = angle - shifter[shiftIndex]
end

function shift_angle(angle, shifter::Polynomial)
    newAngle = angle - shifter(angle)
end

function get_gs_mean(gsVector, ksVector)
    N = sum(gsVector)
    ss = 0.0
    for i in eachindex(gsVector)
        ss += gsVector[i] * ksVector[i] 
    end
    return ss/N
end

function get_rms( gsVector, ksVector )
    if (sum(gsVector) == 0.0 )
        # println("WARNING! Sum of g is 0! Provided empty diagonal. RMS is set to 0!")
        return 0.0
    end
    meanOfSquares = sum( gsVector .* (ksVector ) .^2 ) / sum(gsVector)  # sum[ (gᵢ*( kᵢ)^2 ) ] / sum(gᵢ)
    return sqrt( meanOfSquares )
end

function get_rms( gsVector )
    if (sum(gsVector) == 0.0 )
        # println("WARNING! Sum of g is 0! Provided empty diagonal. RMS is set to 0!")
        return 0.0
    end
    
    return sqrt(1/(length(gsVector))*sum(gsVector .^2))
end


function turn90_cc(_m)    # returns a matrix which elements are turned by 90degrees in counter-clowise direction
    mat = zeros(size(_m)) # initiaite matrix of zeros of same size as _m
    n   = size(_m)[1]     # number of rows of matrix _m
    c   = Vector{Float64}(undef,n) # placeholder for the elemnts of the column

    for i in 0:n-1
        c           = _m[:,end-i]      # take last column first, move to the "left"
        mat[i+1,:]  = reshape(c, 1, n) # last column becomes first row, second last column becomes second row, etc.
    end
    return mat
end


function get_gs_df(_tree, _dϕ, _sign = "p"; _theta = "thetaEmitted", _phi = "thetaEscaped")
    if("weights" in names(_tree))
        fh2d = Hist2D(                                          
            (_tree[!, _theta], _tree[!, _phi]),  
            Weights(_tree.weights),    
            (0:_dϕ:180, 0:_dϕ:180), 
        ) 
    else
        fh2d = Hist2D(                                          
            (_tree[!, _theta], _tree[!, _phi]),  
            (0:_dϕ:180, 0:_dϕ:180), 
        ) 
    end

    df = DataFrame(k = get_k_factors(fh2d) .* _dϕ)  # Gs saved in a dataframe first column is k-factor, next columns correspond to ϕ-slices
    
    for (i, n) in enumerate(1:_dϕ:180)
        cutEdges1 = get_cut_edges(n - 1, 1, _dϕ, _sign)                   # provides the lower and upper cut 
        if("weights" in names(_tree))
            
            sdf = @chain _tree begin                                         # filter out the dataframe
                @subset((cutEdges1[1] .<= _tree[!,_phi] .<= cutEdges1[2]))
                select(_phi, _theta, :weights)
            end
                
            fh2d = Hist2D(
                (sdf[!, _theta], sdf[!, _phi]),
                Weights(sdf.weights),
                (0:_dϕ:180, 0:_dϕ:180),
                )
        else
            sdf = @chain _tree begin                                         # filter out the dataframe
                @subset((cutEdges1[1] .<= _tree[!,_phi] .<= cutEdges1[2]))
                select(_phi, _theta)
            end

            fh2d = Hist2D(
                (sdf[!, _theta], sdf[!, _phi]),
                (0:_dϕ:180, 0:_dϕ:180),
                )
        end
        
    
        colName = string(cutEdges1[1],"-",cutEdges1[2])   # df column name in format minϕ-maxϕ
        df[!,colName] = get_diagonal_sums(fh2d)
    end

    return df
end

function get_gs_df(_tree, _dϕ, _sign = "p")
    fh2d = Hist2D(                                          
    (_tree.thetaEmitted, _tree.thetaEscaped),  
    Weights(_tree.weights),    
    (0:_dϕ:180, 0:_dϕ:180), 
    ) 

    df = DataFrame(k = get_k_factors(fh2d) .* _dϕ)  # Gs saved in a dataframe first column is k-factor, next columns correspond to ϕ-slices
    
    for (i, n) in enumerate(1:_dϕ:180)
        cutEdges1 = get_cut_edges(n - 1, 1, _dϕ, _sign)                   # provides the lower and upper cut 
    
        sdf = @chain _tree begin                                         # filter out the dataframe
            @subset((cutEdges1[1] .<= :thetaEscaped .<= cutEdges1[2]))
            @select(:thetaEscaped, :thetaEmitted, :weights)
        end

        fh2d = Hist2D(
            (sdf[!, :thetaEmitted], sdf[!, :thetaEscaped]),
            Weights(sdf.weights),
            (0:_dϕ:180, 0:_dϕ:180),
        )
    
        colName = string(cutEdges1[1],"-",cutEdges1[2])   # df column name in format minϕ-maxϕ
        df[!,colName] = get_diagonal_sums(fh2d)
    end

    return df
end

function df_to_mat( df::DataFrame )
    mat = Matrix( df[:,2:end] ) # create a matrix from df and drop first column

    nCol = size(mat)[2]   # m
    nRow = size(mat)[1]   # n

    matNew = zeros(nCol,nRow) # create zeros mxn matrix

    for rNew in 1:nCol
        matNew[rNew, :] = reverse(reshape(mat[:,rNew],1,nRow))
    end
    return matNew
end

function get_bin_center(x, nBins)
    dϕ = 180/nBins
    if x == 0.0
        return dϕ/2
    end
    return ceil(x/dϕ)*dϕ - dϕ/2
end

function get_bin_center(x, xMin, xMax, nBins)
    binWidth = (xMax - xMin) / nBins
    if x == 0.0
        return binWidth/2
    end
    return ceil(x/binWidth)*binWidth - binWidth/2
end


function get_rms_set(_tree, _dϕ, _ds = 1)
    RMS = Vector{Vector{Float64}}(undef, Int(180/_dϕ))                   # container for the RMS vectors (each row represnts one slice dϕ)
    S = Vector{Vector{Float64}}(undef, Int(180/_dϕ))                     # container for the shifts vectors (each row represnts one slice dϕ)
    dfG = get_gs_df(_tree, _dϕ, "p")

    for (i, n) in enumerate(1:_dϕ:180)
        cutEdges1 = get_cut_edges(n - 1, 1, _dϕ, "p")                    # provides the lower and upper cut 
        shifts    = -180:_ds:180                               # range of shifts is limited from left and right by empty diagonals, ie. for _dϕ = 0-5, shifts can go only [0-175)
        sdf       = @chain _tree begin                                   # filter out the dataframe
            @select(:thetaEscaped, :thetaEmitted)                       # keeps only the two angles columns
            @subset((cutEdges1[1] .<= :thetaEscaped .<= cutEdges1[2]))  # keeps only rows where ϕ is within the cut edges
        end
    
        rms = Vector{Float64}(undef, length(shifts))
    
        for j in 1:length(shifts)
            @rtransform! sdf :thetaModified = :thetaEscaped .+ shifts[j]
            @rsubset! sdf 0 .<= :thetaModified .<= 180
    
            fh2d = Hist2D(                                          
                (sdf.thetaEmitted, sdf.thetaModified),      
                (0:_dϕ:180, 0:_dϕ:180), 
                ) 
    
            g      = get_diagonal_sums(fh2d)

            rms[j] = get_rms(g, dfG[:,1])
    
        end
        RMS[i] = rms
        S[i]   = shifts

    end
    return RMS, S
end

function get_min_rms(_RMS )
    return [minimum(r) for r in _RMS]
end

function get_min_shifts(_S, _RMS)
    return [_S[i][argmin(_RMS[i])] for i in 1:length(_RMS)]
end

function get_phi_prime(_ϕ, _shifts)
    dϕ = Int(180/length(_shifts))
    bin = get_bin_center(_ϕ, Int(180/dϕ))
    return _ϕ + _shifts[Int(ceil(bin/dϕ))]
end

function get_hist_mean(_bincounts, _bincenters)
    N = sum(_bincounts)
    ss = 0.0
    for i in eachindex(_bincounts)
        ss += _bincounts[i] * _bincenters[i] 
    end
    return ss/N
end

function get_hist_mean(_hist)
    binCounts = bincounts(_hist)
    binCenters = bincenters(_hist)
    N = sum(binCounts)
    ss = 0.0
    for i in eachindex(binCounts)
        ss += binCounts[i] * binCenters[i] 
    end
    return ss/N
end

function get_hist_median(_hist)
    halfPoint = (sum(bincounts(_hist))+1)/2 # get the half value of the observations
    if (halfPoint == 0.0)
        return 0.0
    end

    cs = cumsum(bincounts(_hist))           # vector of cumulative sum of the bin counts (median is found in the bin, where cs[i] is more than halfPoint)       
    idx = 1

    while cs[idx] < halfPoint               # if the cumulative sum of the bin contents is at the edge (ie. halfPoint = 2, cs[2] = 2, we save the next bin -> 3)
        idx += 1
    end

    return bincenters(_hist)[idx]
end


end #MODULE END