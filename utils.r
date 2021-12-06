getPeaks_L = function(file,
                      peak_range = NA,
                      # time alignment
                      ref_file = NA,
                      # get MS1 peaks
                      MS1 = TRUE,
                      # get MS2
                      MS2 = FALSE) {
    require(xcms)
    require(parallel)
    
    # get adjustRtime, result is file.rtcor
    if (!is.na(ref_file)) {
        f.in = c(ref_file, file)
        raw_data = readMSData(
            files = f.in,
            mode = "onDisk",
            msLevel. = NULL,
            centroided = TRUE
        )
        rtcor = adjustRtime(
            raw_data,
            # values of below parameters come from IPO result
            param = ObiwarpParam(
                binSize = 0.9,
                distFun = "cor_opt",
                response = 5.5,
                gapInit = 0.4,
                gapExtend = 2.7,
                factorDiag = 2,
                factorGap = 1,
                localAlignment = FALSE,
                centerSample = 1,
            ),
            msLevel = 1L
        )
        rtcor = cbind(featureData(raw_data)@data, rtcor = rtcor)
        
        file.rtcor = rtcor[rtcor$fileIdx == 2, , drop = FALSE]
        file.rtcor[, "retentionTime"] = file.rtcor[, "rtcor"]
        file.rtcor = file.rtcor[,-ncol(file.rtcor)]
        file.rtcor[, "fileIdx"] = 1
        
        rm(f.in, raw_data)
    }
    
    # --- get peaks according to peak-range ---- #
    if (MS1) {
        peak_range = peak_range[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop =
                                    FALSE]
        mode(peak_range) = "numeric"
        
        file.xcmsRaw = xcmsRaw(filename = file)
        if (!is.na(ref_file)) {
            file.xcmsRaw@scantime = file.rtcor[file.rtcor$msLevel == 1, "retentionTime", drop =
                                                   TRUE]
        }
        peaks = getPeaks(file.xcmsRaw, peakrange = peak_range)
        peaks[, "mz"] = apply(peaks[, c("mzmin", "mzmax"), drop = FALSE], 1, median)
        peaks[, "rt"] = apply(peaks[, c("rtmin", "rtmax"), drop = FALSE], 1, median)
        peaks = cbind(peaks, sample = 1)
        # peaks = subset(peaks, subset=peaks[,"into"]!=0)
    } else {
        peaks = NA
    }
    
    # --- get MS2 ----- #
    if (MS2) {
        feature_table_MS2 = file.rtcor[file.rtcor$msLevel == 2, , drop = FALSE]
        aa = mzR::openMSfile(file)
        allPeaks_MS2 = mzR::peaks(aa)[feature_table_MS2[, "spIdx"]]
    } else {
        allPeaks_MS2 = NA
        feature_table_MS2 = NA
    }
    result = list(
        peaks = peaks,
        allPeaks_MS2 = allPeaks_MS2,
        feature_table_MS2 = feature_table_MS2,
        timeline = rtcor[, c("fileIdx", "retentionTime", "rtcor")]
    )
    return (result)
}


# ------ get peak-range from sample sets ----------- #
findPeaks_L = function(f.in, ref = NA, cwp = cwp) {
    cwp = cwp
    if (!is.na(ref)) {
        f.in = c(ref, f.in)
    }
    pd <-
        data.frame(
            sample_name = sub(
                basename(f.in),
                pattern = ".mzXML",
                replacement = "",
                fixed = TRUE
            ),
            sample_group = rep("sample", length(f.in)),
            stringsAsFactors = FALSE
        )
    
    # Read data and find peaks using centwave
    raw_data = readMSData(
        files = f.in,
        pdata = new("NAnnotatedDataFrame", pd),
        mode = "onDisk",
        msLevel. = NULL,
        centroided = TRUE
    )
    xdata = findChromPeaks(raw_data, param = cwp)
    
    if (is.na(ref)) {
        return(chromPeaks(xdata))
    }
    # ---------- Alignment retention time ------------- #
    ## Using obiwarp algorithm
    if (!is.na(ref)) {
        binSize = 0.9 # from IPO result
        xdata = adjustRtime(
            xdata,
            # values of below parameters come from IPO result
            param = ObiwarpParam(
                binSize = binSize,
                distFun = "cor_opt",
                response = 5.5,
                gapInit = 0.4,
                gapExtend = 2.7,
                factorDiag = 2,
                factorGap = 1,
                localAlignment = FALSE,
                centerSample = 1,
            ),
            msLevel = 1L
        )
        return (xdata)
        
    }
}

# getPeakRange from multiple files or matrix
getPeakRange_multipleFiles = function(f.in = NA,
                                      peaks = NULL,
                                      ref,
                                      cwp,
                                      cores = 2,
                                      absMz = 0.01,
                                      absRt = 15) {
    if (is.null(peaks)) {
        peaks_adjustRT = findPeaks_L(f.in = f.in,
                                     ref = ref,
                                     cwp = cwp)
    } else{
        peaks_adjustRT = peaks
    }
    
    pdp <-
        PeakDensityParam(
            sampleGroups = peaks_adjustRT$sample_group,
            minFraction = 0.25,
            bw = 1,
            binSize = 0.015
        )
    grouped_peaks <- groupChromPeaks(peaks_adjustRT, param = pdp)
    
    features <-
        as.matrix(featureDefinitions(
            grouped_peaks,
            mz = numeric(),
            rt = numeric(),
            ppm = 0,
            type = "any"
        ))
    new_feature_range <- t(apply(features, 1, function(x) {
        c(
            min(grouped_peaks@msFeatureData$chromPeaks[x[length(x)][[1]], "mzmin"]),
            max(grouped_peaks@msFeatureData$chromPeaks[x[length(x)][[1]], "mzmax"]),
            min(grouped_peaks@msFeatureData$chromPeaks[x[length(x)][[1]], "rtmin"]),
            max(grouped_peaks@msFeatureData$chromPeaks[x[length(x)][[1]], "rtmax"])
        )
    }))
    features[, c("mzmin", "mzmax", "rtmin", "rtmax")] <-
        new_feature_range
    colnames(features) <-
        c("mz",
          "mzmin",
          "mzmax",
          "rt",
          "rtmin",
          "rtmax",
          "npeaks",
          "sample",
          "id")
    features[, "id"] <- paste0(features[, "id"][[1]], collapse = ",")
    
    
    return(features)
    
}

convert_num_to_char = function(x, n = 4) {
    x = as.character(x)
    len = nchar(x)
    imp = rep("0", n - len)
    imp = paste(imp, collapse = "")
    char = paste("FT", imp, x, sep = "", collapse = "")
    return(char)
}


alignMS2_to_MS1 = function(allPeaks_MS2,
                           feature_table_MS2,
                           peak_range,
                           cores = 3) {
    require(parallel)
    
    peak_range = peak_range[, c("mzmin", "mzmax", "rtmin", "rtmax")]
    peak_range = cbind(peak_range, id = 1:nrow(peak_range))
    mode(peak_range) = "numeric"
    
    feature_table = feature_table_MS2[, c("precursorMZ", "retentionTime")]
    feature_table = cbind(feature_table, id = 1:nrow(feature_table))
    
    feature_table_order = order(feature_table[, "precursorMZ"])
    feature_table_sorted = feature_table[feature_table_order,]
    
    Match = rep(NA, nrow(feature_table))
    
    binary_search = function(x) {
        #browser()
        beg = 1
        end = nrow(feature_table)
        while (beg <= end) {
            mid = floor((beg + end) / 2)
            if (feature_table_sorted[mid, "precursorMZ"] >= x["mzmin"] &
                feature_table_sorted[mid, "precursorMZ"] <= x["mzmax"]) {
                # found one MS2 whose "precursorMZ" is within x["mzmin"] and x["mzmax"]
                if (feature_table_sorted[mid, "retentionTime"] >= x["rtmin"] &
                    feature_table_sorted[mid, "retentionTime"] <= x["rtmax"] &
                    is.na(Match[feature_table_sorted[mid, "id"]])) {
                    Match[feature_table_sorted[mid, "id"]] <<- x["id"]
                }
                #backward search
                idx = mid - 1
                while (idx >= beg) {
                    if (feature_table_sorted[idx, "precursorMZ"] >= x["mzmin"] &
                        feature_table_sorted[idx, "precursorMZ"] <= x["mzmax"]) {
                        if (feature_table_sorted[idx, "retentionTime"] >= x["rtmin"] &
                            feature_table_sorted[idx, "retentionTime"] <= x["rtmax"] &
                            is.na(Match[feature_table_sorted[idx, "id"]])) {
                            Match[feature_table_sorted[idx, "id"]] <<- x["id"]
                        }
                        idx = idx - 1
                    } else{
                        break
                    }
                }
                # forward search
                idx = mid + 1
                while (idx <= end) {
                    if (feature_table_sorted[idx, "precursorMZ"] >= x["mzmin"] &
                        feature_table_sorted[idx, "precursorMZ"] <= x["mzmax"]) {
                        if (feature_table_sorted[idx, "retentionTime"] >= x["rtmin"] &
                            feature_table_sorted[idx, "retentionTime"] <= x["rtmax"] &
                            is.na(Match[feature_table_sorted[idx, "id"]])) {
                            Match[feature_table_sorted[idx, "id"]] <<- x["id"]
                        }
                        idx = idx + 1
                    } else{
                        break
                    }
                }
                
                return(0)
            } else if (feature_table_sorted[mid, "precursorMZ"] < x["mzmin"]) {
                beg = mid + 1
            } else{
                end = mid - 1
            }
        }
        
        return(0)
    }
    
    apply(peak_range, 1, binary_search)

    # save(Match, file="Match.rData")
    
    result = cbind(feature_table_MS2, MS1 = Match)
    return(list(MS2 = allPeaks_MS2, MS2_to_MS1 = result))
    
}


get_MS2_cor = function(MS2_set1,
                       MS2_set2,
                       direction = c("forward", "reverse"),
                       cores = 2) {
    direction = match.arg(direction, c("forward", "reverse"))
    score = lapply(MS2_set1, function(x) {
        spec.exp = x
        colnames(spec.exp)[1:2] = c("mz", "intensity")
        score1 = lapply(MS2_set2, function(y) {
            if (is.null(y))
                return(0)
            spec.lib = y
            colnames(spec.lib)[1:2] = c("mz", "intensity")
            if (nrow(spec.exp) == 0 | nrow(spec.lib) == 0)
                return(0)
            
            score1 = GetMatchResult(
                spec.exp = spec.exp,
                spec.lib = spec.lib,
                direction = direction
            )[[1]]
            return(score1)
        })
        score2 = unlist(score1)
        return(score2)
    })
    score = do.call(rbind, score)
    return(score)
}


add_adduct = function(lib.meta, lib.adduct) {
    adduct = read.csv(lib.adduct, stringsAsFactors = FALSE)
    lib.mz <- t(sapply(lib.meta[, 'mz'], function(mz) {
        mz <- as.numeric(mz)
        apply(adduct, 1, function(info.adduct) {
            x <- gsub('\\(', '',
                      gsub('M.*', '', info.adduct['adduct']))
            xm <- ifelse(x == '', 1, as.numeric(x))
            xm * mz + as.numeric(info.adduct['mz'])
        })
    }))
    rownames(lib.mz) = lib.meta[, "labid"]
    colnames(lib.mz) <- adduct$adduct
    return(lib.mz)
}

ms2Match = function(db.MS2 = db.MS2,
                        db = file.path(ref_path, "lib", "rFT_meta_lib.RData"),
                        ce = 30,
                        mode = "pos",
                        absMz = 0.015,
                        lib.adduct) {
    load(db)
    
    comp = db$compound
    comp.mode = comp[[mode]]
    comp.mode.ce = sapply(comp.mode, function(x) {
        x[[as.character(ce)]]
    })
    
    meta = db$meta
    meta.lib = meta$compound
    rownames(meta.lib) = meta.lib[, "labid"]
    meta.precursor = add_adduct(lib.meta = meta.lib, lib.adduct = lib.adduct)
    
    MS2_to_MS1 = db.MS2$MS2_to_MS1
    MS1_index = unique(MS2_to_MS1$MS1)
    
    MS1_index = MS1_index[!is.na(MS1_index)]
    
    score = mclapply(MS1_index, function(idx) {
        MS2.idx = which(MS2_to_MS1$MS1 == idx)
        MS2.idx.mz = median(MS2_to_MS1[MS2.idx, "precursorMZ"])
        
        mz.filter = abs(MS2.idx.mz - as.numeric(meta.precursor)) <= absMz
        mz.filter = which(mz.filter)
        axe.col = ceiling(mz.filter / nrow(meta.precursor))
        axe.row = mz.filter - (axe.col - 1) * nrow(meta.precursor)
        
        db.temp = comp.mode.ce[rownames(meta.precursor)[axe.row]]
        
        idx.MS2 = db.MS2$MS2[MS2.idx]
        
        # compare db.MS2 and idx.MS2
        score = get_MS2_cor(MS2_set1 = idx.MS2,
                            MS2_set2 = db.temp,
                            direction = "forward")
        return(list(
            score = score,
            adduct = colnames(meta.precursor)[axe.col]
        ))
    }, mc.cores = 8)
    
    # voting
    ID = lapply(score, function(item, cutoff = 0.5) {
        df = item$score
        if (is.null(df))
            return(c(NA, NA, NA))
        
        score_median = apply(df, 2, function(x) {
            sum(x, na.rm = T) / length(x)
        })
        if (mean(is.na(score_median)) == 1)
            return(c(NA, NA, NA))
        
        max_index = which(score_median == max(score_median, na.rm = T))
        score_max = score_median[max_index]
        if (score_max > cutoff) {
            return(c(
                as.numeric(score_max),
                names(score_max),
                item$adduct[max_index]
            ))
        } else {
            return(c(NA, NA, NA))
        }
    })
    
    ID = do.call(rbind, ID)
    
    # get anno of MS1
    anno = cbind(MS1_index, ID)
    colnames(anno) = c("MS1_index", "score", "labid", "adduct")
    
    anno = anno[!is.na(anno[, "labid"]),]
    metabolite = meta.lib[anno[, "labid"], "name"]
    anno = cbind(anno, metabolite)
    rownames(anno) = 1:nrow(anno)
    
    return(anno)
}

GetMatchResult <- function(spec.exp,
                           spec.lib,
                           weight.int = 1,
                           weight.mz = 0,
                           ppm.ms2match = 30,
                           mz.ppm.thr = 400,
                           ppm.sanity.check = 100,
                           is.sanity.check = FALSE,
                           direction = 'forward',
                           ...) {
    # browser()
    if (is.sanity.check) {
        switch(direction,
               'forward' = {
                   if (any(c(
                       GetDiffMZppm(spec.exp[, 'mz']),
                       GetDiffMZppm(spec.lib[, 'mz'])
                   ) <= ppm.sanity.check)) {
                       stop('Difference between m/z is too small!!')
                   }
               },
               'reverse' = {
                   if (any(GetDiffMZppm(spec.lib) <= ppm.sanity.check)) {
                       stop('Difference between m/z is too small!!')
                   }
               },
               stop('Error setup for parameter: direction!!!'))
    }
    # browser()
    
    spec2match <-
        GetSpec2Match(
            spec.exp,
            spec.lib,
            # GetSpec2Match returns match(spec.exp, spec.exp+spec.lib) for "exp"
            ppm.ms2match = ppm.ms2match,
            # and match(spec.lib, spec.exp+spec.lib) + spec.exp for "lib"
            mz.ppm.thr = mz.ppm.thr,
            direction = direction
        )
    int.weighted.pk  <-
        GetWeightedInt(spec2match$exp, weight.mz, weight.int)
    int.weighted.lib <-
        GetWeightedInt(spec2match$lib, weight.mz, weight.int)
    match.score <-
        GetDotProduct(int.weighted.pk, int.weighted.lib)
    attr(match.score, 'spec') <- spec2match
    attr(match.score, 'spec.compared') <-
        cbind('exp' = int.weighted.pk,
              'lib' = int.weighted.lib)
    return(match.score)
}


GetSpec2Match <- function(spec.exp,
                          spec.lib,
                          ppm.ms2match = 30,
                          mz.ppm.thr = 400,
                          direction = c('reverse', 'forward')) {
    direction = match.arg(direction)
    mz.pool   <- sort(c(spec.exp[, 'mz'], spec.lib[, 'mz']))
    spec.temp <-
        cbind('mz' = mz.pool, 'intensity' = 0)
    
    spec.exp.temp  <-
        MatchFromTemp(spec.exp, spec.temp)
    spec.lib.temp <-
        MatchFromTemp(spec.lib, spec.temp)
    
    
    pk.spec  <-
        MatchSpec(spec.exp.temp,
                  ppm.ms2match = ppm.ms2match,
                  mz.ppm.thr = mz.ppm.thr)
    lib.spec <-
        MatchSpec(spec.lib.temp,
                  ppm.ms2match = ppm.ms2match,
                  mz.ppm.thr = mz.ppm.thr)
    
    if (direction == 'reverse') {
        idx.rm <- which(lib.spec[, 'intensity'] == 0)
        if (length(idx.rm > 0)) {
            pk.spec  <- pk.spec[-idx.rm, , drop = FALSE]
            lib.spec <- lib.spec[-idx.rm, , drop = FALSE]
        }
    }
    
    return(list('exp' = pk.spec, 'lib' = lib.spec))
}

MatchFromTemp <- function(spec, temp) {
    temp[match(spec[, 'mz'], temp[, 'mz']), 'intensity'] <-
        spec[, 'intensity']
    temp
}

MatchSpec <- function(spec,
                      ppm.ms2match = 30,
                      mz.ppm.thr = 400) {
    while (TRUE) {
        mz.diff.ppm <- GetDiffMZppm(spec[, 'mz'], mz.ppm.thr = mz.ppm.thr)
        idx <- which(mz.diff.ppm < ppm.ms2match)
        if (length(idx) > 0) {
            i <- tail(idx, 1)
            j <- which.max(spec[c(i, i + 1), 'intensity'])
            spec[i, 'intensity'] <- spec[i + j - 1, 'intensity']
            i2 <- i + 1
            spec[i, 'mz'] <- spec[i2, 'mz']
            spec <- spec[-i - 1, , drop = FALSE]
        } else {
            break
        }
    }
    return(spec)
}

GetDiffMZppm <- function(mz, mz.ppm.thr = NULL) {
    mz.diff <- diff(mz) / mz[-1] * 1e6
    if (!is.null(mz.ppm.thr)) {
        idx <- which(mz[-1] <= mz.ppm.thr)
        mz.diff[idx] <- mz.diff[idx] * mz[-1][idx] / mz.ppm.thr
    }
    mz.diff
}

GetWeightedInt <- function(spec,
                           weight.mz = 0,
                           weight.int = 1) {
    return(spec[, 'mz'] ^ weight.mz * spec[, 'intensity'] ^ weight.int)
}

GetDotProduct <- function(x, y) {
    return(sum(x * y) / sqrt((sum(x ^ 2) * sum(y ^ 2))))
}
