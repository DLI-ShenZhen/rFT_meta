load_packages = function() {
    require(xcms)
    require(RColorBrewer)
    require(magrittr)
    require(tools)
}
suppressPackageStartupMessages(load_packages())
source("utils.R")
total_core_num = detectCores()
max_used_core_num = ifelse(is.na(total_core_num), 1, floor(total_core_num /
                                                               2))
if (max_used_core_num < 1)
    max_used_core_num = 1

options(stringsAsFactors = FALSE, mc.cores = max_used_core_num)

# set parameters
{
    # values of below parameters came from IPO result
    peakwidth = c(6.12, 38.4)
    ppm = 15
    noise = 0
    snthresh = 10
    mzdiff = -9.99999e-06
    prefilter = c(3, 1000)
    mzCenterFun = "wMean"
    integrate = 1
    fitgauss = FALSE
    
    ## Using CentWave algorithm to find peaks
    cwp = CentWaveParam(
        peakwidth = peakwidth,
        ppm = ppm,
        noise = noise,
        snthresh = snthresh,
        mzdiff = mzdiff,
        prefilter = prefilter,
        mzCenterFun = mzCenterFun,
        integrate = integrate,
        fitgauss = fitgauss
    )
}

# set out dir
out = "rFT_meta_result"
if (!dir.exists(out)) {
    dir.create(out)
}

# read ref_files
ref_path = file.path("database")
ref_sample = file.path(ref_path, "reference.mzXML")
# read files
path = file.path(".")
files = dir(path)[grep("mzXML$", dir(path))]
files = file.path(path, files)

# make reference
if (!file.exists(file.path(ref_path, "ref_peaks.rData"))) {
    features = getPeakRange_multipleFiles(
        f.in = files,
        ref = ref_sample,
        cwp = cwp,
        cores = max_used_core_num,
        absMz = 0.01,
        absRt = 15
    )
    save(features, file = file.path(ref_path, "ref_peaks.rData"))
} else{
    load(file.path(ref_path, "ref_peaks.rData"))
}

if (!file.exists(file.path(ref_path, "ref_MS2.rData"))) {
    # get MS2
    MS2 = list()
    for (file in files) {
        MS2[[file]] = getPeaks_L(
            file = file,
            ref_file = ref_sample,
            peak_range = NA,
            MS1 = FALSE,
            MS2 = TRUE
        )
    }
    
    allPeaks_MS2 = c()
    for (i in MS2) {
        allPeaks_MS2 = c(allPeaks_MS2, i$allPeaks_MS2)
    }
    
    feature_table_MS2 = lapply(1:length(MS2), function(i) {
        sampleName = names(MS2)[i]
        df = MS2[[i]]$feature_table_MS2
        df = df[, c("precursorMZ",
                    "retentionTime",
                    "basePeakMZ",
                    "basePeakIntensity",
                    "spIdx")]
        df = cbind(df, sample = sampleName)
        return(df)
    })
    feature_table_MS2 = do.call(rbind, feature_table_MS2)
    
    save(MS2,
         allPeaks_MS2,
         feature_table_MS2,
         file = file.path(ref_path, "ref_MS2.rData"))
} else{
    load(file.path(ref_path, "ref_MS2.rData"), verbose = T)
}

if (!file.exists(file.path(ref_path, "db.MS2.rData"))) {
    db.MS2 = alignMS2_to_MS1(
        allPeaks_MS2 = allPeaks_MS2,
        feature_table_MS2 = feature_table_MS2,
        peak_range = features,
        cores = 3
    )
    save(db.MS2, file = file.path(ref_path, "db.MS2.rData"))
} else{
    load(file.path(ref_path, "db.MS2.rData"), verbose = T)
}

args <- commandArgs(TRUE)
polarity = args[1]
if (!file.exists(file.path(ref_path, "MS1_anno.rData"))) {
    if (polarity == "pos") {
        lib.adduct = file.path(ref_path, "list", "RP_POS.csv")
    } else {
        lib.adduct = file.path(ref_path, "list", "RP_NEG.csv")
    }
    MS1_anno = ms2Match(
        db.MS2 = db.MS2,
        db = file.path(ref_path, "lib", "rFT_meta_lib.RData"),
        ce = 30,
        mode = polarity,
        absMz = 0.015,
        lib.adduct = lib.adduct
    )
    
    save(MS1_anno, file = file.path(ref_path, "MS1_anno.rData"))
} else{
    load(file = file.path(ref_path, "MS1_anno.rData"))
}


# get intensity
intens = list()
for (file in files) {
    intens[[file]] = getPeaks_L(file = file,
                                ref_file = ref_sample,
                                peak_range = features)
}


## RESULT
result = sapply(intens, function(x) {
    x$peaks[, "into"]
})
colnames(result) = sapply(colnames(result), function(x) {
    strsplit(basename(x), split = ".", fixed = TRUE)[[1]][1]
})

name = sapply(1:nrow(result), convert_num_to_char, n = 5)

MS1_anno = as.matrix(MS1_anno)


anno = matrix(nrow = nrow(result), ncol = 4)
colnames(anno) = colnames(MS1_anno)[2:5]
anno[as.numeric(MS1_anno[, 1]),] = MS1_anno[, 2:5]
anno[is.na(anno)] = "."


result = cbind(name = name, 
               as.data.frame(features[, c("mz",
                                          "mzmin",
                                          "mzmax",
                                          "rt",
                                          "rtmin",
                                          "rtmax")]), 
               as.data.frame(result), anno)
rownames(result) = name

save(result, file = file.path(out, "intensity.rData"))
library(data.table)
fwrite(result, file = file.path(out, "rFT_meta_results.csv"))
