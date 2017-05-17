rm(list=ls())

Sys.setenv(TZ="UTC")
Sys.getenv("TZ")
library("mops")
library("ppso")



# setwd("/home/koehn/Shiva/modeling/echse_prj/hypsoRR/mahanadi/Lisei_middle_reaches_upstream/Calibration")
setwd("e:/data/kra/echse/Calibration_Ann")
periods= data.frame(begin=ISOdatetime(2001,3,2,8,0,0),end=ISOdatetime(2010,12,31,23,0,0))

periods.ignore= data.frame(begin=ISOdatetime(1988,1,1,0,0,0),end=ISOdatetime(1989,12,31,0,0,0))
################################################################################
# SETTINGS #####################################################################
################################################################################
#wdir          = "/home/koehn/Shiva/modeling/echse_prj/hypsoRR/mahanadi/Lisei_middle_reaches_upstream/Calibration/" # wdir war und bleibt auskommentiert
#outdir        = "imd4/Basantpur/out" # out odrner in Calibration ordner
outdir="calibration/out"
#mymodelpath   = "/home/koehn/Shiva/modeling/echse/echse_engines/bin/hypsoRR" #
mymodelpath = "e:/data/kra/echse/Echse/echse_engines/bin/hypsoRR"

#myrangetable  = "imd4/Basantpur/tbl_ranges_big.txt" # calibration ordner
myrangetable = "calibration/tbl_ranges_bigger.txt"

#sim_file      = c(qx_avg="imd4/Basantpur/out/gage_Basantpur.txt")#,end_of_interval="end_of_interval",qx_avg="qx_avg") 
sim_file = c(qx_avg="calibration/out/gage_Basantpur.txt")
# gage_Basantpur.txt habe ich nicht, nur gage_Basantpur.csv -> richtig?oder wird das erstellt?
sim_colTime   = c("end_of_interval")
sim_colValue  = c("qx_avg")
obs_file      = read.table("calibration/flow_Basantpur_0.txt", sep="\t", stringsAsFactors=FALSE, colClasses=c("POSIXct","numeric"), header=TRUE)

obs_colTime   = "end_of_interval"

#mops_log      = "imd4/Basantpur/mcs_and_pso.log"
mops_log = "calibration/mcs_and_pso.log" 
#ppso_log      = "imd4/Basantpur/dds.log"  
ppso_log = "calibration/dds.log" 
#ppso_proj     = "imd4/Basantpur/dds.pro"
ppso_proj ="calibration/dds.pro"  

parameter_fmt = "" # set to "" for ECHSE and "%13.4f" for HBV Nordic
if_fails_clear= c("out", "phiMod.err.html", "phiMod.log") #for ECHSE: c("out", "phiMod.err.html", "phiMod.log")


model_args= c(
  #file_control= "imd4/Basantpur/cnf_imd4_Basantpur",
  file_control="calibration/cnf_imd4_Basantpur_Ann.txt",
  file_err= paste(outdir,"cnf_raingages.err.html",sep="/"),
  file_log= paste(outdir,"cnf_raingages.log",sep="/"),
  format_err= "html", silent= "true",
  #table_objectDeclaration= "/home/koehn/Shiva/modeling/echse_prj/hypsoRR/mahanadi/Lisei_middle_reaches_upstream/data/topocatch/out_csa100_Arc_Gis_splitting2/Basantpur/objDecl.gage_Basantpur.txt",
  table_objectDeclaration="out_csa100_Arc_Gis_splitting/Basantpur/objDecl.gage_Basantpur.txt",
  #cat_hypsoRR_numParamsIndividual= "imd4/Basantpur/paramNum_cat_MCS_template_updated_mult_str_surf.txt",
  cat_hypsoRR_numParamsIndividual= "calibration/paramNum_cat_MCS_template_updated_mult_str_surf.txt",
  #cat_hypsoRR_numParamsShared= "cat_num_shared.txt",
  cat_hypsoRR_numParamsShared="../echse_projekt_Ann/data/params/cat_num_shared.txt",
  outputDirectory= outdir,
  #file_template="imd4/Basantpur/paramNum_cat_MCS_template2_new.txt",
  file_template="calibration/paramNum_cat_MCS_template2_new.txt",
  #file_result="imd4/Basantpur/paramNum_cat_MCS_template_updated_mult_str_surf.txt",
  file_result= "calibration/paramNum_cat_MCS_template_updated_mult_str_surf.txt",
    char_open="{",
  char_close="}"
)


# goodness of fit function: modified nash-sutcliffe efficiency
gof_func      = function(obs, sim) {
  NSE = 1 - mean((sim-obs)^2)/var(obs)
  pBias = (sum(sim-obs)/sum(obs))*100
  #mNSE = (NSE * ((100-abs(pBias))/100))*(-1)
  rmse= sqrt(mean((sim-obs)^2))
  ifelse(NSE>0,mNSE <- (NSE * ((100-abs(pBias))/100))*(-1),mNSE <- (NSE * ((100+abs(pBias))/100))*(-1))
  return(mNSE)
}

gof_func2      = function(obs, sim) {
  NSE = 1 - mean((sim-obs)^2)/var(obs)
  pBias = (sum(sim-obs)/sum(obs))*100
  #mNSE = (NSE * ((100-abs(pBias))/100))*(-1)
  rmse= sqrt(mean((sim-obs)^2))
  ifelse(NSE>0,mNSE <- (NSE * ((100-abs(pBias))/100))*(-1),mNSE <- (NSE * ((100+abs(pBias))/100))*(-1))
  all = c(mNSE=mNSE, pBias=pBias, rmse=rmse, NSE=NSE)
  return(all)
}


################################################################################
# END SETTINGS #################################################################
################################################################################

#setwd(wdir)

# Process arguments
updatePars_byTemplate= function(parameters, model_args) {
  
  stopifnot(all(c("file_template","file_result","char_open","char_close")
                %in% names(model_args)))
  n= update_template(file_template=model_args[10], vect_updates=parameters,
                     char_open=model_args[12], char_close=model_args[13],
                     file_result=model_args[11], overwrite=TRUE)
  if (n == 0)
    stop(paste("Updating of '",model_args[10],"' failed.",sep=""))
  
  return(0)
}



save.output = function(parameters, gof, moreArgs_final) {
  
  if (!dir.exists("calibration/out")) stop("Directory does not exist.")
  files= list.files(path="calibration/out", pattern=".+[.]txt?", full.names=TRUE)
  for (f in files) file.remove(f)
  files= list.files(path="calibration/out", pattern=".+[.]log?", full.names=TRUE)
  for (f in files) file.remove(f)
  files= list.files("calibration/out", pattern=".+[.]html?", full.names=TRUE)
  for (f in files) file.remove(f)
  #    
  return(0)      # debug
}
# bei lisei  calibration/out = imd4/Basantpur/out


objfunc = function(parameters) {
  # format parameters as needed by HBV subroutine PARINN
  if (parameter_fmt!="") {
    parameters = sprintf(parameter_fmt, parameters)
  }
  #   # set parameter names in parameter vector
  attr(parameters, "names")=read.table(file=myrangetable, header=TRUE, sep="\t", colClasses= c("character","numeric","numeric"))[,1] 
  print (paste(Sys.time(),": objfunc called..."))
  # modelError_multiDim updates parameter file, calls the model, reads the results and returns GOF
  error = modelError_multiDim(
    parameters=parameters,
    model_path = mymodelpath,
    model_args = model_args,
    func_first = updatePars_byTemplate, 
    moreArgs_first = model_args, #removeInfo.args,
    func_final =save.output, 
    moreArgs_final = removeOutput.args,
    sim_files = sim_file,
    sim_colsTime= sim_colTime, 
    sim_colsValue= sim_colValue,
    sim_colsep = "\t",
    sim_timeConv = function(x) { as.POSIXct(strptime(x,"%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = "UTC") },
    observed= obs_file,
    obs_colTime= obs_colTime, 
    periods.select=periods,
    periods.ignore=periods.ignore,
    gof_function = gof_func, 
  )
  print (paste("Error:",error))
  flush.console()
  return (error[["qx_avg"]])
}


# Just dummy first/final functions to be used whenever needed
dummyfuncfinal = function(parameters, gof, args) {
  gc(verbose=FALSE)
  return(0)}
dummyfuncfirst = function(parameters, args) {
  gc(verbose=FALSE)
  return(0)}


# Read one catchment (cat) file and catch possible errors
readcat <- function(f) {
  out <- tryCatch(
    {
      read.table(f, sep="\t", stringsAsFactors = FALSE, header=TRUE,
                       colClasses = c("character", "numeric") ) 
    },
    error=function(cond) {
      message(paste("Error in reading file:", f))
      #message(cond)
      message("Continue without that file.")
      return(NA)
    },
    warning=function(cond) {
      message(paste("File caused a warning:", f))
      message("Here's the original warning message:")
      #message(cond)
      message("Continue without that file.")
      # Choose a return value in case of warning
      return(NA)
    },
    finally={
      # Finally: nothing else.
    }
  )    
  return(out)
}


# This function computes the areal average of monthly ETR
#    from the model's cat_*.txt output files.
#    In order to use it as first/final function in mops,
#    it has to use the arguments "parameters, gof, moreArgs" 
process_catfiles = function(parameters, gof, args) {
  
  # moreArgs
  #    pattern="calibration/out/cat_*.txt", 
  #    timecol="end_of_interval", 
  #    varcol="etr",
  #    outfile="calibration/out/etr_monthly.txt"
  
  # Find output files and initialize containers
  print("Post-processing of cat-files with ETR.")
  files = Sys.glob(args["pattern"] )
  files = sort(files)
  vals = c()
  colnames = c()
  # Iterate over files
  for (f in files) {
    ##print(f)
    cat = readcat(f)
    if (!is.data.frame(cat)) {
      next
    }
    
    if (length(vals)==0) {
      vals = cat[,args["varcol"]]
      dtimes = as.POSIXct(cat[,args["timecol"]], format="%Y-%m-%d %H:%M:%S")
    } else {
      vals = cbind(vals, cat[,args["varcol"]])
    }
    colnames = c(colnames, tools::file_path_sans_ext(f)) 
  }
  names(vals) = colnames
  # Units: from m/s (over an hour) to mm
  vals = vals * 1e3 * 3600
  # Compute average value over all subcatchments
  avg = apply(X=vals, MARGIN=1, FUN=mean, na.rm=TRUE)

  # Compute monthly sums
  months = strftime(dtimes, format="%Y-%m-01 00:00:00")
  monthly = aggregate(avg, by=list(months), "sum", na.rm=TRUE)
  names(monthly) = c("month", "etr")
  
  # Write results file
  write.table(monthly, args["outfile"], quote=FALSE, sep = "\t", row.names=FALSE)
  
  gc(verbose=FALSE)
  return(0)
}


# Argument vector for process_catfiles
process_catfiles_args= c(
  pattern="calibration/out/cat_*.txt", 
  timecol="end_of_interval", 
  varcol="etr",
  outfile="calibration/out/etr_monthly.txt"
)

## DEBUGGING
##readcat("calibration/out/cat_460.txt")
##process_catfiles(NULL, NULL, process_catfiles_args)
##plot(as.POSIXct(out$month), out$etr, type="l")


# Compute areal everage ETR from MODIS data
#   This will produce the "observed" data.frame 
#   to be used with error2 in objfunc2.
etr_from_modis = function(infile="modis_etr/data.txt") {
  cats = read.table(infile, sep="\t", stringsAsFactors = FALSE, header=TRUE,
                   na.strings = "nan" )
  # Compute average over all catchment
  avg = apply(X=cats[,2:ncol(cats)], MARGIN=1, FUN=mean, na.rm=TRUE)
  df = data.frame(month=as.POSIXct(cats$datetime), etr=avg)
  gc(verbose=FALSE)
  return(df)
}

# This objective function attempts to consider both discharge and ETR
objfunc2 = function(parameters) {
  # format parameters as needed by HBV subroutine PARINN
  if (parameter_fmt!="") {
    parameters = sprintf(parameter_fmt, parameters)
  }
  #   # set parameter names in parameter vector
  attr(parameters, "names")=read.table(file=myrangetable, header=TRUE, sep="\t", colClasses= c("character","numeric","numeric"))[,1] 
  print (paste(Sys.time(),": objfunc2 called..."))
  # 1st call
  error1 = modelError_multiDim(
    parameters=parameters,
    model_path = mymodelpath,
    model_args = model_args,
    func_first = updatePars_byTemplate,
    moreArgs_first = model_args,
    func_final =process_catfiles,
    moreArgs_final = process_catfiles_args,
    sim_files = sim_file,
    sim_colsTime= sim_colTime,
    sim_colsValue= sim_colValue,
    sim_colsep = "\t",
    sim_timeConv = function(x) { as.POSIXct(
      strptime(x,"%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = "UTC") },
    observed= obs_file,
    obs_colTime= obs_colTime,
    periods.select=periods,
    periods.ignore=periods.ignore,
    gof_function = gof_func2
  )
  gc(verbose=FALSE)
  # 2nd call
  error2 = modelError_multiDim(
    parameters=parameters,
    model_path = "calibration/dummymodel.cmd",
    model_args = NULL,
    func_first = dummyfuncfirst, 
    moreArgs_first = NULL, 
    func_final =save.output, 
    moreArgs_final = removeOutput.args,
    sim_files = c(etr="calibration/out/etr_monthly.txt"),
    sim_colsTime= c("month"), 
    sim_colsValue= c("etr"),
    sim_colsep = "\t",
    sim_timeConv = function(x) { as.POSIXct(
      strptime(x,"%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = "UTC") },
    observed= etr_from_modis(),
    obs_colTime= "month", 
    periods.select=periods,
    periods.ignore=periods.ignore,
    gof_function = gof_func2 
  )
  gc(verbose=FALSE)
  error1 = error1[["qx_avg"]]
  error2 = error2[["etr"]]
  print (paste("NSE1=",error1["NSE"], "; pBias1=", 
               error1["pBias"], "; mNSE=", error1["mNSE"]))
  print (paste("NSE2=",error2["NSE"], "; pBias2=", 
               error2["pBias"], "; mNSE2=", error2["mNSE"]))
  flush.console()
  error = 0.8*error1["mNSE"] + 0.2*error2["mNSE"]
  print (paste("Error:",error))
  gc(verbose=FALSE)
  return (error)
}





###############################################################

# Call the optimisation routine from ppso package, using our objective function
# read paramter ranges
param_bounds = read.table(file=myrangetable, 
                          header=TRUE, 
                          sep="\t", 
                          colClasses= c("character","numeric","numeric"))[,2:3]
rownames(param_bounds) = read.table(file=myrangetable, 
                                    header=TRUE, 
                                    sep="\t", 
                                    colClasses= c("character","numeric","numeric"))[,1]


# TEST CALL (also a way you can just call the model ONCE)
##out = objfunc2(apply(X=param_bounds, MARGIN=1, FUN=mean))


result = optim_dds(
  objective_function = objfunc2,
  number_of_parameters = length(param_bounds[,1]),
  number_of_particles =  1,
  max_number_function_calls= 100,#max_number_function_calls, #### ERstmal weniger läufe z.B. 10 probieren 
  r=0.2,
  abstol = -Inf, 
  reltol = -Inf, 
  max_wait_iterations=50,
  parameter_bounds = param_bounds,
  lhc_init=FALSE, 
  do_plot = NULL, 
  wait_for_keystroke = FALSE,
  logfile=  "calibration/dds.log",  
  projectfile = "calibration/dds.pro",
  load_projectfile = "yes", 
  break_file=NULL, 
  plot_progress=FALSE, 
  tryCall=FALSE)



# # plot progress
# plot_optimization_progress(logfile=  ppso_log,
#                            projectfile = ppso_proj,
#                            progress_plot_filename=NULL,
#                            goodness_plot_filename=NULL, 
#                            cutoff_quantile=0.95, 
#                            verbose=FALSE)
# 
# 
# 
# 
# 
# 
# ######## produce data for plotting
# save.output = function(parameters, gof, moreArgs_final) {
#   
#   return(0)      # debug
# }
# 
# sample_params = t(as.data.frame(read.table("calibration/dds.pro", header=FALSE, skip=1, sep="\t")[,1:12]))
# #rownames(sample_params)<- read.table("imd4/Basantpur/tbl_ranges_big.txt", header=TRUE, sep="\t", colClasses= c("character","numeric","numeric"))[,1]
# #rownames(sample_params)<- read.table("calibration/tbl_ranges_big.txt", header=TRUE, sep="\t", colClasses=c("character", "number", "numeric"))[,1] # funktioniert nicht
# rownames(sample_params)<- read.table("calibration/tbl_ranges_bigger.txt", header=TRUE, sep="\t", colClasses=c("factor", "numeric", "numeric"))[,1]
# gof_func      = function(obs, sim) {
#   NSE = 1 - mean((sim-obs)^2)/var(obs)
#   #NSE2 = (-1. + mean((obs - sim)^2) / var(obs) )
#   pBias = (sum(sim-obs)/sum(obs))*100
#   #mNSE = (NSE * ((100-abs(pBias))/100))*(-1)
#   rmse= sqrt(mean((sim-obs)^2))
#   ifelse(NSE>0,mNSE <- (NSE * ((100-abs(pBias))/100))*(-1),mNSE <- (NSE * ((100+abs(pBias))/100))*(-1))
#   #mNSE = (NSE * ((100-abs(pBias))/100))*(-1)    # modified NS
#   #mNSE = (NSE * ((100-abs(pBias))/100))*(-1)    # modified NS
#   #return(mNSE)
#   return(list(mNSE, NSE, pBias, rmse)) 
#   #return(NSE*(-1))
# }
# 
# test = objfunc(sample_params)
# #big <- read.table("calibration/tbl_ranges_big.txt", header=T, sep="\t")
# #big
# #class(big[,3])
