

# Setup -------------------------------------------------------------------


rm(list=ls())
setwd("~/Server/Excel/")

if(system.file(package="readODS") == "") install.packages("readODS")
library(readODS)


# Data --------------------------------------------------------------------


s = "1"
origin = paste0("/Scenario_", s, "/")                                      #<---------
load_path = paste0("../Scenarios", origin)
out_list = c("out1", "out2", "out3", "out4", "out5", "out6", "out7", "out8", "out9", "out10", 
             "out11", "out12", "out20", "out21", "out22", "out23")
for (i in out_list) {
  print(paste0(load_path, sep = i, ".rda"))
  load(paste0(load_path, sep = i, ".rda"))
}
keep = c("call", "ASB Post", "OR stats")


# Excel -------------------------------------------------------------------


o1 = get(out_list[1])
o2 = get(out_list[2])
init = mapply(rbind, o1[keep], o2[keep])
for (i in 3:length(out_list)) {
  oi = get(out_list[i]) 
  init = mapply(rbind, init, oi[keep])
}

excel_data = cbind(init[["call"]], round(init[["ASB Post"]], 3), round(init[["OR stats"]], 3))
excel_data = data.frame(excel_data)

og = excel_data[0,]
og[1, ] = NaN
og["call"] = "Original" 
og[2:(length(o1[["ASB Pre"]]) + 1)] = round(o1[["ASB Pre"]], 3)

excel_data = rbind(og, excel_data)


excel_name = paste0("/Excel_", s, ".ods")
write_ods(excel_data, path = paste0(load_path, excel_name), append = F, sheet = "Stats")
trt_stats = data.frame(o1[["Trt stats"]])
trt_stats = round(trt_stats, 3)
write_ods(trt_stats, path = paste0(load_path, excel_name), append = T, row_names = T, sheet = "Trt info")


# Graphs ------------------------------------------------------------------


for (i in out_list) {
  tmp = get(i)
  
  png(filename = paste0(load_path, "/Pre_Post_", s, "_", i, ".png"))
  print(tmp[["Pre-Post"]])
  dev.off()
  
  png(filename = paste0(load_path, "/Scatterplot_", s, "_", i, ".png"))
  print(tmp[["Scatterplot"]])
  dev.off()
}











