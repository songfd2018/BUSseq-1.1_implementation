rm(list=ls())
library(ggplot2)
setwd("Scalability")
if(!dir.exists("Image")){
  dir.create("Image")
}


# plot the curve of running time
# cpu time
gene_list <- c(1500, 3000, 4500, 6000, 9000)
cell_list_cpu <- c(1000, 2000, 4000, 10000, 20000)

# by gene
time_by_gene_cpu <- c(0.868, 1.768, 2.416, 3.411, 5.123)
RAM_by_gene_cpu <- c(0.229 , 0.442 , 0.661 , 0.867, 1.291)
dat_by_gene_cpu <- data.frame(GeneNum = gene_list,
                              RunningTime = time_by_gene_cpu,
                              RAMUsage = RAM_by_gene_cpu)

ratio <- round(max(time_by_gene_cpu)/max(RAM_by_gene_cpu))
jpeg(paste0("Image/cpu_RAM_time_by_gene.jpg"),width = 800, height = 600)
ggplot(dat_by_gene_cpu, aes(x = GeneNum)) +
  geom_line(aes(y = RAMUsage * ratio, colour = "RAM"), size = 1.5, linetype="dashed") +
  geom_point(aes(y = RAMUsage * ratio, colour = "RAM"), size = 4, shape = 16) +
  geom_line(aes(y = RunningTime, colour = "Computing Time"), size =1.5) +
  geom_point(aes(y = RunningTime, colour = "Computing Time"), size = 4, shape = 15) +
  scale_y_continuous(sec.axis = sec_axis(~./ ratio, name = "RAM Usage (GB)")) +
  #scale_colour_manual(values = c("#00BFC4", "#F8766D")) +
  labs(y = "Computing time (hour)",
       x = "Gene number",
       colour = "Variable") +
  guides(size = FALSE, linetype = FALSE, colour = guide_legend(override.aes = list(linetype = c("solid","dashed"), shape = c(15,16)))) +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 32), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 32),#, angle = 45))
        axis.title=element_text(size=36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 28,face="bold"),
        legend.key.width = unit(3, "cm"),
        legend.position="top")
dev.off()

# by cell
time_by_cell_cpu <- c(1.768, 3.263, 6.343, 16.426, 29.702 )
RAM_by_cell_cpu <- c(0.442, 0.509, 0.643, 1.046, 1.726)
dat_by_cell_cpu <- data.frame(CellNum = cell_list_cpu,
                              RunningTime = time_by_cell_cpu,
                              RAMUsage = RAM_by_cell_cpu)

ratio <- round(max(time_by_cell_cpu)/max(RAM_by_cell_cpu))
jpeg(paste0("Image/cpu_RAM_time_by_cell.jpg"),width = 800, height = 600)
ggplot(dat_by_cell_cpu, aes(x = CellNum)) +
  geom_line(aes(y = RAMUsage * ratio, colour = "RAM"), size = 1.5, linetype="dashed") +
  geom_point(aes(y = RAMUsage * ratio, colour = "RAM"), size = 4, shape = 16) +
  geom_line(aes(y = RunningTime, colour = "Computing Time"), size =1.5) +
  geom_point(aes(y = RunningTime, colour = "Computing Time"), size = 4, shape = 15) +
  scale_y_continuous(sec.axis = sec_axis(~./ ratio, name = "RAM Usage (GB)")) +
  #scale_colour_manual(values = c("#00BFC4", "#F8766D")) +
  labs(y = "Computing time (hour)",
       x = "Cell number",
       colour = "Variable") +
  guides(size = FALSE, linetype = FALSE, colour = guide_legend(override.aes = list(linetype = c("solid","dashed"), shape = c(15,16)))) +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 32), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 32),#, angle = 45))
        axis.title=element_text(size=36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 28,face="bold"),
        legend.key.width = unit(3, "cm"),
        legend.position="top")
dev.off()

# by storage iterations
iter_list_cpu <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
time_by_iter_cpu <- c(1.711, 1.679, 1.718, 1.778, 1.752, 1.738, 1.726, 1.733, 1.768, 1.783)
RAM_by_iter_cpu <- c(0.062, 0.062, 0.065, 0.069, 0.081,0.100, 0.138, 0.252, 0.442, 0.799)
dat_by_iter_cpu <- data.frame(IterNum = iter_list_cpu,
                              RunningTime = time_by_iter_cpu,
                              RAMUsage = RAM_by_iter_cpu)
ratio <- 1
jpeg(paste0("Image/cpu_RAM_time_by_iter.jpg"),width = 1200, height = 600)
ggplot(dat_by_iter_cpu, aes(x = log(IterNum))) +
  geom_line(aes(y = RAMUsage * ratio, colour = "RAM"), size = 1.5, linetype="dashed") +
  geom_point(aes(y = RAMUsage * ratio, colour = "RAM"), size = 4, shape = 16) +
  geom_line(aes(y = RunningTime, colour = "Computing Time"), size =1.5) +
  geom_point(aes(y = RunningTime, colour = "Computing Time"), size = 4, shape = 15) +
  scale_y_continuous(sec.axis = sec_axis(~./ ratio, name = "log(RAM Usage) (log(GB))")) +
  #scale_colour_manual(values = c("#00BFC4", "#F8766D")) +
  labs(y = "log(Computing time) (log(hour))",
       x = "Logarithm of the iteration number per storage",
       colour = "Variable") +
  guides(size = FALSE, linetype = FALSE, colour = guide_legend(override.aes = list(linetype = c("solid","dashed"), shape = c(15,16)))) +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 32), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 32),#, angle = 45))
        axis.title=element_text(size=36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 28,face="bold"),
        legend.key.width = unit(3, "cm"),
        legend.position="top")
dev.off()

# gpu time
cell_list_gpu <- c(1000, 2000, 4000, 10000, 20000, 50000)

# by gene
time_by_gene_gpu <- c(0.194, 0.365, 0.535, 0.714, 1.057)
RAM_by_gene_gpu <- c(1.304, 1.712, 2.122, 2.532, 3.347)
dat_by_gene_gpu <- data.frame(GeneNum = gene_list,
                              RunningTime = time_by_gene_gpu,
                              RAMUsage = RAM_by_gene_gpu)

ratio <- 1/3
jpeg(paste0("Image/gpu_RAM_time_by_gene.jpg"),width = 800, height = 600)
ggplot(dat_by_gene_gpu, aes(x = GeneNum)) +
  geom_line(aes(y = RAMUsage * ratio, colour = "RAM"), size = 1.5, linetype="dashed") +
  geom_point(aes(y = RAMUsage * ratio, colour = "RAM"), size = 4, shape = 16) +
  geom_line(aes(y = RunningTime, colour = "Computing Time"), size =1.5) +
  geom_point(aes(y = RunningTime, colour = "Computing Time"), size = 4, shape = 15) +
  scale_y_continuous(sec.axis = sec_axis(~./ ratio, name = "RAM Usage (GB)")) +
  #scale_colour_manual(values = c("#00BFC4", "#F8766D")) +
  labs(y = "Computing time (hour)",
       x = "Gene number",
       colour = "Variable") +
  guides(size = FALSE, linetype = FALSE, colour = guide_legend(override.aes = list(linetype = c("solid","dashed"), shape = c(15,16)))) +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 32), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 32),#, angle = 45))
        axis.title=element_text(size=36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 28,face="bold"),
        legend.key.width = unit(3, "cm"),
        legend.position="top")
dev.off()


# by cell
time_by_cell_gpu <- c(0.365, 0.440, 0.588, 0.856, 1.240, 2.455)
RAM_by_cell_gpu <- c(1.712, 1.925, 2.347, 3.624, 5.745, 10.372)
dat_by_cell_gpu <- data.frame(CellNum = cell_list_gpu,
                              RunningTime = time_by_cell_gpu,
                              RAMUsage = RAM_by_cell_gpu)

ratio <- 1/4
jpeg(paste0("Image/gpu_RAM_time_by_cell.jpg"),width = 800, height = 600)
ggplot(dat_by_cell_gpu, aes(x = CellNum)) +
  geom_line(aes(y = RAMUsage * ratio, colour = "RAM"), size = 1.5, linetype="dashed") +
  geom_point(aes(y = RAMUsage * ratio, colour = "RAM"), size = 4, shape = 16) +
  geom_line(aes(y = RunningTime, colour = "Computing Time"), size =1.5) +
  geom_point(aes(y = RunningTime, colour = "Computing Time"), size = 4, shape = 15) +
  scale_y_continuous(sec.axis = sec_axis(~./ ratio, name = "RAM Usage (GB)")) +
  #scale_colour_manual(values = c("#00BFC4", "#F8766D")) +
  labs(y = "Computing time (hour)",
       x = "Cell number",
       colour = "Variable") +
  scale_linetype_manual(values=c("solid", "dashed")) + 
  scale_shape_manual(values = ) + 
  guides(size = FALSE, linetype = FALSE, colour = guide_legend(override.aes = list(linetype = c("solid", "dashed"), shape = c(15,16)))) +
  theme(axis.text.x = element_text(face = "bold", #color = "#993333", 
                                   size = 32), #angle = 45),
        axis.text.y = element_text(face = "bold", #color = "blue", 
                                   size = 32),#, angle = 45))
        axis.title=element_text(size=36,face="bold"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.text=element_text(size=24,face="bold"),
        legend.title = element_text(size = 28,face="bold"),
        legend.key.width = unit(3, "cm"),
        legend.position="top")
dev.off()