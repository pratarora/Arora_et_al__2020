#
#
# This is a program to calculate various parameters for E-cadherin
# quantifications-----
#
#
#
# "rawdata" contains data read from the csv file(provided by user). It
# should contain the following columns:
#
# 1) Genotype
# 2) Cell name
# 3) Mean
# 4) Length
# 5) Boundary -- If you are doing cell wise analysis(put same letter
# in the whole column in case analysis is cell wise and by boundaries)
#
#
# "stack.all" contains stack wise calculated data which includes:
# 1) Height at a particular stack
# 2) Normalized Height
# 3) Normalized Ecad
# 4) Absolute Ecad at each stack
# "cell.all" contains cell wise data. Which includes the following:
#
# 1) The total Height of the Cell
# 2) Total Ecad in the cell
#
# "genotype.all" contains
# genotype wise data which includes the
# following:
#
# 1) Mean apical perimeter
# 2) Median apical perimeter
# 3) Mean Total Ecad per cell
# 4) Median Total Ecad per cell
# 5) Mean Height of the cell
# 6) Median Height of the cell
# 7) First quantile of Height of cell
# 8) Third quantile of height of cell
#
# "locfit.model" contains the locfit data and the height at which the
# levels were predicted.
#
# "apical.perimeter" contains the apical most stack used for apical
# perimeter analysis
#
#
#
#
#clears memory--------
rm(list = ls())
#sets a working directory-----------
#This is where all your graphs will be saved
setwd(
"C:/Users/prata/Desktop/Test for paper/"
)
#libraries required for the program to work.-----------------
# List of libraries required for the program to run
pkg <- c("scales", "locfit","tidyr","purrr","dplyr","ggplot2","broom",
"data.table",
"RCurl","stringr",
"grid",
"gridExtra",
"RColorBrewer")
# Check if the libraries are installed and install libraries not present
new.pkg <- pkg[!(pkg %in% installed.packages())]
if (length(new.pkg)) {
install.packages(new.pkg)
}
rm("pkg","new.pkg")
#load libraries in the memory
library(scales)
library(locfit)
library(tidyr)
library(purrr)
library(dplyr)
library(ggplot2)
library(data.table)
library(RCurl)
library(stringr)
library(grid)
library(gridExtra)
library(RColorBrewer)
#stops warnings from being displayed
options(warn = -1)
# read raw data------------------
# the file should have the following components:
# Genotype
# Cell name Area
# Mean
# Boundary
#Folder and file names where your input file is saved
Folderadd <-
"C:/Desktop/Test for paper/"
filename <-
"Figure 1- source data 1"
format <- ".csv"
fulladdress <- paste(Folderadd, filename, format, sep = "")
rawdata <-
read.csv(fulladdress,
header = TRUE,
sep = ",") %>% na.omit()
# data arranging--------------------
# To have Genotype and Cell name, Boundary order as in the file (Otherwise
# it comes alphabetically which may not be desirable)
rawdata$Genotype <-
factor(rawdata$Genotype, levels = unique(rawdata$Genotype))
rawdata$Cell.name <-
factor(rawdata$Cell.name, levels = unique(rawdata$Cell.name))
rawdata$Boundary <-
factor(rawdata$Boundary, levels = unique(rawdata$Boundary))
#group by genotype and cell name calculate absolute ecad for each stack,
# and nest the data to help calculations on groups
by_cellname <- rawdata %>%
group_by(Genotype, Cell.name) %>%
mutate(AbsoluteEcad = Mean * Area) %>%
nest()
#function to find height
rawheight <- function(df) {
height <- vector("double", nrow(df))
for (i in seq_along(df$Mean)) {
height[i] <- (i - 1) * 0.28 #the height of the z stack
}
height
}
#function to normalize height
normheightfun <- function(df) {
nheight <- vector("double", nrow(df))
maxh <- nrow(df)
minh <- 1
for (i in seq_along(df$Mean)) {
nheight[i] <- ((i) - minh) / (maxh - minh)
}
nheight
}
#function to normalize ecad levels
normecadfun <- function(df) {
necad <- vector("double", nrow(df))
maxe <- max(df$Mean)
mine <- min(df$Mean)
for (i in seq_along(df$Mean)) {
necad[i] <- (df$Mean[i] - mine) * 100 / (maxe - mine)
}
necad
}
#dataframe which now has height per stack, normalized height, normalized
ecad per genotype-cell name
stack.all <- by_cellname %>%
mutate(
Stackheight = data %>% map(rawheight),
Normalizedheight = data %>% map(normheightfun),
NormalizedEcad = data %>% map(normecadfun),
Uniquecellname = paste(Genotype, Cell.name)
) %>% unnest()
#to know the total ecad in a cell and height of a cell
cell.all <- stack.all %>%
group_by(Genotype, Cell.name) %>%
summarize(Total_Ecad = sum(AbsoluteEcad),
Height_of_cell = max(Stackheight))
#to fit linear regression to find the abnormal polarized distribution-----------------
#function for linear regression
#locfit -------------------------------------------
# locfit will fit the data cell wise and predict values at given cell
heights
#nest the data in stack.all with normalized height to use in locfit for
prediction
locfit.input <- stack.all %>%
group_by(Genotype, Cell.name) %>%
nest()
# locfit.input$data[[1]]
#cell heights chosen to interpolate ecad levels
predictvalues <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#To calculate the y-axis limits in locfit graphs
ylimmax <- max(rawdata$Mean) + 250
ylimmin <- min(rawdata$Mean) - 100
#function to fit locfit
locfit_model <- function(df) {
fit <- locfit(Mean ~ Normalizedheight, data = df)
plot(
fit,
get.data = TRUE,
ylim = c(ylimmin, ylimmax),
xlim = c(0, 1),
xlab = "Normalized Cell Height",
ylab = "Ecad Level(in AU)"
)
par(new = TRUE)
#To superimpose graphs
predictions <- predict(fit, predictvalues)
}
#to plot the locfit graphs
#to start a new graph page and close all old ones,
# so that new graphs don't get superimposed on the old ones
plot.new()
dev.off()
plot.new()
#add locfit data to grouped cells and add cell height(on which prediction
# is made) to each cell
locfit.model <- locfit.input %>%
mutate(fitvalues = data %>% map(locfit_model)) %>%
unnest(fitvalues) %>%
group_by(Genotype, Cell.name) %>%
mutate(Cell_height = predictvalues)
locfit.model <- locfit.model %>% group_by(Genotype, Cell.name)
#to save graph for locfit data
dev.copy(
svg,
file = paste0(filename, "_locfitgraphs.svg"),
width = 4,
height = 4
)
dev.off ()
#for creating groups for each cell for ggplot graph of lines -------------------------------------------------
graph.lines <- stack.all %>% mutate (Uniquecellname = paste(Genotype,
Cell.name))
#to find out the apical perimeter of the cell
apical.perimeter <- graph.lines %>%
group_by(Genotype, Cell.name) %>%
slice(1) %>% select(Genotype, Cell.name, Length)
cell.all <- cell.all %>% full_join(apical.perimeter, by = c("Genotype",
"Cell.name"))
cell.all <- rename(cell.all, "Apical Perimeter" = Length)
# to calculate Genotype wise parameters-------------------------
genotype.peri <- apical.perimeter %>%
group_by(Genotype) %>%
summarise(Mean_peri = mean(Length),
Median_peri = median(Length))
genotype.summary <- cell.all %>% group_by(Genotype) %>%
summarise(
Mean_TotalEcad = mean(Total_Ecad),
Median_TotalEcad = median(Total_Ecad),
Mean_Heightofcell = mean(Height_of_cell),
Median_Heightofcell = median(Height_of_cell),
First_quantile_Height = quantile(Height_of_cell, prob = 0.25),
Third_quantile_Height = quantile(Height_of_cell, prob = 0.75),
count = n()
)
#final dataframe for all calculated parameters by merging the above two data
# frames
genotype.all <- genotype.peri %>% full_join(genotype.summary, by =
"Genotype")
#writing csv--------------------------------------------------------
#write csv with height, normalized height,absolute ecad in each stack
write.csv(
stack.all,
file = paste0(filename, "_stackwise_calculations.csv"),
row.names = FALSE
)
#write csv with cell level parameters
write.csv(
cell.all,
file = paste0(filename, "_cellwise_calculations.csv"),
row.names = FALSE
)
#write csv with locfit fitdata
write.csv(locfit.model,
file = paste0(filename, "_locfitdata.csv"),
row.names = FALSE)
#write csv with genotype level parameters
write.csv(
genotype.all,
file = paste0(filename, "_genotypewise_calculations.csv"),
row.names = FALSE
)
#graphs--------------------------------------------------------------
#to plot boxplots for Total Ecad vs Genotype-------------
p <-
ggplot(data = cell.all, aes(x = Genotype, y = Total_Ecad, fill =
Genotype)) +
theme(
plot.title = element_text(hjust = 0.5),
text = element_text(size = 14),
legend.position = "none",
axis.text.x = element_text(angle = 45, hjust = 1)
) +
labs(title = "Total E-cadherin", x = "Genotype", y = "E-cadherin
Levels(in AU)") +
scale_y_continuous(labels = comma) +
theme(
plot.title = element_text(hjust = 0.5),
legend.position = "none",
text = element_text(size = 14)
) +
geom_boxplot(outlier.shape = NA, alpha = 0.7) +
geom_point(position = position_jitter(width = .1),
shape = 21,
alpha = 0.4)#+
# coord_cartesian(ylim = c(0, 800000))
print(p)
dev.copy(
svg,
file = paste0(filename, "_Total Ecad.svg"),
width = 4,
height = 4
)
dev.off ()
#to plot boxplots for Height of cell vs Genotype--------------
q <-
q <-
ggplot(data = cell.all, aes(x = Genotype, y = Height_of_cell, fill =
Genotype)) +
theme(
plot.title = element_text(hjust = 0.5),
text = element_text(size = 14),
legend.position = "none",
axis.text.x = element_text(angle = 45, hjust = 1)
) +
labs(title = "Height of Cells",
x = "Genotype",
y = expression(paste("Height(in ", mu, "m)"))) +
scale_y_continuous(labels = comma) +
theme(
plot.title = element_text(hjust = 0.5),
legend.position = "none",
text = element_text(size = 14)
) +
geom_boxplot(outlier.shape = NA, alpha = 0.7) +
geom_point(
position = position_jitterdodge(
jitter.height = 0,
jitter.width = 0.5,
dodge.width = .1
),
shape = 21,
alpha = 0.4
)#+
# coord_cartesian(ylim = c(0, 8))
print(q)
dev.copy(
svg,
file = paste0(filename, "_totalcellheight.svg"),
width = 4,
height = 4
)
dev.off ()
#to plot boxplots for Height of cell vs locfit data------------------
r <-
ggplot(data = locfit.model,
aes(
x = factor(Cell_height),
y = fitvalues,
fill = Genotype
)) +
theme(plot.title = element_text(hjust = 0.5), text = element_text(size
=
14)) +
labs(title = "E-cadherin levels at various Cell Heights", x =
"Normalized Height", y = "E-cadherin Levels(in AU)") +
geom_boxplot(position = position_dodge(width = 0.8), outlier.shape =
NA) +
scale_y_continuous(labels = comma) +
geom_point(
alpha = 0.4,
position = position_jitterdodge(jitter.width = .2, dodge.width =
0.8),
shape = 21
)#+
# coord_cartesian(ylim = c(0, 1750))
print(r)
dev.copy(
svg,
file = paste0(filename, "_locfit_boxplot.svg"),
width = 4,
height = 4
)
dev.off ()
#to plot lineplots for cell height vs Ecad levels------------------
s <-
ggplot(
data = graph.lines,
aes(
x = Stackheight,
y = Mean,
group = Uniquecellname,
fill = Genotype,
color = Genotype
)
) +
theme(plot.title = element_text(hjust = 0.5), text = element_text(size
=
14)) +
labs(title = "Height vs Ecad Levels",
x = "Height(in " ~ (mu * m),
y = "E-cadherin Levels(in AU)") +
scale_y_continuous(labels = comma) +
geom_point(position = position_jitter(width = .1),
shape = 21,
alpha = 0.4) +
geom_line(alpha = 0.4)#+
# coord_cartesian(ylim = c(0, 1750))
print(s)
dev.copy(
svg,
file = paste0(filename, "_Rawheight.svg"),
width = 4,
height = 4
)
dev.off ()
#to plot lineplots for Normalized cell height vs Ecad levels----------------
u <-
ggplot(
data = graph.lines,
aes(
x = Normalizedheight,
y = Mean,
group = Uniquecellname,
fill = Genotype,
color = Genotype
)
) +
theme(plot.title = element_text(hjust = 0.5), text = element_text(size
=
14)) +
labs(title = "Normalized Height vs Ecad Levels", x = "Normalized
Height", y = "E-cadherin Levels(in AU)") +
geom_point(shape = 21, alpha = 0.4) +
scale_y_continuous(labels = comma) +
geom_line(alpha = 0.4) +
#coord_cartesian(ylim = c(0, 1750))#+
facet_grid(facets = . ~ Boundary)
print(u)
dev.copy(
svg,
file = paste0(filename, "_Normalized_height.svg"),
width = 4,
height = 4
)
dev.off ()
#to plot lineplots for Normalized cell height vs Normalized Ecad levels-
-----------
t <-
ggplot(
data = graph.lines,
aes(
x = Normalizedheight,
y = NormalizedEcad,
group = Uniquecellname,
fill = Genotype
# color = Genotype
)
) +
theme(
plot.title = element_text(hjust = 0.5),
text = element_text(size = 14),
legend.position = "none",
axis.text.x = element_text(angle = 45, hjust = 1)
) +
scale_y_continuous(labels = comma) +
labs(title = "E-Cadherin Distribution", x = "Normalized Height", y =
"Normalized E-cadherin Levels") +
# stat_smooth(aes(
#x = Normalizedheight,
#y = NormalizedEcad,
#group = Genotype,
#fill = Genotype,
#color = Genotype
# ),method='lm')+
geom_bin2d(bins = 15, alpha = 0.3) +
geom_point(shape = 21, alpha = 0.15) + facet_grid(. ~ Boundary) #+
# stat_density_2d(aes(colour = Genotype), geom="density_2d", alpha= 0.2)
print(t)
dev.copy(
svg,
file = paste0(filename, "_Normalized height Normalized Ecad.svg"),
width = 4,
height = 4
)
dev.off ()
#to plot apical perimeter vs genotype------------
b <-
ggplot(data = apical.perimeter, aes(x = Genotype, y = Length, fill =
Genotype)) +
theme(
plot.title = element_text(hjust = 0.5),
text = element_text(size = 14),
legend.position = "none",
axis.text.x = element_text(angle = 45, hjust = 1)
) +
labs(title = "Apical Perimeter of Cells",
x = "Genotype",
y = "Length(in" ~ (mu * m)) +
geom_boxplot(outlier.shape = NA, alpha = 0.7) +
scale_y_continuous(labels = comma) +
theme(
plot.title = element_text(hjust = 0.5),
legend.position = "none",
text = element_text(size = 14)
) +
geom_point(position = position_jitter(width = .1),
shape = 21,
alpha = 0.4)#+
# coord_cartesian(ylim = c(40, 120))
print(b)
dev.copy(
svg,
file = paste0(filename, "_apical perimeter.svg"),
width = 4,
height = 4
)
dev.off ()
# multiplot graphs-------------------------------------
# To plot multiple graphs in one page
grid.newpage()
pushViewport(viewport(layout = grid.layout(17, 6)))
vplayout <- function(x, y)
viewport(layout.pos.row = x, layout.pos.col = y)
print(r, vp = vplayout(1:6, 2:6))
print(p, vp = vplayout(7:11, 1:2))
print(q, vp = vplayout(7:11, 3:4))
print(b, vp = vplayout(7:11, 5:6))
print(t, vp = vplayout(12:17, 2:5))
dev.copy(
svg,
file = paste0(filename, "multiplot_all.svg"),
width = 10,
height = 10
)
dev.off()
