library(readxl)

# Read Excel Files
probes <- read_excel("./data/probesample.xlsx")
beta_data <- read_excel("./data/beta.xlsx")

# Creating data matrix
probe_ids <- sub("_.*", "", beta_data$`probe set`)
data.m <- as.matrix(beta_data[, -1])
rownames(data.m) <- probe_ids
# save(data.m, file = "data.m.RData")

# Remove rows with any NA values
complete_rows <- complete.cases(data.m)
data.m <- data.m[complete_rows, ]

# Update probe_ids to match the filtered data.m
probe_ids <- probe_ids[complete_rows]
rownames(data.m) <- probe_ids

print(paste("Removed", sum(!complete_rows), "rows containing NA values"))
print(paste("Final dimensions of data.m:", nrow(data.m), "rows by", ncol(data.m), "columns"))

### DoBMIQ.R
source("./BMIQ_1.4.R")

print("Read Success.")


index <- which(probes$name %in% rownames(data.m))
index <- index[match(rownames(data.m),probes$name[index])]

print("Index created successfully.")

probe_targetid <- probes$targetid[index]
type_ids <- substr(probe_targetid, nchar(probes$targetid) - 1, nchar(probes$targetid) - 1)

###
type1.idx <- which(type_ids==1);
type2.idx <- which(type_ids==2);
design.v <- type_ids


pdf("Profiles.pdf",width=4,height=3);
for(s in 1:ncol(data.m)){
   plot(density(data.m[type1.idx,s]));
   d.o <- density(data.m[type2.idx,s]);
   points(d.o$x,d.o$y,type="l",col="red");
print(s);
}
dev.off();


for(s in 1:ncol(data.m)){
  beta.v <- data.m[,s];
  bmiq.o <- BMIQ(beta.v,design.v,sampleID=s);
  tmp.v <- bmiq.o$nbeta;
  save(tmp.v,file=paste("bmiq",s,".Rd",sep=""));
  print(paste("Done BMIQ for sample ",s,sep=""));
}

bmiq.m <- data.m;
rm(data.m);
for(s in 1:ncol(bmiq.m)){
  load(paste("bmiq",s,".Rd",sep=""));
  bmiq.m[,s] <- tmp.v;
  print(s);
}


save(bmiq.m,file="bmiq.Rd");

source_directory <- "./"
results_directory <- "./results/"

file_types <- "\\.pdf$|\\.Rd$"
list_of_files <- list.files(source_directory, pattern=file_types, full.names=TRUE)

# Move each file to the target directory
for (file in list_of_files) {
  file_name <- basename(file)  # Extract the file name
  target_path <- file.path(results_directory, file_name)  # Create target path
  file.rename(file, target_path)  # Move the file
}

# Confirmation message
cat("Moved", length(list_of_files), "files to", results_directory, "\n")

stop()

### select CpGs (remove CHs) - Wasn't able to do yet...
annoindex <- which(anno450k.m[,1] %in% rownames(bmiq.m))
annoindex <- annoindex[match(rownames(bmiq.m),anno450k.m[annoindex])]
anno.m <- anno450k.m[annoindex,]

cg.idx <- grep("cg",anno.m[,1]);
anno2.m <- anno.m[cg.idx,];
bmiq2.m <- bmiq.m[cg.idx,];
rm(bmiq.m);
save(bmiq2.m,file="BMIQ_RESULT.Rd")
save(anno2.m,file="anno2.m.Rd")
