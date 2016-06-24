#!/usr/bin/Rscript

# Copyright (C) 2013,2014 Ole Tange, Mike DeGiorgio, Anna-Sapfo
# Malaspinas, Jose Victor Moreno-Mayar, Yong Wang and Free Software
# Foundation, Inc.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# For Interative use:
#    eval(parse(text=readLines("dist/src/bammds_plot.r")))
# sort( sapply(ls(),function(x){object.size(get(x))}))

legend_file <- "tmp/Han_700000reads_hg19,Wollstein.hg19,.legend.filled.csv";
asd_file <- "tmp/Han_700000reads_hg19,Wollstein.hg19,.asd";
mds_file <- "/tmp/mike.pdf"
xvector <- 1;
yvector <- 2;
option_mds <- T;
option_pca <- T;
option_summary <- T;
option_compact <- F;
PlotTitle <- "Worldwide" #used only if option_compact

arg <- commandArgs(trailingOnly = TRUE);
asd_file <- arg[1];
legend_orig <- arg[2];
legend_file <- arg[3];
mds_file <- arg[4];
xvector <- as.integer(arg[5]);
yvector <- as.integer(arg[6]);
option_mds <- as.logical(arg[7]);
option_pca <- as.logical(arg[8]);
option_summary <- as.logical(arg[9]);
Closest<-NA

if(is.na(xvector)) {
  message("bammds_plot.r: Warning: xvector undefined. Using 1.")
  xvector <- 1;
}

if(is.na(yvector)) {
  message("bammds_plot.r: Warning: yvector undefined. Using 2.")
  yvector <- 2;
}

## Read in data
read_data_files <- function(asd_file,legend_file) {
  ## 3.9s:
  ##  template_line <- read.table(nrow=1,file=asd_file);
  ##  asd_orig <<- scan(file=asd_file,what=template_line);
  ## slow:
  ## f <- file(asd_file_unc)
  ## system.time(bigdf <- sqldf("select * from f", dbname = tempfile(), file.format = list(sep=" ",header = F, row.names = F)))
  ## slowest:
  ## asd_orig <<- read.table(asd_file);

  if(option_summary) {
    nrows <- -1L;
  } else {
    ## Read a single row to find the number of rows to read if --no-summary
    first_line <- read.table(nrow=1,file=asd_file);
    nrows <- length(first_line)/2;
  }
  require_install("data.table");
  asd_orig <- as.data.frame(data.table::fread(asd_file, sep="auto", sep2="auto", nrows=nrows,
                                              header="auto", na.strings="NA",
                                              stringsAsFactors=FALSE, verbose=FALSE, autostart=30L,
                                              skip=-1L, select=NULL, colClasses=NULL,
                                              integer64=getOption("datatable.integer64")))
  ## Remove the last column if uneven (It is an artifact)
  asd_orig <<- asd_orig[,1:(ncol(asd_orig)-ncol(asd_orig)%%2)]

  legends <<- read.csv(legend_file);
  ## Prepend # to the hex color if it is missing.
  c <- as.character(legends$Color);
  c[substr(c,1,1)!="#"] <- paste("#",c[substr(c,1,1)!="#"],sep="");
  legends$Color <<- c;
}

## Identify individuals removed in legend file.
identify_individuals_removed_in_legend_file <- function() {
  ## Remember number of unfiltered individuals
  ## width/2 because there are 2 values per individual
  num_indv_unfiltered <<- dim(asd_orig)[2]/2;  
  ## Remove population and empty lines
  l <- legends[ ! legends$Individual %in% c("*","") ,];
  ## Sort them by the order column
  legend_sort <<- l[order(l$order),];
  ## Identify the lines with sample = -1
  keep_indv <<- legend_sort$sample != -1;
}

## Unpack data from asd_orig
unpack_data_from_asd_orig <- function() {
  ## asd_orig contains a row for each individual containing:
  ##   (asd, count) for all individuals,
  ##   For every marker: 1 if defined for this individual
  n <- num_indv_unfiltered;
  asds_counts <- asd_orig[1:n,];
  # all_markers <<- as.big.matrix(t(data.matrix(asd_orig[(n+1):nrow(asd_orig),]))[c(T,F),],type="char");
  if(option_summary) {
    all_markers <<- t(data.matrix(asd_orig[(n+1):nrow(asd_orig),]))[c(T,F),];
  }
  ## asds_counts is asd[1], count[1], asd[2], count[2], ...
  all_asds <- asds_counts[c(T,F)];
  all_counts <<- asds_counts[c(F,T)];
  ## Distance per defined snp
  all_asd <<- all_asds/all_counts;
  ## Empty asd_orig to save memory
  asd_orig <<- "Emptied";
  all_asds <- "Emptied";
  asds_counts <- "Emptied";
  gc();
}

## Identify individuals with no shared markers
identify_individuals_with_no_shared_markers <- function() {
  n <- num_indv_unfiltered;
  ## Every individual must share some data with
  ## all other individuals.
  ## If it does not, mark the individual for removal
  
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      if(all_counts[i,j] <= 0) {
        if(keep_indv[i]) {
	  message(paste(sep="",
                        "Individual ",i, "(", legend_sort[i,c("Individual")], ")",
                        " shares no markers with ",
                        j, "(", legend_sort[j,c("Individual")], "). Removing ", i));
	  keep_indv[i] <<- F;
        }
      }
    }
  }
}

## Remove individuals
remove_individuals <- function() {
  ## Remove marked individuals from variables
  legend_indv <<- legend_sort[keep_indv,];
  asd <<- all_asd[keep_indv,keep_indv];
  if(option_summary) {
    ## markers <<- as.big.matrix(all_markers[keep_indv,],type="char");
    markers <<- all_markers[keep_indv,];
    all_markers <<- "Emptied";
    gc();
  }
}

## Population legends is all legends except individual and empty lines
compute_legend_pop <- function() {
  l_pop <- legends[ legends$Individual %in% c("*",""), ]
  l_pop2 <- l_pop[ ! l_pop$Population %in% c("*","") ,];
  legend_pop <<- l_pop2[ l_pop2$sample %in% c(0,1) ,];
}

compute_sample_reference_list <- function() {
  ## Which individuals are marked as samples (BAM-files)?
  samplelist <<- legend_indv$sample == 1;
  ## Which populations are marked as samples (BAM-files)?
  samplepoplist <<- legend_pop$sample == 1;
  ## Which individuals are references (tped-files)?
  referencelist <<- legend_indv$sample == 0;
  ## Which populations are references (tped-files)?
  referencepoplist <<- legend_pop$sample == 1;
}

perform_scaling <- function() {
  ## Perform the scaling
  ##
  ## k is the number of dimensions. If you have n points, you need n-1
  ## dimensions to set them. eig indicates that we want the
  ## eigenvalues. add indicates that we want to add a constant to all the
  ## distances so the eigenvalues after the transformation are going to
  ## be positive.

  if(length(asd) < 3 || dim(asd)[1] < 3) {
    message("Error: Fewer than 3 individuals. Cannot plot.");
    quit("no",1);
  }
  classical.mds <- cmdscale(d = asd, k = (nrow(asd)-1),eig=TRUE, add=T); 
  
  ## Extract vectors
  mds.x <<- classical.mds$points[,xvector]; #coords for first dimension picked
  mds.y <<- classical.mds$points[,yvector]; #coords for second dimension picked
  r <<- rank(mds.x); #rank of the matrix...independent eigenvectors
  
  lambdas.norm <<- 100*classical.mds$eig/sum(classical.mds$eig); # get the % explained by each dimension: TO CHECK

  pca <- prcomp(asd);
  ## Extract vectors
  pca.points <- pca$x;
  pca.x <<- pca.points[,xvector]
  pca.y <<- pca.points[,yvector]
}

## Convert "B" => ascii(B) = 66
Letters2Nums <- function(x) {
	y <- x
	for(i in seq_along(x)) {
		if(identical(grep("[0-9]", x[i]), integer(0))) {
			y[i] <- as.numeric(charToRaw(x[i]))
		} else {
			y[i] <- as.numeric(x[i])
		}
	}
	return(as.numeric(y))
}

## Determine the output file format and open accordingly
open_plot_file <- function(mds_file) {
  if (option_compact){ width_inches <- 3.34;
                   height_inches<- 3.34
  }
  else{width_inches <- 6.26;
  height_inches <- 3.5;
  }
  resolution <- 300
  if(grepl(".pdf$", mds_file, ignore.case = T)) {
    pdf(mds_file, height=height_inches, width=width_inches,useDingbats=F);
  } else if(grepl(".png$", mds_file, ignore.case = T)) {
    save_file <- gsub("(....)$", ".%1d\\1", mds_file)
    png(save_file, height=height_inches, width= width_inches,units = 'in',res=resolution);
    mds_file <<- gsub("(....)$", ".*\\1", mds_file)
  } else if(grepl(".svg$", mds_file, ignore.case = T)) {
    save_file <- gsub("(....)$", ".%1d\\1", mds_file)
    svg(save_file, height=height_inches, width=width_inches);
    mds_file <<- gsub("(....)$", ".*\\1", mds_file)
  } else if(grepl(".jpg$", mds_file, ignore.case = T)) {
    save_file <- gsub("(....)$", ".%1d\\1", mds_file)
    jpeg(save_file, height=height_inches, width= width_inches,units = 'in',res=resolution);
    mds_file <<- gsub("(....)$", ".*\\1", mds_file)
  } else if(grepl(".jpeg$", mds_file, ignore.case = T)) {
    save_file <- gsub("(.....)$", ".%1d\\1", mds_file)
    jpeg(save_file, height=height_inches, width= width_inches,units = 'in',res=resolution);
    mds_file <<- gsub("(.....)$", ".*\\1", mds_file)
  } else if(grepl(".tif$", mds_file, ignore.case = T)) {
    save_file <- gsub("(....)$", ".%1d\\1", mds_file)
    tiff(save_file, height=height_inches, width= width_inches,units = 'in',res=resolution);
    mds_file <<- gsub("(....)$", ".*\\1", mds_file)
  } else if(grepl(".tiff$", mds_file, ignore.case = T)) {
    save_file <- gsub("(.....)$", ".%1d\\1", mds_file)
    tiff(save_file, height=height_inches, width= width_inches,units = 'in',res=resolution);
    mds_file <<- gsub("(.....)$", ".*\\1", mds_file)
  } else if(grepl(".csv$", mds_file, ignore.case = T)) {
    df <- data.frame(legend_indv[,c("Population", "pop_label",
                                    "Individual", "indv_label")],
                     mds.x, mds.y, pca.x, pca.y,
                     legend_indv[,"Color"],
                     legend_indv[,"pch"]);
    names(df) <- c("Population", "pop_label", "Individual",
                   "indv_label", "mds_x", "mds_y", "pca_x", "pca_y",
                   "Color","Symbol");
    write.table(df,file=mds_file,sep=",",row.names=F)
  } else {
    ## Unknown format
    error <- paste("Unknown plot format:",mds_file);
    write(error, stderr());
    quit("no",1);
  }
}

close_plot <- function(mds_file) {
  write(paste(sep="", "Plot located in: ", mds_file),stdout());
  dummy <- dev.off();
}

## Plot populations
plot_populations <- function(pca_mds,x,y) {
  ## 1 column if compact option: The MDS-plot and the legend
  if (option_compact){
  par(mfrow=c(1,1),mar=c(1, 1, 1.5, 0.5),mgp=c(0,0,0));
  }
  else{
  ## Split the paper into 2 columns: The MDS-plot and the legend
  par(mfrow=c(1,2),mar=c(1, 1, 1.5, 0.5),mgp=c(0,0,0));
  }
  
  indv_cols <- c("Population", "sample","Individual","order");
  pop_cols <- c("Population", "pop_label", "Color", "pch", "cex");
  
  # Get the pch and the color from population and the rest from the individual.
  indv_in_pop_unsorted <- merge(legend_indv[,indv_cols],legend_pop[,pop_cols],by="Population",);
  indv_in_pop <- indv_in_pop_unsorted[order(indv_in_pop_unsorted$order),];
                   
  ## Plot MDS for populations
  {
    label_pop <- as.character(legend_pop[samplepoplist,"pop_label"]);
    if(!option_compact){PlotTitle <- paste(pca_mds,collapse=' ',paste(label_pop,collapse=" "));}
        
    ## Set the limits and title
    plot(c(),c(), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)),
         xaxt="n",yaxt="n",
         xlab=paste("dimension ",xvector,"(",round(lambdas.norm[xvector],2),"%)"),
         ylab=paste("dimension ",yvector,"(",round(lambdas.norm[yvector],2),"%)"),
         main=PlotTitle,
         cex.main=0.8,cex.lab=0.7);
    cex <- legend_indv$cex;
    pch <- as.character(indv_in_pop$pch);
    ## Indexes of [n] and pch points.
    circle_points <- grep("\\[.*\\]", indv_in_pop$pch);

if(length(circle_points)==1){
	getDist<-function(centroid, x1, y1){
		pair<-rbind(centroid[1:2], c(x1,y1))
		return(as.numeric(dist(pair)))
	}
	
tablita<-data.frame(x=mds.x[-1], y=mds.y[-1], p=as.character(indv_in_pop$Population)[-1])
	centroids<-data.frame(x=tapply(tablita$x, tablita$p, mean), y=tapply(tablita$y, tablita$p, mean), p=names(tapply(tablita$x, tablita$p, mean)))
	x1<-mds.x[1]
	y1<-mds.y[1]
	Closest<<-names(which.min(apply(centroids, 1, getDist, x1, y1)))
	#message(paste("Closest population:", Closest))
}


    pch_points <- grep("\\[.*\\]", indv_in_pop$pch, invert=T);
    ## Replace the [1] .. [n] with a dummy space in pch. pch is not used for plotting [1]
    pch[circle_points] <- " "
    ## Convert pch single letters to the corresponding pch value (a => 97, B => 66).
    pch <- Letters2Nums(pch);
    colors <- as.character(indv_in_pop[,"Color"]);
    
    ## Plot point with pch != [n]: Use the pch as point
    if(length(pch_points) != 0) {
      points(x[pch_points], y[pch_points], cex=cex[pch_points], pch=pch[pch_points], col=colors[pch_points]);
    }
    
    ## Plot points with pch == [n]: Text with a circle around.
    if(length(circle_points) != 0) {
      ## Circle_text of [1] = 1
      circle_text <- gsub("\\[|\\]","",indv_in_pop$pch);
      if (Letters2Nums(circle_text[1]) > 25){ #if the first of the ancient points is not a number between 1 and 25
      text(x[circle_points], y[circle_points], circle_text[circle_points], col=rgb(0,0,0), cex=cex[circle_points],
           xlim=c(min(x),max(x)), ylim=c(min(y),max(y)));
      points(x[circle_points], y[circle_points], pch=1, col=colors[circle_points], cex=cex[circle_points]+2);
      }else{ # otherwise we assume that none of [] numbers should be circled...(this should be improved)
          points(x[circle_points], y[circle_points], pch=Letters2Nums(circle_text), col=colors[circle_points], cex=cex[circle_points]);
      }
          
    }
  }

  ## Remove duplicate pop_labels (useful for plotting continents)
  ## Check duplicates on pop_label,Color,pch
  ## legends <- legend_pop[!duplicated(legend_pop[,c('pop_label','Color','pch')]),]
  ## Check duplicates on pop_label only
  legends <- legend_pop[!duplicated(legend_pop[,c('pop_label')]),]
  
  ## Plot legend for populations
  plot_legend(legends$pop_label,legends$Color,legends$pch);
}

## Plot individuals
plot_individuals <- function(pca_mds,x,y) {
  ## 1 column if compact option: The MDS-plot and the legend
  if (option_compact){
  par(mfrow=c(1,1),mar=c(1, 1, 1.5, 0.5),mgp=c(0,0,0));
  }
  else{
  ## Split the paper into 2 columns: The MDS-plot and the legend
  par(mfrow=c(1,2),mar=c(1, 1, 1.5, 0.5),mgp=c(0,0,0));
  } 
  ## Plot MDS for individuals
  
    labelindivs <- as.character(legend_indv[samplelist,"indv_label"]);
  
    PlotTitle <- paste(pca_mds,labelindivs,collapse=' ');
    ## Set the limits and title
    plot(c(),c(), xlim=c(min(x),max(x)), ylim=c(min(y),max(y)),
         xaxt="n", yaxt="n", 
         xlab=paste("dimension ",xvector,"(",round(lambdas.norm[xvector],2),"%)"),
         ylab=paste("dimension ",yvector,"(",round(lambdas.norm[yvector],2),"%)"),
         main=PlotTitle,
         cex.main=0.8,cex.lab=0.7);
    cex <- legend_indv$cex;
    pch <- as.character(legend_indv$pch);
    ## Indexes of [n] and pch points.
    circle_points <- grep("\\[.*\\]", legend_indv$pch);
    pch_points <- grep("\\[.*\\]", legend_indv$pch, invert=T);
    ## Replace the [1] .. [n] with a dummy space in pch. pch is not used for plotting [1]
    pch[circle_points] <- " "
    ## Convert pch single letters to the corresponding pch value (a => 97, B => 66).
    pch <- Letters2Nums(pch);
    colors <- as.character(legend_indv[,"Color"]);
    
    ## Plot point with pch != [n]: Use the pch as point
    if(length(pch_points) != 0) {
      points(x[pch_points], y[pch_points], cex=cex[pch_points], pch=pch[pch_points], col=colors[pch_points]);
    }
    
    ## Plot points with pch == [n]: Text with a circle around.
    if(length(circle_points) != 0) {
      ## Circle_text of [1] = 1
      circle_text <- gsub("\\[|\\]","",legend_indv$pch);
      if (Letters2Nums(circle_text[1]) > 25){ #if the first of the ancient points is not a number between 1 and 25
      text(x[circle_points], y[circle_points], circle_text[circle_points], col=rgb(0,0,0), cex=cex[circle_points],
           xlim=c(min(x),max(x)), ylim=c(min(y),max(y)));
      points(x[circle_points], y[circle_points], pch=1, col=colors[circle_points], cex=cex[circle_points]+2);
      }else{ # otherwise we assume that none of [] numbers should be circled...(this should be improved)
          points(x[circle_points], y[circle_points], pch=Letters2Nums(circle_text), col=colors[circle_points], cex=cex[circle_points]);
      }
      
    ### Plot points with pch == [n]: Text with a circle around.
    #if(length(circle_points) != 0) {
    #  ## Circle_text of [1] = 1
    #  circle_text <- gsub("\\[|\\]","",legend_indv$pch);
    #  text(x[circle_points], y[circle_points], circle_text[circle_points], col=rgb(0,0,0), cex=cex[circle_points],
    #       xlim=c(min(x),max(x)), ylim=c(min(y),max(y)));
    #  points(x[circle_points], y[circle_points], pch=1, col=colors[circle_points], cex=cex[circle_points]+2);
    #}
  }
  
  ## Plot legend for individuals
  plot_legend(legend_indv$indv_label,colors,legend_indv$pch);
}

plot_legend <- function(l,colors,pch) {
  ## The max number of Ms that will fit in a N column legend
  ## 8 col = MM = 2
  ## 7 col = MMM = 3
  ## 6 col = MMMM = 4
  ## 5 col = MMMMMM = 6
  ## 4 col = MMMMMMMM = 8
  ## 3 col = MMMMMMMMMMM = 11
  ## 2 col = MMMMMMMMMMMMMMMMMM = 18
  number_of_ms <- c(99, 18, 11, 8, 6, 4, 3, 2);
                    
  label <- as.character(l);
  ## In the legend we want [n] to be shown as n. Letters2Nums(pch) 
  pch <- gsub("\\[|\\]", "", pch);
  pch <- Letters2Nums(pch);

  ## 19 legends per column
  rows = 22

  ## Choose number of columns
  n = 0;
  while(n < length(label)) {
    for(columns in 8:1) {
      chunk_size <- rows*columns;
      start <- 1 + n;
      end <- min(length(label),start + chunk_size - 1);
      label_chunk <- label[start:end];
      max_label_len <- max(as.numeric(lapply(label_chunk, function(x) { nchar(as.character(x)) })))
      if(max_label_len < number_of_ms[columns]) {
        ## Do the plot
        colors_chunk <- as.character(colors[start:end]);
        pch_chunk <- pch[start:end];       
        
        if(option_compact){       
            legend("center", legend=label_chunk, col=colors_chunk, pch=pch_chunk, cex=0.65, ncol=ceiling(length(label_chunk)/rows),bty="n");
        }
        else{
            plot.new(); 
            plot.window(xlim = c(0, 5), ylim = c(0, 2));
            title(main = "");
            legend("topleft", legend=label_chunk, col=colors_chunk, pch=pch_chunk, cex=0.65, ncol=ceiling(length(label_chunk)/rows));    
        }
        ## Continue
        n = n + chunk_size;
        ## Next n
        break;
      }
    }
  }
}
  
# Compute summary information
summary_information <- function(asd_file, legend_file, mds_file) { 
  sample_names <- c(as.character(legend_indv[samplelist,"indv_label"]), "Reference panel");
  num_markers <- dim(markers)[2];
  
  ## Remove markers that are not shared with any sample
  ## To save memory (does not work for single sample, and
  ## memory effect seems minimal).
  ## markers <<- markers[,colSums(markers[samplelist,]) >= 1]

  sample_markers <- markers[samplelist,];
  # The values are 0 or 1 so the merged is simply the max of each row.
  ref_merged <- apply(markers[referencelist,], 2, max);
  # size: 3 x 24000
  sample_ref_markers <- rbind(sample_markers,ref_merged);

  overlapping_markers <- sample_ref_markers %*% t(sample_ref_markers);
  last <- dim(overlapping_markers)
  overlapping_markers[last[1],last[2]] <- num_markers;
  rownames(overlapping_markers) <- sample_names;
  colnames(overlapping_markers) <- sample_names;

if(is.na(Closest)){
  write(paste(sep="", "SUMMARY INFORMATION\n###################\n",
              "Number of markers: ", num_markers, "\n",
              "Overlapping markers between the samples and the reference panel"),stdout());
}else{
  write(paste(sep="", "SUMMARY INFORMATION\n###################\nClosest population: ", Closest, "\n", 
              "Number of markers: ", num_markers, "\n",
              "Overlapping markers between the samples and the reference panel"),stdout());
}
  print(overlapping_markers);

  require_install("gridExtra");
  par(mfrow=c(1,1));
  plot.new();
  title(main = "Summary information")
  grid.table(overlapping_markers,
             gpar.coretext=gpar(fontsize = 8),
             gpar.coltext = gpar(fontsize = 8),
             gpar.rowtext = gpar(fontsize = 8)
             );
  grid.text(paste(sep="", "Date: ", date(), "\n",
                  "ASD-file:\n", asd_file, "\n",
                  "Legend-file:\n", legend_orig, "\n",
                  "Plot in:\n", mds_file, "\n"),0.5,0.75,gp=gpar(fontsize=5));
if(is.na(Closest)){
  grid.text(paste(sep="", "Number of markers: ", num_markers, "\n",
                  "If you use this tool for a publication please cite \nMalaspinas A-S, Tange O, Moreno-Mayar JV, Rasmussen M, DeGiorgio M, Wang Y, et al.\nbammds: a tool for assessing the ancestry of low-depth whole-genome data\nusing multidimensional scaling (MDS). Bioinformatics. 2014 Jun 28;btu410.\n"),0.5,0.2,
            gp=gpar(fontsize=8));
}else{
  grid.text(paste(sep="", "Number of markers: ", num_markers, "\nClosest population in the first plot (page 1): ", Closest, "\n", 
                  "If you use this tool for a publication please cite\nMalaspinas A-S, Tange O, Moreno-Mayar JV, Rasmussen M, DeGiorgio M, Wang Y, et al.\nbammds: a tool for assessing the ancestry of low-depth whole-genome data\nusing multidimensional scaling (MDS). Bioinformatics. 2014 Jun 28;btu410.\n"),0.5,0.2,
            gp=gpar(fontsize=8));

}



}

require_install <- function(pkg) {
  local({r <- getOption("repos");
         r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})
  list.of.packages <- c(pkg)
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  dir <- unlist(strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep))[1L]
  dir.create(dir, recursive = TRUE, showWarnings = FALSE);
  .libPaths( c( .libPaths(), dir) );
  if(length(new.packages)) install.packages(new.packages, dir);
  require(pkg, character.only=T, quiet=TRUE);
}

all <- function() {
  read_data_files(asd_file, legend_file);
  identify_individuals_removed_in_legend_file();
  unpack_data_from_asd_orig();
  identify_individuals_with_no_shared_markers();
  remove_individuals();
  compute_legend_pop();
  compute_sample_reference_list();
  perform_scaling();
  open_plot_file(mds_file);
  if(option_mds) {
    plot_populations("MDS:",mds.x,mds.y);
  }
  if(option_pca) {
    plot_populations("PCA:",pca.x,pca.y);
  }
  if(option_mds&!option_compact) {
    plot_individuals("MDS:",mds.x,mds.y);
  }
  if(option_pca&!option_compact) {
    plot_individuals("PCA:",pca.x,pca.y);
  }
  if(option_summary) {
    summary_information(asd_file, legend_file, mds_file);
  }
  close_plot(mds_file);
}

all();


