denaturation.curves <- function(RFU_data,
                                normalization = FALSE,
                                number_of_standards = 1,
                                sample_size_standards = 1,
                                sample_number = 1,
                                temp_range_min = 65,
                                temp_range_max = 95,  
                                xlab = "temperature", 
                                ylab = "RFU",
                                col_standards = "black",
                                col_samples = "forestgreen",
                                lwd_standards = 1.5,
                                lwd_samples = 0.75,
                                lty = 1,
                                xlim = c(temp_range_min, temp_range_max),
                                ...) {
  
  if(max(sample_number) > (ncol(RFU_data)-1-number_of_standards)){print("ERROR: select sample numbers within your dataset")} else {
    if(number_of_standards == 0 || sample_size_standards == 0){print("ERROR: please select at least 1 standard")} else { 
      if(number_of_standards < sample_size_standards){print("ERROR: total number of standards ('number_of_standards') must be lower or equal the sample size of the standards ('sample_size_standards')")} else {
        
        data_1 <- subset(RFU_data, RFU_data[1] >= temp_range_min & RFU_data[1] <= temp_range_max)
        data_range01_data <- apply(data_1[2:ncol(data_1)], MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
        data <- as.data.frame(cbind(data_1[1], data_range01_data))
        
        if(normalization == T){data_1 <- subset(RFU_data, RFU_data[1] >= temp_range_min & RFU_data[1] <= temp_range_max)
        data_range01_data <- apply(data_1[2:ncol(data_1)], MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
        data <- as.data.frame(cbind(data_1[1], data_range01_data))
        } else {if(normalization == FALSE){
          data <- subset(RFU_data, RFU_data[1] >= temp_range_min & RFU_data[1] <= temp_range_max)
        } else{
          print("ERROR: Choose either 'TRUE' or 'FALSE'; default is 'TRUE'")
        }
        }
        counter = 1 #counting the loops
        for(i in seq(2,(number_of_standards+1), sample_size_standards)){
          if (i == 2) {
            data_means_1 <- rowMeans(data[i:(i+(sample_size_standards-1))])
            assign(paste(names(data)[i],sep=""),data_means_1)
            data_means_2 <- data.frame(get(names(data)[i]))
            colnames(data_means_2)[i-1] <- names(data)[i]
          } else {
            data_means_1 <- rowMeans(data[i:(i+(sample_size_standards-1))])
            assign(paste(names(data)[i],sep=""),data_means_1)
            data_means_2 <- cbind(data_means_2, data_means_1)
            colnames(data_means_2)[{
              counter = counter + 1
            }] <- names(data)[i]
          }
        }
        
        data_means_2 <- cbind(data[1],data_means_2)
        
        if(sample_number[1] == 0){
          matplot(data_means_2[1], 
                  data_means_2[,2:ncol(data_means_2)], 
                  lty = lty, 
                  type='l', xlab=xlab, ylab=ylab,
                  xlim = xlim,
                  las=1, cex.lab=1.25,lwd=lwd_standards,
                  xaxt = "n", col = col_standards)
        } else {
          matplot(data_means_2[1], 
                  data_means_2[,2:ncol(data_means_2)], 
                  lty = lty, 
                  type='l', xlab=xlab, ylab=ylab,
                  xlim = xlim,
                  las=1, cex.lab=1.25,lwd=lwd_standards,
                  ylim= c(min(cbind(data_means_2[,2:ncol(data_means_2)],data[,c(sample_number+number_of_standards+1)])*0.99), 
                          max(cbind(data_means_2[,2:ncol(data_means_2)],data[,c(sample_number+number_of_standards+1)])*1.01)),
                  xaxt = "n", col = col_standards)
          
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar)) 
          par(new=T)
          
          matplot(data[1], 
                  data[,c(sample_number+number_of_standards+1)], 
                  lty =  lty, 
                  type='l', xlab=xlab, ylab=ylab,
                  xlim = xlim,
                  las=1, cex.lab=1.25,lwd=lwd_samples,
                  ylim= c(min(cbind(data_means_2[,2:ncol(data_means_2)],data[,c(sample_number+number_of_standards+1)])*0.99), 
                          max(cbind(data_means_2[,2:ncol(data_means_2)],data[,c(sample_number+number_of_standards+1)])*1.01)),
                  
                  xaxt = "n", col = col_samples)
        }
        axis(1, at = seq(temp_range_min, temp_range_max, by = 0.1), las=1, tck=-0.02, labels = FALSE)
        axis(1, at = seq(temp_range_min, temp_range_max, by = 1), las=1, tck=-0.03, lwd = 1.25)
        axis(3, at = seq(temp_range_min, temp_range_max, by = 0.1), las=1, tck=-0.02, labels = FALSE)
        axis(3, at = seq(temp_range_min, temp_range_max, by = 1), las=1, tck=-0.03, lwd = 1.25)
      }
    }
  }
}
