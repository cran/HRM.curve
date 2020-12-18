difference.curves <- function(RFU_data,
                              derivative = "minus_1st_derivative",
                              normalization = FALSE,
                              standard_diff_curve = 1,
                              number_of_standards = 1,
                              sample_size_standards = 1,
                              sample_number = 1,
                              temp_range_min = 65,
                              temp_range_max = 95,  
                              draw_peaks = FALSE,
                              xlab = "temperature", 
                              ylab = "difference RFU",
                              col_standards = "black",
                              col_samples = "forestgreen",
                              lwd_standards = 1.5,
                              lwd_samples = 0.75,
                              lty = 1,
                              lty_peaks_standards = 1,
                              lty_peaks_samples = 3,
                              xlim = c(temp_range_min, temp_range_max),
                              ...) {
  if(max(sample_number) > (ncol(RFU_data)-1-number_of_standards)){print("ERROR: select sample numbers within your dataset")} else {
    if(number_of_standards == 0 || sample_size_standards == 0){print("ERROR: please select at least 1 standard")} else { 
      if(number_of_standards < sample_size_standards){print("ERROR: total number of standards ('number_of_standards') must be lower or equal the sample size of the standards ('sample_size_standards')")} else {
        if(standard_diff_curve > number_of_standards ||standard_diff_curve >  (number_of_standards/sample_size_standards)){print("ERROR: please select a standard for the difference curve")} else {
          
          melting_curves_1st_neg_deriv <- function(data){
            counter = 1 #counting the loops
            for(i in seq(1, ncol(data), 1)){
              if (i == 1) {
                diff_RFUs <-(-diff(data[,i]))
                assign(paste(names(data)[i],sep=""),diff_RFUs)
                diff_RFUs_1 <- data.frame(get(names(data)[i]))
                colnames(diff_RFUs_1)[i] <- names(data)[i]
              } else {
                diff_RFUs_2 <-(-diff(data[,i]))
                assign(paste(names(data)[i],sep=""),diff_RFUs_2)
                diff_RFUs_1 <- cbind(diff_RFUs_1, diff_RFUs_2)
                colnames(diff_RFUs_1)[{
                  counter = counter + 1
                }] <- names(data)[i]
              }
            }
            diff_RFUs_1 <- as.data.frame(cbind(data[-1,1], diff_RFUs_1[,-1]))
            colnames(diff_RFUs_1) <- c("Temperature", colnames(data[-1]))
            diff_RFUs_1
          }
          
          melting_curves_2nd_deriv <- function(data){
            counter = 1 #counting the loops
            for(i in seq(1, ncol(data), 1)){
              if (i == 1) {
                diff_RFUs <- (diff(diff(data[,i])))
                assign(paste(names(data)[i],sep=""),diff_RFUs)
                diff_RFUs_1 <- data.frame(get(names(data)[i]))
                colnames(diff_RFUs_1)[i] <- names(data)[i]
              } else {
                diff_RFUs_2 <- (diff(diff(data[,i])))
                assign(paste(names(data)[i],sep=""),diff_RFUs_2)
                diff_RFUs_1 <- cbind(diff_RFUs_1, diff_RFUs_2)
                colnames(diff_RFUs_1)[{
                  counter = counter + 1
                }] <- names(data)[i]
              }
            }
            diff_RFUs_1 <- as.data.frame(cbind(data[-c(1:2),1], diff_RFUs_1[,-1]))
            colnames(diff_RFUs_1) <- c("Temperature", colnames(data[-1]))
            diff_RFUs_1
          }
          
          data_1 <- subset(RFU_data, RFU_data[1] >= temp_range_min & RFU_data[1] <= temp_range_max)
          data_range01_data <- apply(data_1[2:ncol(data_1)], MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
          data <- as.data.frame(cbind(data_1[1], data_range01_data))
          
          if(normalization == T){data_1 <- subset(RFU_data, RFU_data[1] >= temp_range_min & RFU_data[1] <= temp_range_max)
          data_range01_data <- apply(data_1[2:ncol(data_1)], MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
          data <- as.data.frame(cbind(data_1[1], data_range01_data))
          } else {if(normalization == FALSE){
            data <- subset(RFU_data, RFU_data[1] >= temp_range_min & RFU_data[1] <= temp_range_max)
          } else{
            print("Choose either 'TRUE' or 'FALSE'; default is 'TRUE'")
          }
          }
          
          if(derivative == "minus_1st_derivative"){
            matrix_melt <- melting_curves_1st_neg_deriv(data)
          } else {
            if(derivative == "2nd_derivative"){
              matrix_melt <- melting_curves_2nd_deriv(data)
            } else {print("Choose either the 'minus_1st_derivative' or '2nd_derivative' derivative; default is 'minus_1st_derivative'")}
          }
          
          counter = 1 #counting the loops
          for(i in seq(2,(number_of_standards+1), sample_size_standards)){
            if (i == 2) {
              data_means_1 <- rowMeans(matrix_melt[i:(i+(sample_size_standards-1))])
              assign(paste(names(matrix_melt)[i],sep=""),data_means_1)
              data_means_2 <- data.frame(get(names(matrix_melt)[i]))
              colnames(data_means_2)[i-1] <- names(matrix_melt)[i]
            } else {
              if(sample_size_standards == 1){
                data_means_1 <- rowMeans(matrix_melt[i])
                assign(paste(names(matrix_melt)[i],sep=""),data_means_1)
                data_means_2 <- cbind(data_means_2, data_means_1)
                colnames(data_means_2)[{
                  counter = counter + 1
                }] <- names(matrix_melt)[i]
              }else{
                data_means_1 <- rowMeans(matrix_melt[i:(i+(sample_size_standards-1))])
                assign(paste(names(matrix_melt)[i],sep=""),data_means_1)
                data_means_2 <- cbind(data_means_2, data_means_1)
                colnames(data_means_2)[{
                  counter = counter + 1
                }] <- names(matrix_melt)[i]
              }
            }
          }
          
          diff_STDs_1 <- ((data_means_2[,standard_diff_curve]) - data_means_2[,1:ncol(data_means_2)])
          diff_STDs_2 <- cbind(matrix_melt[1],diff_STDs_1)
          diff_samples <- ((data_means_2[,standard_diff_curve]) - matrix_melt[,c(sample_number+number_of_standards+1)])
          
          if(sample_number[1] == 0){
            matplot(diff_STDs_2[1], 
                    diff_STDs_2[,2:ncol(diff_STDs_2)], 
                    lty = lty, 
                    type='l', xlab=xlab, ylab=ylab,
                    xlim = xlim,
                    las=1, cex.lab=1.25,lwd=lwd_standards,
                    xaxt = "n", col = col_standards)
          } else {
            matplot(diff_STDs_2[1], 
                    diff_STDs_2[,2:ncol(diff_STDs_2)], 
                    lty = lty,
                    xlim = xlim,
                    type='l', xlab=xlab, ylab=ylab,
                    ylim= c(min(cbind(diff_STDs_2[,-1],diff_samples)*1.05), 
                            max(cbind(diff_STDs_2[,-1],diff_samples))*1.05),
                    las=1, cex.lab=1.25,lwd=lwd_standards,
                    xaxt = "n", col = col_standards)
            
            oldpar <- par(no.readonly = TRUE)
            on.exit(par(oldpar)) 
            par(new=T)
            
            matplot(diff_STDs_2[1], 
                    diff_samples, 
                    lty =  lty,
                    xlim = xlim,
                    type='l', xlab=xlab, ylab=ylab,
                    ylim= c(min(cbind(diff_STDs_2[,-1],diff_samples)*1.05), 
                            max(cbind(diff_STDs_2[,-1],diff_samples))*1.05),
                    las=1, cex.lab=1.25,lwd=lwd_samples,
                    xaxt = "n", col = col_samples)
            
            if(draw_peaks == T){
              #max
              for(i in seq(2, ncol(diff_STDs_2), 1)){
                max_points <-  diff_STDs_2[with(diff_STDs_2, order(-ave(diff_STDs_2[,i], diff_STDs_2[1], FUN = max), -diff_STDs_2[,i])), ]
                if(max_points[1,i] == 0){}else{
                  segments(x0 = max_points[1,1],
                           y0 = max(max_points)*10,
                           x1 = max_points[1,1],
                           y1 = max_points[1,i], lty = lty_peaks_standards,
                           col = rep(col_standards, 10000)[i-1])
                }
              }
              
              diff_samples_temp <- cbind(diff_STDs_2[1], diff_samples)
              
              for(i in seq(2,ncol(diff_samples_temp),1)){
                max_points <-  diff_samples_temp[with(diff_samples_temp, order(-ave(diff_samples_temp[,i], diff_samples_temp[1], FUN = max), -diff_samples_temp[,i])), ]
                segments(x0 = max_points[1,1],
                         y0 = max(max_points)*10,
                         x1 = max_points[1,1],
                         y1 = max_points[1,i], lty = lty_peaks_samples,
                         col = rep(col_samples, 10000)[i-1])
                
              }
              #min
              for(i in seq(2, ncol(diff_STDs_2), 1)){
                max_points <-  diff_STDs_2[with(diff_STDs_2, order(ave(diff_STDs_2[,i], diff_STDs_2[1], FUN = min), diff_STDs_2[,i])), ]
                if(max_points[1,i] == 0){}else{
                  segments(x0 = max_points[1,1],
                           y0 = max(max_points)-1000,
                           x1 = max_points[1,1],
                           y1 = max_points[1,i], lty = lty_peaks_standards,
                           col = rep(col_standards, 10000)[i-1])
                }
              }
              
              diff_samples_temp <- cbind(diff_STDs_2[1], diff_samples)
              
              for(i in seq(2,ncol(diff_samples_temp),1)){
                max_points <-  diff_samples_temp[with(diff_samples_temp, order(ave(diff_samples_temp[,i], diff_samples_temp[1], FUN = min), diff_samples_temp[,i])), ]
                segments(x0 = max_points[1,1],
                         y0 = max(max_points)-1000,
                         x1 = max_points[1,1],
                         y1 = max_points[1,i], lty = lty_peaks_samples,
                         col = rep(col_samples, 10000)[i-1])
              }
            }
          }
          axis(1, at = seq(temp_range_min, temp_range_max, by = 0.1), las=1, tck=-0.02, labels = FALSE)
          axis(1, at = seq(temp_range_min, temp_range_max, by = 1), las=1, tck=-0.03, lwd = 1.25)
          axis(3, at = seq(temp_range_min, temp_range_max, by = 0.1), las=1, tck=-0.02, labels = FALSE)
          axis(3, at = seq(temp_range_min, temp_range_max, by = 1), las=1, tck=-0.03, lwd = 1.25)
        }
      }
    }
  }
}
