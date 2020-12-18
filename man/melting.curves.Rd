\name{Melting curves}
\alias{melting.curves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Melting curves for HRM curve analysis
}
\description{
Function draws melting curves for selected standards and samples. The number of standards is specified by \link[HRM.curve:melting.curves]{number_of_standards}. If all standards have the same number of replicates (balanced design), the mean for each standard can be drawn by specifying \link[HRM.curve:melting.curves]{sample_size_standards}; however, \strong{replicates should always be drawn individually first}. The sample numbers of samples for which difference curves will be drawn can be specified as a vector of integers by specifying \link[HRM.curve:melting.curves]{sample_number}. The function allows to normalize relative fluorescence values (RFUs) by scaling between 0 to 1 (see \link[HRM.curve:melting.curves]{normalization}) within the selected temperature range (\link[HRM.curve:melting.curves]{temp_range_min}, and \link[HRM.curve:melting.curves]{temp_range_max}).
Please cite Schiwek et. al (2020) <https://doi.org/10.3390/pathogens9040270> if you use this function in your work.
}
\usage{
melting.curves(RFU_data,
               derivative = "minus_1st_derivative",
               normalization = FALSE,
               number_of_standards = 1,
               sample_size_standards = 1,
               sample_number = 1,
               temp_range_min = 65,
               temp_range_max = 95,
               draw_peaks = FALSE,
               xlab = "temperature", 
               ylab = "-d(RFU)",
               col_samples = "forestgreen",
               col_standards = "black",
               lwd_standards = 1.5,
               lwd_samples = 0.75,
               lty = 1,
               lty_peaks_standards = 1,
               lty_peaks_samples = 3,
               xlim = c(temp_range_min, temp_range_max),
               ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RFU_data}{
Data frame containing RFU data. Data must be numeric and arranged columnwise.
The \strong{first column } must contain the \strong{temperature steps}.
\strong{Standards} must start in the \strong{second column}.
\strong{Samples} must start \strong{following the standards}.
An example of input data can be found here: <https://github.com/LukasBeule/HRM/blob/main/HRM_example_data.txt>
}
  \item{derivative}{
Character, selection of the derivative of the RFU data. For the negative first derivative, enter "minus_1st_derivative”, for the second derivative, enter "2nd_derivative”. Default is "minus_1st_derivative”.
}
  \item{normalization}{
Logical, RFUs within the selected temperature range (\link[HRM.curve:melting.curves]{temp_range_min}, and \link[HRM.curve:melting.curves]{temp_range_max}) of each sample are scaled (0, 1). Default is FALSE.
}
  \item{number_of_standards}{
Integer, total number of standards. Default is 1.
}
  \item{sample_size_standards}{
Integer, if all standards have the same number of replicates (balanced design) and \link[HRM.curve:denaturation.curves]{sample_size_standards} is > 1, the mean of all standards will be plotted instead of each individual standard. Default is 1.
}
  \item{sample_number}{
Vector of integer(s), specifying sample number(s) of samples that will be analyzed. Default is 1.}
  \item{temp_range_min, temp_range_max}{
Numeric, minimum (\link[HRM.curve:melting.curves]{temp_range_min}) and maximum (\link[HRM.curve:melting.curves]{temp_range_min}) temperature for analysis. Default \link[HRM.curve:melting.curves]{temp_range_min} is 65 and \link[HRM.curve:melting.curves]{temp_range_max} is 95.
}
  \item{draw_peaks}{
Logical, draws vertical lines from the maximum (and minimum if \link[HRM.curve:melting.curves]{derivative} is "2nd_derivative"") of the melting curves to the x-axis. Default is FALSE.
}
  \item{xlab, ylab, col_standards, col_samples, lwd_standards, lwd_samples, lty, lty_peaks_standards, lty_peaks_samples, xlim, ...}{
Graphical parameters.
}
}
\details{
\link[HRM.curve]{melting.curves} draws melting curves for selected standards and samples.
}
\value{
\link[HRM.curve]{melting.curves} returns a plot of difference curves.
}
\author{
Lukas Beule, Simon Schiwek, Charlotte Rodemann, Petr Karlovsky
}
\references{
Schiwek, S, Beule, L, Vinas, M, Pfordt, A, von Tiedemann, A, & Karlovsky, P (2020). High-Resolution Melting (HRM) Curve Assay for the Identification of Eight Fusarium Species Causing Ear Rot in Maize. Pathogens, 9(4), 270.

<https://doi.org/10.3390/pathogens9040270>
}
\seealso{
See also \link[HRM.curve]{denaturation.curves} and \link[HRM.curve]{difference.curves} for further HRM curve analysis.
}
\examples{
##Input data must be numeric and arranged columnwise.
##Input data should not contain any categorial data.
##An example dataset containing input data with 8 standards
##with 3 technical replicates each and 10 samples can be found here:
##<https://github.com/LukasBeule/HRM/blob/main/HRM_example_data.txt>

##load example RFU data:
RFU_data <- read.table("https://raw.githubusercontent.com/LukasBeule/HRM/main/HRM_example_data.txt",
                       sep="\t", dec=".", header=TRUE)

##generate non-normalized melting curves from 81 to 92 degrees Celcius 
##using negative first derivative
melting.curves(RFU_data,
               derivative = "minus_1st_derivative",
               normalization = FALSE,
               number_of_standards = 24,
               sample_size_standards = 1,
               sample_number = c(1:10),
               temp_range_min = 81,
               temp_range_max = 92,
               col_samples = "#000000",
               col_standards = c(rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), each = 3))
               )
                                      
##generate normalized melting curves from 81 to 92 degrees Celcius 
##using negative first derivative and draw_peak = TRUE
melting.curves(RFU_data,
               derivative = "minus_1st_derivative",
               normalization = TRUE,
               number_of_standards = 24,
               sample_size_standards = 1,
               sample_number = c(1:10),
               temp_range_min = 81,
               temp_range_max = 92,
               draw_peaks = TRUE,
               col_samples = "#000000",
               col_standards = c(rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), each = 3))
               )
                                      
##generate non-normalized melting curves from 81 to 92 degrees Celcius 
##using second derivative and draw_peak = TRUE
melting.curves(RFU_data,
               derivative = "2nd_derivative",
               normalization = FALSE,
               number_of_standards = 24,
               sample_size_standards = 1,
               sample_number = c(1:10),
               temp_range_min = 81,
               temp_range_max = 92,
               draw_peaks = TRUE,
               ylab = "d2(RFU)",
               col_samples = "#000000",
               col_standards = c(rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), each = 3))
               )
}