\name{Denaturation curves}
\alias{denaturation.curves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Denaturation curves for HRM curve analysis
}
\description{
Function draws denaturation curves for selected standards and samples. The number of standards is specified by \link[HRM.curve:denaturation.curves]{number_of_standards}. If all standards have the same number of replicates (balanced design), the mean for each standard can be drawn by specifying \link[HRM.curve:denaturation.curves]{sample_size_standards}; however, \strong{replicates should always be drawn individually first}. The sample numbers of samples for which difference curves will be drawn can be specified as a vector of integers by specifying \link[HRM.curve:denaturation.curves]{sample_number}. The function allows to normalize relative fluorescence values (RFUs) by scaling between 0 and 1 (see \link[HRM.curve:denaturation.curves]{normalization}) within the selected temperature range (\link[HRM.curve:denaturation.curves]{temp_range_min}, and \link[HRM.curve:denaturation.curves]{temp_range_max}).
Please cite Schiwek et. al (2020) <https://doi.org/10.3390/pathogens9040270> if you use this function in your work.
}
\usage{
denaturation.curves(RFU_data,
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
                    xlim,
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
  \item{normalization}{
Logical, RFUs within the selected temperature range (\link[HRM.curve:denaturation.curves]{temp_range_min}, and \link[HRM.curve:denaturation.curves]{temp_range_max}) of each sample are scaled (0,
1). Default is FALSE.
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
Numeric, minimum (\link[HRM.curve:denaturation.curves]{temp_range_min}) and maximum (\link[HRM.curve:denaturation.curves]{temp_range_min}) temperature for analysis. Default \link[HRM.curve:denaturation.curves]{temp_range_min} is 65 and \link[HRM.curve:denaturation.curves]{temp_range_max} is 95.
}
  \item{xlab, ylab, col_samples, col_standards, lwd_standards, lwd_samples, lty, xlim, ...}{
Graphical parameters.
}
}
\details{
\link[HRM.curve]{denaturation.curves} draws denaturation curves for the selected standards and samples.
}
\value{
\link[HRM.curve]{denaturation.curves} returns a plot of denaturation curves.
}
\author{
Lukas Beule, Simon Schiwek, Charlotte Rodemann, Petr Karlovsky
}
\references{
Schiwek, S, Beule, L, Vinas, M, Pfordt, A, von Tiedemann, A, & Karlovsky, P (2020). High-Resolution Melting (HRM) Curve Assay for the Identification of Eight Fusarium Species Causing Ear Rot in Maize. Pathogens, 9(4), 270.

<https://doi.org/10.3390/pathogens9040270>
}
\seealso{
See also \link[HRM.curve]{melting.curves} and \link[HRM.curve]{difference.curves} for further HRM curve analysis.
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

##generate non-normalized denaturation curves
denaturation.curves(RFU_data,
                    normalization = FALSE,
                    number_of_standards = 24,
                    sample_size_standards = 1,
                    sample_number = c(1:10),
                    col_samples = "#000000",
                    col_standards = c(rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), each = 3))
                    )    

##generate normalized denaturation curves from 81 to 92 degrees Celcius
denaturation.curves(RFU_data,
                    normalization = TRUE,
                    number_of_standards = 24,
                    sample_size_standards = 1,
                    sample_number = c(1:10),
                    temp_range_min = 81,
                    temp_range_max = 92,  
                    col_samples = "#000000",
                    col_standards = c(rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), each = 3))
                    )     
}