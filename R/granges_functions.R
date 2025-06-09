### Inter GRanges overlaps and merging overlapping regions

#' Check if there are any overlapping regions within a GRanges object
#'
#' @param gr A GRanges object.
#' @param ignore.strand A logical value indicating whether to ignore strand when checking for overlaps.
#' Default is TRUE.
#' @return FALSE if there are no overlapping regions in gr and TRUE if there are.
#' @export
check_self_overlaps = function(gr, ignore.strand = TRUE){
  
  # Find which regions in gr overlap
  overlaps = data.frame(findOverlaps(gr, gr, ignore.strand = ignore.strand))
  
  # Check if any region overlaps any other region besides itself
  overlaps = dplyr::filter(overlaps, queryHits != subjectHits)
  
  # Return TRUE or FALSE depending on if any overlaps found
  return(nrow(overlaps) > 0)

} 

#' Merge overlapping regions in a GRanges object
#'
#' @param gr A GRanges object.
#' @return A GRanges object.
#' @export
merge_overlapping_regions = function(gr){
  
  # Merge overlapping regions in gr and return
  return(reduce(gr, min.gapwidth = 0))

} 

#' Merge touching and overlapping regions in a GRanges object
#'
#' @param gr A GRanges object.
#' @return A GRanges object.
#' @export
merge_touching_regions = function(gr){
  
  # Merge overlapping regions in gr and return
  return(reduce(gr, min.gapwidth = 1))

} 

### Number of bases covered by a GRanges and the size of the intersection between two GRanges

#' Calculate the number of unique bases covered by all regions in a GRanges object
#'
#' @param gr A GRanges object.
#' @return An numeric value.
#' @export
count_covered_bases = function(gr){
  
  return(sum(width(reduce(gr, ignore.strand = T))))

}

#' Calculate the number of bases in the intersection of two GRanges objects
#'
#' @param gr1 A GRanges object.
#' @param gr2 A GRanges object.
#' @param ignore.strand TRUE or FALSE indicating whether strand should be ignored when calculating intersections. Default is TRUE.
#' @param overlap_measure One of "absolute", "proportion" or "jaccard" indicating whether to calculate 
#' the absolute size of the intersection in base pairs, the proportion base pairs of gr1 overlapping gr2 
#' or the Jaccard index of the intersection in terms of base pairs. Default value is "absolute".
#' @return A numeric value
#' @export
calculate_regions_intersections <- function(gr1, gr2, ignore.strand = TRUE, overlap_measure = "absolute"){
  
  # Check that inputs have the correct data type
  stopifnot(is(gr1, "GRanges"), is(gr2, "GRanges"), 
    S4Vectors::isTRUEorFALSE(ignore.strand), is(overlap_measure, "character"))
  
  # Check allowed value provided for overlap_measure
  match.arg(overlap_measure, c("absolute", "proportion", "jaccard"))
  
  # Create GRanges with the intersection and union of gr1 and gr2
  intersection <- GenomicRanges::intersect(gr1, gr2, ignore.strand = ignore.strand)
  union <- c(gr1, gr2)
  
  # Calculate proportion, Jaccard index or absolute overlap depending on overlap_measure
  if(overlap_measure == "proportion"){
    return(count_covered_bases(intersection)/count_covered_bases(gr1))
  } else if(overlap_measure == "jaccard"){
    return(count_covered_bases(intersection)/count_covered_bases(union))
  } else {
    return(count_covered_bases(intersection))
  }

}

### Expanding GRanges

#' Adjust GRanges
#'
#' Adjust ranges in a GRanges object upstream and downstream by specified numbers of bases. Takes strand into account.
#'
#' @param gr A GRanges object
#' @param upstream An integer specifying the the number of bases to adjust each region upstream of its start. Default is 0 bases
#' @param Downstream An An integer specifying the number of bases to adjust each region downstream of its end. Default is 0 bases
#' @return A GRanges object. 
#' @export
adjust_gr = function(gr, upstream = 0, downstream = 0) {
  
  # Check for each range if it's on the negative or positive strand
  strand_is_minus = as.logical(strand(gr) == "-")
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  
  # Adjust ranges based on whether they are on the positive or negative strand
  start(gr)[on_plus] = start(gr)[on_plus] - upstream
  start(gr)[on_minus] = start(gr)[on_minus] - downstream
  end(gr)[on_plus] = end(gr)[on_plus] + downstream
  end(gr)[on_minus] = end(gr)[on_minus] + upstream
  return(gr)
} 

### Filtering GRanges

#' Filter GRanges
#'
#' @param gr A GRanges object
#' @param column The name of a metadata column of gr to filter on.
#' @param values A vector of values of column to filter for.
#' @param invert Whether to filter for ranges without these values for the specified column. Default is FALSE.
#' @return A GRanges object. 
#' @export
filter_gr = function(gr, column, values, invert = FALSE){
  
  # Check that column is a metadata column of gr
  if(!column %in% names(mcols(gr))){stop("Provided metadata column name is not present in gr")}
  
  # Filter gr for ranges where column is in values, inverting the results if specified
  if(!invert){
    gr_filtered = gr[mcols(gr)[[column]] %in% values]
  } else {
    gr_filtered = gr[!mcols(gr)[[column]] %in% values]
  }
  
  return(gr_filtered)
  
}

#' Return the regions in a query GRanges which overlap a subject GRanges, keeping the metadata from query regions
#'
#' @param gr1 A GRanges object.
#' @param gr2 A GRanges object.
#' @return A GRanges object. 
#' @export
bedtools_intersect = function(query_gr, subject_gr){
  
  # Create a data.frame from query_gr and subtract 1 from start so that it is 0-based
  query_gr_df = data.frame(query_gr)
  query_gr_df$start = query_gr_df$start - 1
  
  # Create a temporary BED file and write query_gr_df there
  query_gr_tempbed = tempfile()
  data.table::fwrite(query_gr_df, query_gr_tempbed, sep = "\t", col.names = F, quote = F)
  
  # Reduce subject_gr to unique regions
  subject_gr = reduce(subject_gr, ignore.strand = T)
  
  # Create a data.frame from subject_gr and subtract 1 from start so that it is 0-based
  subject_gr_df = data.frame(subject_gr)
  subject_gr_df$start = subject_gr_df$start - 1
  
  # Create a temporary BED file and write subject_gr_df there
  subject_gr_tempbed = tempfile()
  data.table::fwrite(subject_gr_df, subject_gr_tempbed, sep = "\t", col.names = F, quote = F)
  
  # Run bedtools intersect and save to a temporary file
  intersect_bed = tempfile()
  system(paste("bedtools intersect -a", query_gr_tempbed, "-b", subject_gr_tempbed, ">", intersect_bed))
  
  # Read in the result of bedtools intersect and convert to a GRages and return
  intersect_df = data.table::fread(intersect_bed, col.names = names(query_gr_df))
  intersect_gr = makeGRangesFromDataFrame(intersect_df, starts.in.df.are.0based = T, keep.extra.columns = T)
  return(intersect_gr)
  
}
