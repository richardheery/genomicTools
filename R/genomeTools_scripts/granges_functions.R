#' Calculate distances of query GRanges upstream or downstream of subject GRanges
#'
#' @param query_gr A GRanges object
#' @param subject_gr A GRanges object. Upstream and downstream refer to the strand of subject_gr
#' @return A numeric vector on distances
#' @export
stranded_distance = function(query_gr, subject_gr){
  
  # Check that query_gr and subject_gr are of the correct length
  if(length(query_gr) != length(subject_gr)){
    if(length(query_gr) > 1 && length(subject_gr) > 1){
      stop("query_gr and subject_gr must either have the same length or one of them must have length 1")
    }
  }
  
  # Initialize a vector of zeros with length equal to the greater length of query_gr or subject_gr
  d = rep(0, max(length(query_gr), length(subject_gr)))
  
  # Find the length of the gap between two GRanges, setting to 0 if they overlap
  d[end(query_gr) < start(subject_gr)] = (end(query_gr) - start(subject_gr))[end(query_gr) < start(subject_gr)]
  d[start(query_gr) > end(subject_gr)] = (start(query_gr) - end(subject_gr))[start(query_gr) > end(subject_gr)]
  
  # If ranges are on different strands, set distance as NA
  d[as.vector(seqnames(query_gr)) != as.vector(seqnames(subject_gr))] = NA
  
  # Set the strand of subject_gr as 1 if it is on the "+" strand and -1 if it is on the "-" strand
  subject_strand = as.numeric(paste0(strand(subject_gr), 1)) 
  
  # Get the signed distance of query_gr from subject_gr and return
  d = d * subject_strand
  return(d)

}

#' Find locations of ranges relative to a reference positions
#'
#' @param gr A GRanges object
#' @param reference_positions A GRanges object. Each range should have width 1. Upstream and downstream refer to reference_positions
#' @param gr_category A category that ranges belong to
#' @return A GRanges object
#' @export
relative_ranges = function(gr, reference_positions, gr_category = NULL){
  
  # Check that all reference positions have width 1
  if(!all(width(reference_positions) == 1)){stop("All reference_positions should have width 1")}

  # Get distances from start and end of ranges in gr from reference_positions
  relative_start = genomeTools::stranded_distance(resize(gr, 1, fix = "start"), reference_positions)
  relative_end = genomeTools::stranded_distance(resize(gr, 1, fix = "end"), reference_positions)
  
  # Create an IRanges with the relative distances
  relative_iranges = IRanges::IRanges(pmin(relative_start, relative_end), pmax(relative_start, relative_end), category = gr_category)
  
  # Convert IRanges to GRanges with "relative" as seqnames
  relative_granges = GenomicRanges::GRanges(seqnames = "relative", ranges = relative_iranges)
  
  # Add metadata from gr to relative_granges
  mcols(relative_granges) = mcols(gr)
  
  # Return relative_granges
  return(relative_granges)
  
}

# Recreate absoulte ranges from relative ranges
#'
#' @param relative_ranges An IRanges object
#' @param reference_positions A GRanges object. Upstream and downstream refer to the strand of subject_gr
#' @return A GRanges object
#' @export
absolute_ranges = function(relative_ranges, reference_positions){

  # Get strand direction and start for reference_positions 
  reference_seqnames = seqnames(reference_positions)
  reference_strand = strand(reference_positions)
  reference_starts = start(reference_positions)
  
  # Get absolute starts and ends
  absolute_starts = start(relative_ranges) + reference_starts
  absolute_ends = end(relative_ranges) + reference_starts
  
  # Create a GRanges 
  gr = GRanges(seqnames = reference_seqnames, strand = reference_strand, 
    ranges = IRanges(start = absolute_starts, end = absolute_ends)
)
  
  return(gr)
}

#' Create relative bins 
#'
#' @param bin_start Start of bins
#' @param bin_end End of bins
#' @param bin_step Size of bins
#' @return A GRanges object
#' @export
create_relative_bins = function(bin_start, bin_end, bin_step){
  
  bin_iranges = IRanges(
    start = seq(bin_start, bin_end - bin_step, bin_step), 
    end = seq(bin_start + bin_step, bin_end, bin_step) - 1
    )
  
  bin_granges = GRanges(seqnames = "relative", ranges = bin_iranges)
  bin_granges$bin_center = start(bin_granges) + width(bin_granges)/2
  
  
  return(bin_granges)

}

#' Create relative bins and count the overlaps of relative ranges with relative bins
#'
#' @param relative_ranges A GRanges object returned by relative_ranges
#' @param bin_start Start of bins
#' @param bin_end End of bins
#' @param bin_step Size of bins
#' @param category A category associated with relative_ranges
#' @return A data.frame with the the number of ranges of each category in each bin
#' @export
bin_relative_ranges = function(relative_ranges, bin_start, bin_end, bin_step, category){
  
  bin_ranges = IRanges(
    start = seq(bin_start, bin_end - bin_step, bin_step), 
    end = seq(bin_start + bin_step, bin_end, bin_step) - 1
    )
  
  if(!is.null(category)){
    ranges_list = 
    split(relative_ranges, category)
  } else {ranges_list = list(relative_ranges)}
  
  overlap_df = data.frame(lapply(ranges_list, function(x)
    setNames(countOverlaps(bin_ranges, ranges(x)), start(bin_ranges) + width(bin_ranges)/2)))
  
  overlap_df = tibble::rownames_to_column(overlap_df, "bin_center")
  overlap_df$bin_center = as.numeric(overlap_df$bin_center)
  
  return(overlap_df)

}

#' Calculate the number of bases in the intersection of all pairs of GRanges from a list
#'
#' @param grl A list of GRanges objects
#' @param ignore.strand A logical value indicating whether strand should be ignored when calculating intersections. Default is TRUE.
#' @param as.proportion A logical value indicating whether the overlap should be returned as a proportion of the GRanges indicated by the columns. Default value is FALSE. 
#' @return A matrix with the size of the overlaps between all pairs of GRanges in grl
#' @export
calculate_regions_intersection_list = function(grl, ignore.strand = T, as.proportion = F){
  
  overlap_matrix = matrix(NA, nrow = length(grl), ncol = length(grl), dimnames = list(names(grl), names(grl)))
  
  for(i in names(grl)){
    for(j in names(grl)){
      overlap_matrix[j, i] = genomeTools::calculate_regions_intersections(grl[[i]], grl[[j]], ignore.strand = ignore.strand, as.proportion = as.proportion)
    }
  }
  
  return(overlap_matrix)

}

#' Calculate the Jaccard index of the intersections of all pairs of GRanges from a list
#'
#' @param grl A list of GRanges objects
#' @param ignore.strand A logical value indicating whether strand should be ignored when calculating intersections. Default is TRUE.
#' @param as.proportion A logical value indicating whether the overlap should be returned as a proprtion of gr1. Default value is FALSE. 
#' @return A matrix with the size of the overlaps between all pairs of GRanges in grl
#' @export
calculate_regions_jaccard_list = function(grl, ignore.strand = T){
  
  jaccard_matrix = matrix(NA, nrow = length(grl), ncol = length(grl), dimnames = list(names(grl), names(grl)))
  
  for(i in names(grl)){
    for(j in names(grl)){
      jaccard_matrix[i, j] = genomeTools::calculate_regions_jaccard(grl[[i]], grl[[j]], ignore.strand = ignore.strand)
    }
  }
  
  return(jaccard_matrix)

}

#' #' Calculate the number of ranges in one GRanges object overlapping another
#' #'
#' #' @param gr1 A GRanges object
#' #' @param gr1 A GRanges object
#' #' @param ignore.strand A logical value indicating whether strand should be ignored when calculating overlap. Default is TRUE.
#' #' @param as.proportion A logical value indicating whether the overlap should be returned as a proportion of gr1. Default value is FALSE. 
#' #' @return An numeric value
#' #' @export
#' calculate_regions_overlaps = function(gr1, gr2, ignore.strand = T, as.proportion = F){
#'   
#'   overlap = length(IRanges::subsetByOverlaps(gr1, gr2, ignore.strand = ignore.strand))
#'   
#'   if(as.proportion){
#'     return(overlap/length(gr1))
#'   } else {
#'     return(overlap)
#'   }
#' 
#' }

#' #' Calculate the number of ranges of each GRanges object in a list overlapping the others
#' #'
#' #' @param grl A list of GRanges objects
#' #' @param ignore.strand A logical value indicating whether strand should be ignored when calculating overlaps. Default is TRUE.
#' #' @param as.proportion A logical value indicating whether the overlap should be returned as a proprtion of gr1. Default value is FALSE. 
#' #' @return A matrix with the size of the overlaps between all pairs of GRanges in grl
#' #' @export
#' calculate_regions_overlap_list = function(grl, ignore.strand = T, as.proportion = F){
#'   
#'   overlap_matrix = matrix(NA, nrow = length(grl), ncol = length(grl), dimnames = list(names(grl), names(grl)))
#'   
#'   for(i in names(grl)){
#'     for(j in names(grl)){
#'       overlap_matrix[i, j] = genomeTools::calculate_regions_overlaps(grl[[i]], grl[[j]], ignore.strand = ignore.strand, as.proportion = as.proportion)
#'     }
#'   }
#'   
#'   return(overlap_matrix)
#' 
#' }

#' Calculate the proportion of each region in a query set overlapping a subject set
#'
#' @param query_gr A GRanges object. The proportion of each region overlapping regions from subject_gr will be calcualted. 
#' @param subject_gr A GRanges object
#' @param ignore.strand A logical value indicating whether to ignore strand when calculating overlaps. Default is TRUE. 
#' @return A matrix with the size of the overlaps between all pairs of GRanges in grl
#' @export
individual_proportion_overlaps = function(query_gr, subject_gr, ignore.strand = T){
  
  # Set names for query_gr if they are missing
  if(is.null(names(query_gr))){
    names(query_gr) = paste0("region_", seq_along(query_gr))
  }
  
  # Merge overlapping regions of subject_gr
  subject_gr = reduce(subject_gr, ignore.strand = ignore.strand)
  
  # Find proportion of regions in query_gr involved in each overlap with subject_gr
  overlaps_df = with(data.frame(findOverlaps(query_gr, subject_gr, ignore.strand = ignore.strand)), 
    data.frame(
      name = names(query_gr[queryHits]), 
      overlap = width(pintersect(query_gr[queryHits], subject_gr[subjectHits]))/width(query_gr[queryHits]))
    )
  
  # Get the sum of the proportions for each region from query_gr
  overlaps_df = data.frame(dplyr::summarise(dplyr::group_by(overlaps_df, name), 
    overlap = sum(overlap)))
  
  # Add back in non-overlapping regions
  overlaps_df$name = factor(overlaps_df$name, levels = names(query_gr))
  overlaps_df = tidyr::complete(overlaps_df, name, fill = list(overlap = 0))
  
  # Put in the same order as query_gr and return
  overlaps_df = overlaps_df[match(names(query_gr), overlaps_df$name), ]
  
}

#' Calculate the number of ranges of each GRanges object in a list overlapping the others with a minimum overlap proportion threshold
#'
#' @param grl A list of GRanges objects
#' @param ignore.strand A logical value indicating whether strand should be ignored when calculating overlaps. Default is TRUE.
#' @param overlap_threshold Minimum overlap proprtion of the query ranges with the subject to consider them overlapping. Default is 0.5. 
calculate_regions_overlap_list = function(grl, ignore.strand = T, overlap_threshold = 0.5){
  
  overlap_matrix = matrix(NA, nrow = length(grl), ncol = length(grl), dimnames = list(names(grl), names(grl)))
  
  for(i in names(grl)){
    for(j in names(grl)){
      overlap_matrix[i, j] = 
        sum(genomeTools::individual_proportion_overlaps(grl[[i]], grl[[j]], ignore.strand = ignore.strand)$overlap > overlap_threshold)/length(grl[[i]])
    }
  }
  
  return(overlap_matrix)
  
}

#' Perform a Chi-squared test for the overlap of CpG sites in a set of query regions with a set of subject regions comapred to a set of control regions
#' 
#' @param query_regions A GRanges object
#' @param subject_regions A GRanges object 
#' @param control_regions A GRanges object 
#' @param background_cpgs A GRanges object with the location of all genomic CpG sites
#' @param alternative The alternative hypothesis. Must be one of "two.sided", "greater" or "less". "greater" is default value.
#' @return The results of the chi-squared test
#' @export
query_subject_cpg_overlap_test = function(query_regions, subject_regions, control_regions, background_cpgs, alternative = "greater"){
  
  # Check that alternative is one of the permitted values
  match.arg(arg = alternative, choices = c("two.sided", "greater", "less"), several.ok = F)
  
  # Find CpG sites overlapping the query, subject and control regions
  query_regions_cpgs = subsetByOverlaps(background_cpgs, query_regions, ignore.strand = T)$name
  subject_regions_cpgs = subsetByOverlaps(background_cpgs, subject_regions, ignore.strand = T)$name
  control_regions_cpgs = subsetByOverlaps(background_cpgs, control_regions, ignore.strand = T)$name
  
  # Find the number of CpGs in common between query_regions and subject_regions and control_regions and subject_regions
  query_subject_overlap = intersect(query_regions_cpgs, subject_regions_cpgs)
  control_subject_overlap = intersect(control_regions_cpgs, subject_regions_cpgs)
  
  # Perform a chi-square test
  chi_test_result = prop.test(
    x = c(length(query_subject_overlap), length(control_subject_overlap)), 
    n = c(length(query_regions_cpgs), length(control_regions_cpgs)), 
    correct = F, alternative = alternative
  )
  
  results_df = data.frame(broom::tidy(chi_test_result))
  results_df$query_subject_overlap = length(query_subject_overlap)
  results_df$control_subject_overlap = length(control_subject_overlap)
  
  return(results_df)
  
}

hg38_ucsc_seqinfo = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[
  GenomeInfoDb::extractSeqlevelsByGroup(species = "Homo_sapiens", style = "UCSC", group = "all"), ]

#' A list of colour vectors
#'
#' @format A Seqinfo object for hg38 for UCSC
"hg38_ucsc_seqinfo"

hg19_ucsc_seqinfo = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[
  GenomeInfoDb::extractSeqlevelsByGroup(species = "Homo_sapiens", style = "UCSC", group = "all"), ]

#' A list of colour vectors
#'
#' @format A Seqinfo object for hg19 for UCSC
"hg19_ucsc_seqinfo"
