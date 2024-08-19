#' Liftover a GRanges object
#'
#' Liftover a GRanges from one genome build to another using a provided liftover chain file
#'
#' @param genomic_region A GRanges object.
#' @param chain A "Chain" object to be used with rtracklayer::liftOver.
#' @param remove_non_mapping A logical value indicating whether to remove regions in the source genome.
#' which could not be mapped to the target genome. Default is TRUE. 
#' @param remove_one_to_many_mapping A logical value indicating whether to remove regions in the source genome.
#' which map to multiple regions in the target genome. Default is TRUE.
#' @param remove_many_to_one_mapping A logical value indicating whether to remove different regions in the source genome.
#' which map to the same region in the target genome i.e. remove regions which overlap after the liftover. Default is TRUE.
#' @param target_regions_overlap_filter An optional GRanges object used to filter the rowRanges by overlaps after liftover, 
#' for example CpG sites from the target genome. Any regions which do not overlap target_regions_overlap_filter will be removed.  
#' @param flatten_ranges A logical value indicating whether to convert the lifted over ranges from a 
#' GRangesList to GRanges if all remaining source regions can be uniquely mapped to the target genome. 
#' @return A GRanges object
#' @export
liftover_granges = function(genomic_regions, chain, remove_non_mapping = T, remove_one_to_many_mapping = T, 
  remove_many_to_one_mapping = T, target_regions_overlap_filter = NULL, flatten_ranges = T){
  
  # Liftover genomic_regions using specified liftover chain file
  liftover_ranges = rtracklayer::liftOver(genomic_regions, chain)
  
  # Put seqlevels of liftover_ranges in the same order as genomic_region
  seqlevels(liftover_ranges) = seqlevels(genomic_regions)
  
  # Initialize selected regions to all liftover_ranges
  selected_ranges = seq_along(liftover_ranges)
  
  # Get the number of regions in the target genome each region in the source genome matches to
  mappings_count = lengths(liftover_ranges)
  
  # Remove non-mapping regions from selected_ranges if specified
  if(remove_non_mapping){
    non_mapping_regions = which(mappings_count == 0)
    selected_ranges = setdiff(selected_ranges, non_mapping_regions)
  }
  
  # Remove one-to-many mapping regions from selected_ranges if specified
  if(remove_one_to_many_mapping){
    one_to_many_mapping_regions = which(mappings_count > 1)
    selected_ranges = setdiff(selected_ranges, one_to_many_mapping_regions)
  }
  
  # Remove many-to-one mapping regions from selected_ranges if specified
  if(remove_many_to_one_mapping){
    self_overlaps = countOverlaps(liftover_ranges, liftover_ranges)
    many_to_one_mapping_regions = which(self_overlaps > 1)
    selected_ranges = setdiff(selected_ranges, many_to_one_mapping_regions)
  }
  
  # Identify regions which overlap target_regions_overlap_filter is provided
  if(!is.null(target_regions_overlap_filter)){
    target_regions_overlaps = countOverlaps(liftover_ranges, target_regions_overlap_filter)
    non_target_overlaps = which(target_regions_overlaps < 1)
    selected_ranges = setdiff(selected_ranges, non_target_overlaps)
  }
  
  # Subset liftover_ranges for selected ranges
  liftover_ranges = liftover_ranges[selected_ranges]
  
  # Flatten liftover_ranges if specified and there are no one-to-many mappings
  if(flatten_ranges & length(unlist(liftover_ranges)) == length(liftover_ranges)){
    liftover_ranges = unlist(liftover_ranges)
  }
  
  # Return liftover_ranges 
  return(liftover_ranges)
  
}