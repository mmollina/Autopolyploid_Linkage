#####
## Map Construction using polymapR
#####
build_map_polymapR<-function(dat, 
                             n.chr = 1, 
                             LOD_bridge = 3, 
                             LOD_assign_LG = 3, 
                             LOD_assign_HM = 3)
{
  segregating_data <- convert_marker_dosages(dosage_matrix = dat, ploidy = 6)
  #####
  ## Parent 1
  #####
  SN_SN_P1 <- linkage(dosage_matrix = segregating_data,
                      markertype1 = c(1,0),
                      target_parent = "P1",
                      other_parent = "P2",
                      ploidy = 6,
                      ncores = 16,
                      pairing = "random")
  
  SN_SN_P1_coupl <- SN_SN_P1[SN_SN_P1$phase == "coupling",] # select only markerpairs in coupling
  
  P1_homologues_1 <- cluster_SN_markers(linkage_df = SN_SN_P1, 
                                        LOD_sequence = c(1:20), 
                                        LG_number = n.chr,
                                        ploidy = 6,
                                        parentname = "P1",
                                        plot_network = FALSE,
                                        plot_clust_size = FALSE)
  
  pz<-names(which(sapply(P1_homologues_1, function(x) length(table(x$cluster))) == 6*n.chr)[1])
  length(table(P1_homologues_1[[pz]]$cluster))
  
  SN_DN_P1 <- linkage(dosage_matrix = segregating_data, 
                      markertype1 = c(1,0),
                      markertype2 = c(2,0),
                      target_parent = "P1",
                      other_parent = "P2",
                      ploidy = 6,
                      pairing = "random")
  
  LGHomDf_P1_1 <- bridgeHomologues(cluster_stack = P1_homologues_1[[pz]], 
                                   linkage_df = SN_DN_P1, 
                                   LOD_threshold = LOD_bridge, 
                                   automatic_clustering = TRUE, 
                                   LG_number = n.chr,
                                   parentname = "P1")
  
  
  SN_SS_P1 <- linkage(dosage_matrix = segregating_data, 
                      markertype1 = c(1,0),
                      markertype2 = c(1,1),
                      target_parent = "P1",
                      other_parent = "P2",
                      ploidy = 6,
                      pairing = "random")
  
  P1_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P1,
                                          LG_hom_stack = LGHomDf_P1_1,
                                          SN_colname = "marker_a",
                                          unassigned_marker_name = "marker_b",
                                          phase_considered = "coupling",
                                          LG_number = n.chr,
                                          LOD_threshold = LOD_assign_LG,
                                          ploidy = 6)
  
  head(P1_SxS_Assigned)
  P1_DxN_Assigned <- assign_linkage_group(linkage_df = SN_DN_P1,
                                          LG_hom_stack = LGHomDf_P1_1,
                                          SN_colname = "marker_a",
                                          unassigned_marker_name = "marker_b",
                                          phase_considered = "coupling",
                                          LG_number = n.chr,
                                          LOD_threshold = LOD_assign_LG,
                                          ploidy = 6)
  
  marker_assignments_P1 <- homologue_lg_assignment(dosage_matrix = segregating_data,
                                                   assigned_list = list(P1_SxS_Assigned, 
                                                                        P1_DxN_Assigned),
                                                   assigned_markertypes = list(c(1,1), c(2,0)),
                                                   LG_hom_stack = LGHomDf_P1_1,
                                                   target_parent = "P1",
                                                   other_parent = "P2",
                                                   ploidy = 6,
                                                   pairing = "random",
                                                   convert_palindrome_markers = FALSE,
                                                   LG_number = n.chr,
                                                   LOD_threshold = LOD_assign_HM,
                                                   write_intermediate_files = FALSE)
  
  #####
  ## Parent 2
  #####
  SN_SN_P2 <- linkage(dosage_matrix = segregating_data,
                      markertype1 = c(1,0),
                      target_parent = "P2",
                      other_parent = "P1",
                      ploidy = 6,
                      ncores = 16,
                      pairing = "random")
  
  SN_SN_P2_coupl <- SN_SN_P2[SN_SN_P2$phase == "coupling",] # select only markerpairs in coupling
  
  P2_homologues_1 <- cluster_SN_markers(linkage_df = SN_SN_P2_coupl, 
                                        LOD_sequence = c(1:20), 
                                        LG_number = n.chr,
                                        ploidy = 6,
                                        parentname = "P2",
                                        plot_network = FALSE,
                                        plot_clust_size = FALSE)
  qz<-names(which(sapply(P2_homologues_1, function(x) length(table(x$cluster))) == 6*n.chr)[1])
  length(table(P2_homologues_1[[qz]]$cluster))
  
  SN_DN_P2 <- linkage(dosage_matrix = segregating_data, 
                      markertype1 = c(1,0),
                      markertype2 = c(2,0),
                      target_parent = "P2",
                      other_parent = "P1",
                      ploidy = 6,
                      pairing = "random")
  
  
  LGHomDf_P2_1 <- bridgeHomologues(cluster_stack = P2_homologues_1[[qz]], 
                                   linkage_df = SN_DN_P2, 
                                   LOD_threshold = LOD_bridge,
                                   automatic_clustering = TRUE, 
                                   LG_number = n.chr,
                                   parentname = "P2")
  
  table(LGHomDf_P2_1$LG, LGHomDf_P2_1$homologue)
  
  SN_SS_P2 <- linkage(dosage_matrix = segregating_data, 
                      markertype1 = c(1,0),
                      markertype2 = c(1,1),
                      target_parent = "P2",
                      other_parent = "P1",
                      ploidy = 6,
                      pairing = "random")
  
  P2_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P2,
                                          LG_hom_stack = LGHomDf_P2_1,
                                          SN_colname = "marker_a",
                                          unassigned_marker_name = "marker_b",
                                          phase_considered = "coupling",
                                          LG_number = n.chr,
                                          LOD_threshold = LOD_assign_LG,
                                          ploidy = 6)
  
  head(P2_SxS_Assigned)
  P2_DxN_Assigned <- assign_linkage_group(linkage_df = SN_DN_P2,
                                          LG_hom_stack = LGHomDf_P2_1,
                                          SN_colname = "marker_a",
                                          unassigned_marker_name = "marker_b",
                                          phase_considered = "coupling",
                                          LG_number = n.chr,
                                          LOD_threshold = LOD_assign_LG,
                                          ploidy = 6)
  
  marker_assignments_P2 <- homologue_lg_assignment(dosage_matrix = segregating_data,
                                                   assigned_list = list(P2_SxS_Assigned, 
                                                                        P2_DxN_Assigned),
                                                   assigned_markertypes = list(c(1,1), c(2,0)),
                                                   LG_hom_stack = LGHomDf_P2_1,
                                                   target_parent = "P2",
                                                   other_parent = "P1",
                                                   ploidy = 6,
                                                   pairing = "random",
                                                   convert_palindrome_markers = FALSE,
                                                   LG_number = n.chr,
                                                   LOD_threshold = LOD_assign_HM,
                                                   write_intermediate_files = FALSE)
  
  #####
  ## Merging Parents 1 and 2
  #####
  
  marker_assignments <- check_marker_assignment(marker_assignments_P1,marker_assignments_P2)
  
  all_linkages_list_P1 <- finish_linkage_analysis(marker_assignment = marker_assignments$P1,
                                                  dosage_matrix = segregating_data,
                                                  target_parent = "P1",
                                                  other_parent = "P2",
                                                  convert_palindrome_markers = FALSE,
                                                  ploidy = 6,
                                                  pairing = "random",
                                                  LG_number = n.chr) 
  
  all_linkages_list_P2 <- finish_linkage_analysis(marker_assignment = marker_assignments$P2,
                                                  dosage_matrix = segregating_data,
                                                  target_parent = "P2",
                                                  other_parent = "P1",
                                                  convert_palindrome_markers = TRUE, # convert 3.1 markers
                                                  ploidy = 6,
                                                  pairing = "random",
                                                  LG_number = n.chr)
  
  linkages <- list()
  for(lg in names(all_linkages_list_P1)){
    linkages[[lg]] <- rbind(all_linkages_list_P1[[lg]], all_linkages_list_P2[[lg]])
  }
  
  #####
  ## MDS algorithm
  #####
  integrated.maplist <- MDSMap_from_list(linkages)
  #####
  ## Creating a phased map
  #####
  phased.maplist <- create_phased_maplist(maplist = integrated.maplist,
                                          dosage_matrix.conv = segregating_data,
                                          N_linkages = n.chr,
                                          ploidy = 6,
                                          marker_assignment.1 = marker_assignments$P1,
                                          marker_assignment.2 = marker_assignments$P2)
  return(list(integrated = integrated.maplist, phased = phased.maplist))
}

