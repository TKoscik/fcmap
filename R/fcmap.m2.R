fcmap.m2 <- function(epi.nii,
                     cortex,
                     mask_ROI,
                     save_file_name,
                     cluster = NULL,
                     cluster_limit = NULL) {

  epi_dims <- nii.dims(epi.nii)
  Cortex <- read.nii.volume(cortex, 1)
  ROI <- read.nii.volume(mask_ROI, 1)
  domain_number <- max(Cortex)

  domain <- as.list(1:domain_number)
  names(domain) <- paste("Domain", 1:length(domain), sep = ".")

  dictionary <- data.frame(matrix(NA, length(Cortex), 3))
  colnames(dictionary) <- c("Row","Column","Level")

  dictionary$Row <- rep(1:epi_dims[1], times = epi_dims[2])
  dictionary$Column <- rep(1:epi_dims[2], each = epi_dims[1])
  dictionary$Level <- rep(1:epi_dims[3], each = (epi_dims[1]*epi_dims[2]))

  Cortex_spots <- which(Cortex > 0)
  ROI_spots <- which(ROI > 0)

  output_domains <- data.frame(matrix(NA, epi_dims[4], domain_number))
  output_ROI <- data.frame(matrix(NA, epi_dims[4], length(ROI_spots)))

  for (i in 1:epi_dims[4]){
    epi_i <- read.nii.volume(epi.nii, vol.num = i)
    output_ROI[i,] <- epi_i[ROI_spots]

    for (j in 1:domain_number){
      domain[[j]] <- which(Cortex == j)
      output_domains[i,j] <- mean(epi_i[domain[[j]]], na.rm = T)
    }
  }

  #look at correlations
  output_cor <- data.frame(matrix(NA, length(ROI_spots), domain_number+1))

  for (i in 1:length(ROI_spots)){
    voxel_i <- output_ROI[,i]

    if (sum(voxel_i, na.rm = T) == 0){
      output_cor[i,length(output_cor[1,])] <- 0
    } else {
      for (j in 1:domain_number){
        output_cor[i,j] <- unname(unlist(rcorr(voxel_i, output_domains[,j], type = "pearson"))[2])
      }

      z <- unlist(output_cor[i,1:domain_number])
      output_cor[i,length(output_cor[1,])] <- which(max(z, na.rm = T) == z)
    }
  }
  #create two files to write into, one with cluster and one w/o. See if works.
  #Clustering Part (optional)
  if (is.null(cluster) == F){
    if (is.null(cluster_limit) == T){
      cluster_limit <- 12
    }

    cluster_output <- data.frame(matrix(NA, length(ROI_spots), 4))
    colnames(cluster_output) <- c("Position","Initial","Cluster","Final")
    cluster_output$Position <- ROI_spots
    cluster_output$Initial <- output_cor[,length(output_cor[1,])]

    unique_mapped_domains <- unique(output_cor[,length(output_cor[1,])])
    number_mapped_domains <- length(unique_mapped_domains)

    domain.voxel.num <- data.frame(matrix(NA,number_mapped_domains,2))
    domain.voxel.num[,1] <- unique_mapped_domains
    colnames(domain.voxel.num) <- c("Domain","Num_Voxels")

    for (i in 1:number_mapped_domains){
      domain.voxel.num[i,2] <- length(which(domain.voxel.num[i,1] == output_cor[,length(output_cor[1,])]))
    }

    clustering_domains <- domain.voxel.num %>% filter(Num_Voxels > cluster_limit) #spot with 10, 12, 12 didn't work

    for (i in 1:length(clustering_domains$Domain)){
      domain.i <- clustering_domains$Domain[i]
      voxels.i <- which(output_cor[,length(output_cor[1,])] == domain.i)
      my_data.i <- data.frame(matrix(NA,length(voxels.i),4))
      colnames(my_data.i) <- c("X","Y","Z","New_Cluster")

      for (j in 1:length(voxels.i)){
        spot.j <- cluster_output$Position[voxels.i[j]]
        my_data.i$X[j] <- dictionary$Row[spot.j]
        my_data.i$Y[j] <- dictionary$Column[spot.j]
        my_data.i$Z[j] <- dictionary$Level[spot.j]
      }

      cluster_data.i <- NbClust(data = my_data.i[,1:3], distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index ="all")
      my_data.i$New_Cluster <- cluster_data.i$Best.partition
      if (i == 1){
        cluster_output$Cluster[voxels.i] <- my_data.i$New_Cluster
      } else {
        counter_max <- max(cluster_output$Cluster, na.rm = T)
        cluster_output$Cluster[voxels.i] <- (my_data.i$New_Cluster + counter_max)
      }
    }
    ###Now I redo WTA, but instead of voxel it is cluster
    domains_clustered <- max(cluster_output$Cluster, na.rm = T)
    cluster_TS_means <- data.frame(matrix(NA, epi_dims[4], domains_clustered))
    for (i in 1:domains_clustered){
      cluster.i <- which(cluster_output$Cluster == i)
      if (length(cluster.i) == 1){
        cluster_TS_means[,i] <- output_ROI[,cluster.i]
      } else {
        cluster_series.i <- output_ROI[,cluster.i]
        cluster_TS_means[,i] <- rowMeans(cluster_series.i)
      }
    }
    ###Look at correlations for the clusters to the cortex
    output_cluster_cor <- data.frame(matrix(NA, domains_clustered, domain_number+1))

    for (i in 1:domains_clustered){
      cluster_i <- cluster_TS_means[,i]

      for (j in 1:domain_number){
        output_cluster_cor[i,j] <- unname(unlist(rcorr(cluster_i, output_domains[,j], type = "pearson"))[2])
      }

      z <- unlist(output_cluster_cor[i,1:domain_number])
      output_cluster_cor[i,(domain_number+1)] <- which(max(z, na.rm = T) == z)

      cluster.i.new <- which(cluster_output$Cluster == i)
      cluster_output$Final[cluster.i.new] <- which(max(z, na.rm = T) == z)
    }
    #Adding in unclustered
    position_interest_unclustered <- which(is.na(cluster_output$Final) == T)
    voxels_unclustered <- cluster_output %>% filter(is.na(Final) == T)
    voxels_unclustered$P.interest <- position_interest_unclustered
    unique_matched_domains <- na.omit(unique(cluster_output$Final))
    num_matched_domains <- length(unique_matched_domains)

    unclustered_cor <- data.frame(matrix(NA, length(position_interest_unclustered), num_matched_domains + 1))

    for (i in 1:length(position_interest_unclustered)){
      unclustered_voxel.i <- unlist(output_ROI[voxels_unclustered$P.interest[i]])

      for (j in 1:num_matched_domains){
        voxels_domain.j <- which(cluster_output$Final == unique_matched_domains[j])
        voxels_domain.j_TS <- output_ROI[,voxels_domain.j]
        voxels_domain.j_TS_mean <- rowMeans(voxels_domain.j_TS, na.rm = T)

        unclustered_cor[i,j] <- unname(unlist(rcorr(unclustered_voxel.i, voxels_domain.j_TS_mean, type = "pearson"))[2])
      }

      z <- unlist(unclustered_cor[i,1:num_matched_domains])
      unclustered_cor[i,(num_matched_domains+1)] <- which(max(z, na.rm = T) == z)
    }
    cluster_output$Final[position_interest_unclustered] <- unclustered_cor[,length(unclustered_cor[1,])]
    output_cor[,length(output_cor[1,])] <- cluster_output$Final
  }

  file.copy(mask_ROI, save_file_name)

  for (i in 1:length(ROI_spots)){
    voxel.i <- ROI_spots[i]
    cor.x <- dictionary$Row[voxel.i]
    cor.y <- dictionary$Column[voxel.i]
    cor.z <- dictionary$Level[voxel.i]

    write.nii.voxel(save_file_name, c(cor.x, cor.y, cor.z, 1), output_cor[i,length(output_cor[1,])])
  }
}
