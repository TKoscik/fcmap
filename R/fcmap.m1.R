fcmap.m1 <- function(epi.nii,
                     network_cortex,
                     network_ROI,
                     save_file_name) {

  epi_dims <- nii.dims(epi.nii)
  network_cortex <- read.nii.volume(network_cortex, 1)
  network_ROI <- read.nii.volume(network_ROI, 1)

  domains_cortex <- max(network_cortex)
  domainC <- as.list(1:domains_cortex)
  names(domainC) <- paste("DomainC", 1:length(domains_cortex), sep = ".")

  domains_ROI <- max(network_ROI) #I want to WTA these domains to the cortex domains
  domainROI <- as.list(1:domains_ROI)
  names(domainROI) <- paste("DomainROI", 1:length(domains_ROI), sep = ".")

  output_cortex <- data.frame(matrix(NA, epi_dims[4], domains_cortex))
  output_ROI <- data.frame(matrix(NA, epi_dims[4], domains_ROI))

  for (i in 1:epi_dims[4]){
    epi_i <- read.nii.volume(epi.nii, vol.num = i)

    for (j in 1:domains_cortex){
      domain.j <- which(network_cortex == j)
      output_cortex[i,j] <- mean(epi_i[domain.j], na.rm = T)
    }

    for (k in 1:domains_ROI){
      domain.k <- which(network_ROI == k)
      output_ROI[i,k] <- mean(epi_i[domain.k], na.rm = T)
    }
  }

  output_corr <- data.frame(matrix(NA, domains_ROI, domains_cortex+1))

  for (i in 1:domains_ROI){
    ROI.i <- output_ROI[,i]

    for (j in 1:domains_cortex){
      cortex.j <- output_cortex[,j]
      output_corr[i,j] <- unname(unlist(rcorr(ROI.i, cortex.j, type = "pearson"))[2])
    }

    z <- unlist(output_corr[i,1:domains_cortex])
    output_corr[i,length(output_corr[1,])] <- which(max(z, na.rm = T) == z)
  }

  dictionary <- data.frame(matrix(NA, length(network_cortex), 3))
  colnames(dictionary) <- c("Row","Column","Level")

  dictionary$Row <- rep(1:epi_dims[1], times = epi_dims[2])
  dictionary$Column <- rep(1:epi_dims[2], each = epi_dims[1])
  dictionary$Level <- rep(1:epi_dims[3], each = (epi_dims[1]*epi_dims[2]))

  file.copy(network_ROI, save_file_name)

  for (i in 1:domains_ROI){
    domain.i <- which(network_ROI == i)
    WTA.i <- output_corr[i,length(output_corr[1,])]

    for (j in 1:length(domain.i)){
      voxel.j <- domain.i[j]
      cor.x <- dictionary$Row[voxel.j]
      cor.y <- dictionary$Column[voxel.j]
      cor.z <- dictionary$Level[voxel.j]

      write.nii.voxel(save_file_name, c(cor.x, cor.y, cor.z, 1), WTA.i)
    }
  }
}
