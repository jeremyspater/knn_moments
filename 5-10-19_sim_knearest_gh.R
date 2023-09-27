#Github version of methods-paper code (27 sep '23); extraneous comments removed
#Function to calculate K-nearest

library(geosphere)

#Function to calculate K-nearest
k_sim = function(dat, k){
  dat_out = dat #output data frame
  
  #all_mat = distm(dat)
  #just use euclidean; we are using meters
  #all_mat = dist(dat) %>% as.matrix()
  all_mat = dist(dat[,c('long','lat')]) %>% as.matrix() #changed 5-15-19 because old version was causing problems for clustered neighborhoods
  diag(all_mat) = NA
  #starting 5-21-18: if k > kmax, set to NA
  kmax = max(rowSums(!is.na(all_mat))) #number of non-NA entries in distance matrix; equal to # of non-na-geocodes -1
  if(k > kmax){ #if k is higher than total number of neighbors, set to na
    dat_out[,paste0('Nearest_',k,'_Same')] = NA # Set to na
  }else{
    dat_out[,paste0('Nearest_',k,'_Same')] = NA # Initialize
    for(i in 1:nrow(dat_out)){ #loop over ppl in neighborhood
      
      #add logic 5-21-18 to make these NA if the geocode is NA
      if(is.na(dat_out$lat[i])){
        
        dat_out[,paste0('Nearest_',k,'_Same')] = NA

      }else{ #if we have a geocode
        Prox1R = all_mat[i,] #ith person's distances to others in neighborhood
        dat_out[i,paste0('Nearest_',k,'_Same')] = sum(dat_out[order(Prox1R)[1:k],'minor'] == dat_out$minor[i], na.rm = T) #count number of nearest-K who are same minority status
      }
    }
  }
  return(dat_out)
}

