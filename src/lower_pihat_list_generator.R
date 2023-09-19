# Created to aid in the automation of the GDC pipeline process
# Idea is to go into the file that is output by step 6 regarding 
# look at the file pihat_min0.2_in_founders.genome and choose the 
# one/ones with the lower z0/z1
# Then record their FID and IID and create a .txt that will be used to remove them 
# from the genomics data 


pihat_dat <- read.csv("pihat_min0.2_in_founders.genome", sep="")

# Comparing Z0 and Z1, noting which is lower...
FID_to_rm = c()
IID_to_rm = c()

for (j in 1:nrow(pihat_dat)) {
  FID_to_rm[j] = ifelse(pihat_dat$Z0[j] < pihat_dat$Z1[j], pihat_dat$FID1[j],
                             pihat_dat$FID2[j])
  
  IID_to_rm[j] = ifelse(pihat_dat$Z0[j] < pihat_dat$Z1[j], pihat_dat$IID1[j],
                             pihat_dat$IID2[j])
}

inds_to_rm = cbind(FID=FID_to_rm, IID=IID_to_rm)

# Trying to remove duplicates
(inds_to_rm = unique.data.frame(inds_to_rm))

write.table(inds_to_rm, "0.2_low_call_rate_pihat.txt", col.names = TRUE, row.names = FALSE,
            quote = FALSE)






###### Start of future improvements thoughts #########

# # Saving output as 0.2_low_call_rate_pihat.txt 
# # This is a tab separated file with column headers FID IID
# 
# # Testing with randomly imputed values for these two individuals 20 times to see if
# # It runs and works after vectorization and not just 1 value... 
# 
# Z0 = c()
# Z1 = c()
# IID1 = c()
# IID2 = c()
# FID1 = c()
# FID2 = c()
# set.seed(1)
# 
# for (i in 1:rdunif(1, 10, 100)) {
#   Z0[i] = runif(1, min = 0, max = 1)
#   Z1[i] = runif(1, min = 0, max = 1)
#   IID1[i] = pihat_dat$IID1
#   IID2[i] = pihat_dat$IID2
#   FID1[i] = pihat_dat$FID1
#   FID2[i] = pihat_dat$FID2
# }
# 
# test_1 = data.frame(FID1, IID1, FID2, IID2, Z0, Z1) # 
# 
# for (j in 1:nrow(test_1)) {
#   FID_to_rm_test[j] = ifelse(test_1$Z0[j] < test_1$Z1[j], test_1$FID1[j],
#                       test_1$FID2[j])
# 
#   IID_to_rm_test[j] = ifelse(test_1$Z0[j] < test_1$Z1[j], test_1$IID1[j],
#                       test_1$IID2[j])
# }
# 
# (inds_to_rm_test = cbind(FID=FID_to_rm_test, IID=IID_to_rm_test))
# 
# write.table(inds_to_rm_test, "0.2_low_call_rate_pihat_test.txt", col.names = TRUE, row.names = FALSE,
#             quote = FALSE)
# 


# # For a later model/update ... 
# Will need to figure out how to set the working directory to the current location
# of the person working on it
#
# # Making these as functions to be used in tapply... 
# FID_to_rm_func = function(dat1, j) {
#   FID_being_removed = ifelse(dat1$Z0[j] < dat1$Z1[j], dat1$FID1[j],
#                              dat1$FID2[j])
#   return(FID_being_removed)
# }
# 
# 
# IID_to_rm_func = function(dat1, j) {
#   IID_being_removed = ifelse(dat1$Z0[j] < dat1$Z1[j], dat1$IID1[j],
#                              dat1$IID2[j])
#   return(IID_being_removed)
# }
# 
# FID_to_rm_func(test_1, 1) # It works!
# IID_to_rm_func(test_1, 1) # It works!
# 
# 
# # Trying tapply...
# k=1:nrow(test_1)
# tapply(X= k, FID_to_RM[k]=FID_to_rm_func(test_1, k), 
#        IID_to_RM[k] = IID_to_rm_func(test_1, k))
# 
# (inds_to_rm_tapply = cbind(FID=FID_to_RM, IID=IID_to_RM))