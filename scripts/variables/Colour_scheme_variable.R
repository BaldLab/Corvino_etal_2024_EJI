####################
# Colour scheme
####################
sample.cols <- c("lightcyan", "lightblue2", "lightsalmon", "lightsalmon3")

condition.cols <- c("royalblue", "salmon2")
names(condition.cols) <- c("US", "Stimulated")

clust.cols <- c("#E41A1C", # Naive_like_1_CM
                "#A6761D", # Naive_like_2_SC
                "#8DA0CB", # Naive_like_3
                "#666666", # Cytotoxic
                "#A6D854", # Type_I_IFN
                "#984EA3", # Stimulated_1
                "#1B9E77", # Stimulated_exhausted
                "#D95F02", # Exhausted_1
                "#7570B3", # Exhausted_2
                "#E7298A", # TRM
                "#E6AB02", # gd_T_g9d2
                "#8DD3C7", # gd_T_non_g9d2
                "#FF7F00", # MAIT 
                "#E78AC3") # Proliferative

clust.names <- c("Naive_like_1_CM",
                 "Naive_like_2_SC",
                 "Naive_like_3",
                 "Cytotoxic",
                 "Type_I_IFN",
                 "Stimulated_1",
                 "Stimulated_exhausted",
                 "Exhausted_1",
                 "Exhausted_2",
                 "TRM",
                 "gd_T_g9d2",
                 "gd_T_non_g9d2",
                 "MAIT",
                 "Proliferative")

names(clust.cols) <- clust.names


# For visualisation use batlow where possible

#install.packages('scico')
#library("scico")

#scico_palette_show()
#batlow.pal <- scico(100, palette = 'batlow')
usefulfunctions::Get.batlow()



# Custom color scheme for clusters
keep.logic <- names(clust.cols) %in% c("gd_T_g9d2", "gd_T_non_g9d2", "MAIT")
innate.cols <- clust.cols[keep.logic]
CD8.cols <- clust.cols[!keep.logic]

# remove un-necessary variables
rm(keep.logic, clust.names)
