
#
# cl <- umap_clustering(ascn,
#                       minPts = max(round(0.05 * length(unique(ascn$cell_id))), 2))
# ascn_switch <- rephase_alleles(ascn, minfrac = 0.05)
#
# pdf("~/Downloads/heatmap_2.pdf", width = 25)
# plotHeatmap(ascn, clusters = cl$clustering,
#             plotcol = "state", plottree = F, reorderclusters = T)
# plotHeatmap(ascn, clusters = cl$clustering,
#             plotcol = "state_phase", plottree = F, reorderclusters = T)
# plotHeatmap(ascn_switch, clusters = cl$clustering,
#             plotcol = "state_phase", plottree = F, reorderclusters = T)
# dev.off()
#
ascn_switch1 <- rephase_alleles(ascn, minfrac = 0.025, field = "Min")
ascn_switch2 <- rephase_alleles(ascn, method = "loh")

dim(create_segments(ascn, state_phase))
dim(create_segments(ascn_switch1, state_phase))
dim(create_segments(ascn_switch2, state_phase))

pdf("~/Downloads/heatmap_5.pdf", width = 25)
plotHeatmap(ascn, clusters = cl$clustering,
            plotcol = "state", plottree = F, reorderclusters = T)
plotHeatmap(ascn, clusters = cl$clustering,
            plotcol = "state_phase", plottree = F, reorderclusters = T)
plotHeatmap(ascn_switch1, clusters = cl$clustering,
            plotcol = "state_phase", plottree = F, reorderclusters = T)
plotHeatmap(ascn_switch2, clusters = cl$clustering,
            plotcol = "state_phase", plottree = F, reorderclusters = T)
dev.off()
#


