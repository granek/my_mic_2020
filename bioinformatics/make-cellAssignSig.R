# Celltype signature from cellassign pkg (Zhang et al, Nature Methods 2019)

cassignsig <- structure(list(Group = c(
  "B cells", "B cells", "B cells", "B cells",
  "B cells", "B cells", "T cells", "T cells", "T cells", "T cells",
  "T cells", "T cells", "T cells", "T cells", "Cytotoxic T cells",
  "Cytotoxic T cells", "Cytotoxic T cells", "Cytotoxic T cells",
  "Cytotoxic T cells", "Cytotoxic T cells", "Cytotoxic T cells",
  "Cytotoxic T cells", "Cytotoxic T cells", "Cytotoxic T cells",
  "Cytotoxic T cells", "Monocyte/Macrophage", "Monocyte/Macrophage",
  "Monocyte/Macrophage", "Monocyte/Macrophage", "Monocyte/Macrophage",
  "Monocyte/Macrophage", "Monocyte/Macrophage", "Monocyte/Macrophage",
  "Monocyte/Macrophage", "Epithelial cells", "Epithelial cells",
  "Epithelial cells", "Epithelial cells", "Myofibroblast", "Myofibroblast",
  "Myofibroblast", "Myofibroblast", "Myofibroblast", "Vascular smooth muscle cells",
  "Vascular smooth muscle cells", "Vascular smooth muscle cells",
  "Vascular smooth muscle cells", "Vascular smooth muscle cells",
  "Vascular smooth muscle cells", "Vascular smooth muscle cells",
  "Vascular smooth muscle cells", "Vascular smooth muscle cells",
  "Endothelial cells", "Endothelial cells", "Endothelial cells",
  "Endothelial cells", "Endothelial cells", "Endothelial cells",
  "Endothelial cells", "Endothelial cells"
), Marker = c(
  "VIM",
  "MS4A1", "CD79A", "PTPRC", "CD19", "CD24", "VIM", "CD2", "CD3D",
  "CD3E", "CD3G", "CD28", "CD4", "PTPRC", "VIM", "CD2", "CD3D",
  "CD3E", "CD3G", "CD28", "CD8A", "GZMA", "PRF1", "NKG7", "PTPRC",
  "VIM", "CD14", "FCGR3A", "CD33", "ITGAX", "ITGAM", "CD4", "PTPRC",
  "LYZ", "EPCAM", "CDH1", "CLDN3", "CLDN4", "VIM", "ACTA2", "COL1A1",
  "COL3A1", "SERPINH1", "VIM", "ACTA2", "MYH11", "PLN", "MYLK",
  "MCAM", "COL1A1", "COL3A1", "SERPINH1", "VIM", "EMCN", "CLEC14A",
  "CDH5", "PECAM1", "VWF", "MCAM", "SERPINH1"
), Weight = c(
  2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
)), row.names = c(
  NA,
  -60L
), class = "data.frame")

head(cassignsig)

write.csv(cassignsig, file = "/hpc/group/chsi-mic-2022/data/cellAssignSig.csv", row.names = FALSE)

