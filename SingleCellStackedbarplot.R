# Stacked Barplot as a function


#’ Stacked proportion barplot for any two categorical metadata in a Seurat object
#’
#’ @param seurat_obj A Seurat object
#’ @param Celltype_var Unquoted name of the metadata column to use for fill (e.g. Celltype)
#’ @param samples_var   Unquoted name of the metadata column to use on the x-axis (e.g. orig.ident)

stacked_bar_proportion <- function(seurat_obj,
                                   Celltype_var,
                                   samples_var,
                                   width = 0.4) {
  # capture user’s arguments
  Celltype_var <- rlang::ensym(Celltype_var)
  samples_var    <- rlang::ensym(samples_var)
  
  # 1) build counts + proportions
  df_counts <- seurat_obj@meta.data %>%
    dplyr::count(!!Celltype_var, !!samples_var) %>%
    dplyr::group_by(!!samples_var) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    as.data.frame()
  
  # 2) extract the same colours that Seurat would use for `Celltype_var`
  p_dim <- DimPlot(seurat_obj, group.by = rlang::as_string(Celltype_var))
  gb <- ggplot2::ggplot_build(p_dim)
  pd <- gb$data[[1]]
  gc <- unique(pd[, c("group","colour")])
  levels_fill <- levels(seurat_obj[[rlang::as_string(Celltype_var)]][,1])
  gc$Ident <- levels_fill[gc$group]
  col_map <- setNames(as.character(gc$colour), gc$Ident)
  
  # 3) make the plot
  library(ggplot2)
  p <- ggplot(df_counts, aes(x = !!samples_var, y = prop, fill = !!Celltype_var)) +
    geom_bar(stat = "identity",
             position = "fill",
             width = width,
             color = "black") +
    scale_fill_manual(values = col_map) +
    scale_y_reverse(
      labels = function(x) paste0(round((1 - x) * 100), "%"),
      breaks = seq(0, 1, 0.2)
    ) +
    theme(
      axis.title    = element_blank(),
      axis.text.x   = element_text(size = rel(1.2), face = "bold", color = "black"),
      axis.text.y   = element_text(size = rel(1.2), face = "bold", color = "black"),
      panel.background = element_blank(),
      axis.ticks     = element_blank(),
      legend.title   = element_blank(),
      text           = element_text(face = "bold", colour = "black")
    ) +
    guides(fill = guide_legend(ncol = 1))
  
  return(p)
}
