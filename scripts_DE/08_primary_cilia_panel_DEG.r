############### Check for specific primary cilia genes #################

# Load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)


########## Plot genes of interest #############3
gene_info <- read_tsv("/users/genomics/jmartinez/a_primary_cilia_project/99_general/Reiter2017_BigCategories.txt")
head(gene_info)


# 1. Load DESeq2 results
res_p8 <- read.csv("/users/genomics/jmartinez/a_primary_cilia_project/06_fc/ctrl_vs_p8_DESeq2_results1.csv")
res_p15 <- read.csv("/users/genomics/jmartinez/a_primary_cilia_project/06_fc/ctrl_vs_p15_DESeq2_results1.csv")
res_p8_vs_p15 <- read.csv("/users/genomics/jmartinez/a_primary_cilia_project/06_fc/p8_vs_p15_DESeq2_results1.csv")

head(res_p8)
head(res_p15)
head(res_p8_vs_p15)


# 2. Merge DESeq2 results with your gene_info
annotated_genes_p8 <- res_p8 %>%
  inner_join(gene_info, by = c("GeneSymbol" = "Gene"))
annotated_genes_p15 <- res_p15 %>%
  inner_join(gene_info, by = c("GeneSymbol" = "Gene"))
annotated_genes_p8_vs_p15 <- res_p8_vs_p15 %>%
  inner_join(gene_info, by = c("GeneSymbol" = "Gene"))

nrow(res_p8)
nrow(annotated_genes_p8)

nrow(res_p15)
nrow(annotated_genes_p15)

nrow(res_p8_vs_p15)
nrow(annotated_genes_p8_vs_p15)

# 3. Add significance info (optional)
annotated_genes_p8 <- annotated_genes_p8 %>%
  mutate(significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 0.5, "Yes", "No"))

# 4. Volcano Plot with annotation by Function or Group
pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/volcano_ctrl_vs_p8_specific.pdf", width = 10, height = 6)
ggplot(annotated_genes_p8, aes(x = log2FoldChange, y = -log10(padj), color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = GeneSymbol), size = 3.5, max.overlaps = 20) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot of Selected Genes - p8",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value",
       color = "Group") +
  theme(legend.position = "bottom")
dev.off()

######################### Super volcano function ########################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(colorspace)
library(RColorBrewer)

create_volcano_plot_selected_genes_grouped <- function(
    annotated_df,
    title,
    padj_threshold = 0.05,
    log2fc_threshold = 0.5,
    custom_group_colors = NULL,
    point_size = 3.5,
    label_size = 3.2,
    label_max_overlaps = 20,
    max_labels = NULL
) {
    required_cols <- c("GeneSymbol", "log2FoldChange", "padj", "Group")
    if (!all(required_cols %in% names(annotated_df))) {
        stop("Missing required columns: ", paste(setdiff(required_cols, names(annotated_df)), collapse = ", "))
    }

    plot_df <- annotated_df %>%
        filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
        mutate(
            is_significant = padj < padj_threshold,
            Group = factor(Group)
        )

    if (nrow(plot_df) == 0) {
        return(ggplot() + labs(title = title, subtitle = "No data to plot") + theme_void())
    }

    # Prepare labels data frame
    labels_df <- plot_df %>%
        arrange(padj, desc(abs(log2FoldChange)))  # Order by significance and effect size

    if (!is.null(max_labels)) {
        labels_df <- labels_df %>% head(max_labels)  # Limit to max_labels
    }

    unique_groups <- levels(plot_df$Group)

    base_colors <- custom_group_colors
    if (is.null(base_colors)) {
        palette_name <- "Set2"
        base_colors <- setNames(
            if (length(unique_groups) <= brewer.pal.info[palette_name, "maxcolors"]) {
                brewer.pal(length(unique_groups), palette_name)
            } else {
                scales::hue_pal()(length(unique_groups))
            },
            unique_groups
        )
    } else {
        missing_groups <- setdiff(unique_groups, names(base_colors))
        if (length(missing_groups) > 0) {
            base_colors <- c(base_colors, setNames(rep("grey50", length(missing_groups)), missing_groups))
        }
    }

    plot_df <- plot_df %>%
        mutate(
            color = base_colors[as.character(Group)],
            alpha = ifelse(is_significant, 1, 0.3)
        )

    p <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(
            aes(fill = color, alpha = alpha, shape = Group),
            size = point_size
        ) +
        scale_fill_identity() +
        scale_alpha_identity() +
        scale_shape_manual(name = "Group", values = rep(21:25, length.out = length(unique_groups))) +
        geom_text_repel(
            data = labels_df,
            aes(label = GeneSymbol),
            size = label_size, max.overlaps = label_max_overlaps,
            box.padding = 0.35, point.padding = 0.3,
            segment.color = 'grey50', segment.size = 0.3,
            force = 1.5
        ) +
        geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey30", linewidth = 0.6) +
        geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "grey30", linewidth = 0.6) +
        labs(
            title = title,
            x = bquote(~Log[2]~ "Fold Change"),
            y = bquote(-~Log[10]~ "(Adjusted P-value)")
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 11),
            legend.position = "bottom",
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
        ) +
        guides(shape = guide_legend(
            override.aes = list(
                fill = base_colors[levels(plot_df$Group)],
                alpha = 1, size = 4
            )
        ))

    return(p)
}

########### Execution of Volcano ###############

# Define custom group colors
my_group_colors <- c(
  "Ciliary trafficking" = "#E41A1C",
  "Ciliogenesis" = "#4DAF4A",
  "Motile cilium structure" = "#377EB8",
  "Non-motile cilium structure" = "#984EA3"
)

# Create and save plot
pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/volcano_ctrl_vs_p8_specific_enhanced.pdf", width = 10, height = 8)
print(create_volcano_plot_selected_genes_grouped(
  annotated_df = annotated_genes_p8,
  title = "Volcano Plot - ctrl vs p8",
  custom_group_colors = my_group_colors
))
dev.off()
head(annotated_genes_p8, 10)

# Create and save plot
pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/volcano_ctrl_vs_p15_specific_enhanced.pdf", width = 10, height = 8)
print(create_volcano_plot_selected_genes_grouped(
  annotated_df = annotated_genes_p15,
  title = "Volcano Plot - ctrl vs p15",
  custom_group_colors = my_group_colors,
  max_labels = 30
))
dev.off()

# Create and save plot
pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/volcano_p8_vs_p15_specific_enhanced.pdf", width = 10, height = 8)
print(create_volcano_plot_selected_genes_grouped(
  annotated_df = annotated_genes_p8_vs_p15,
  title = "Volcano Plot - p8 vs p15",
  custom_group_colors = my_group_colors
))
dev.off()