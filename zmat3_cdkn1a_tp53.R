# ========================================================================================================================================= #
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 05/18/2023
# Description: This script analyzes the Dependency Map CRISPR screening data based on the presence of TP53 mutations in order to determine
# if ZMAT3 and CDKN1A prototypes are conserved pan-cancer
# ========================================================================================================================================= #

# load required packages ----
# Package names
packages <-
  c(
    "scales",
    "ggplot2",
    "readxl",
    "reshape2",
    "plyr",
    "dplyr",
    "data.table",
    "stringr",
    "ggpubr",
    "janitor",
    "tidyverse",
    "magrittr",
    "ggridges",
    "cowplot",
    "ggExtra",
    "wesanderson",
    "ggrepel",
    "tidyverse",
    "viridis",
    "gridExtra",
    "ggblend"
  )

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

`%ni%` <- Negate(`%in%`)
options(scipen = 999)

# make directories for data ----
dir.create("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53")
dir.create(
  "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data"
)
dir.create("~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results")

# ========================= #
# Cell line information ----
# ========================= #
# download cell line information for DepMap lines (e.g. tissue of origin, mutations, etc.)
# cell line mutation file
download.file("https:/ndownloader.figshare.com/files/27902118",
              destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsSomaticMutations.csv")
depmap_lines_mutations <-
  read.csv(
    "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsSomaticMutationsProfile.csv"
  )

# select all TP53 mutations
depmap_lines_mutations_tp53 <- depmap_lines_mutations |>
  subset(HugoSymbol == "TP53") |>
  select(
    VariantType,
    VariantInfo,
    DNAChange,
    ProteinChange,
    HugoSymbol,
    CCLEDeleterious,
    CosmicHotspot,
    ProfileID
  )

# plot the number of cell lines with and without mutations
mut_status <- depmap_lines_mutations |>
  select(ProfileID) |>
  unique() |>
  mutate(
    tp53_status = case_when(
      ProfileID %in% depmap_lines_mutations_tp53$ProfileID ~ "Mut",
      TRUE ~ "WT"
    )
  )

# plot the distribution of mutated vs wt cell lines
p1 <- ggplot(mut_status, aes(x = tp53_status, fill = tp53_status)) +
  geom_bar() +
  scale_fill_manual(values = c("WT" = "#b2182b", "Mut" = "#2166ac")) +
  theme_cowplot() +
  ylab("Number of cell lines") +
  xlab(NULL) +
  guides(fill = guide_legend(title = "TP53 status")) +
  labs(title = "All cell lines in DepMap\n(N = 2,320)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(
  plot = p1,
  filename = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/cell_line_status.pdf",
  dpi = 300,
  width = 5,
  height = 3.5,
  units = "in"
)

# add the cell line identifier using a different file
download.file(url = "https:/figshare.com/ndownloader/files/40449635",
              destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsProfiles.csv")
depmap_cell_line_names = read_csv(
  "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsProfiles.csv"
)

depmap_lines_mutations_tp53 = left_join(depmap_lines_mutations_tp53, depmap_cell_line_names, by = "ProfileID")

# remove unnecessary intermediates
rm(depmap_lines_mutations)

# ============================== #
# Copy number status for p53 ----
# ============================== #
# download the cell line copy number file
download.file(
  "https:/depmap.org/portal/download/all/?releasename=DepMap+Public+23Q2&filename=OmicsCNGene.csv",
  destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsCNGene.csv"
)
depmap_lines_cna_sub <-
  read.csv(
    "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsCNGene.csv"
  ) |>
  select(X, TP53..7157.) |> # only select cell line names and the TP53 copy number info
  set_colnames(c("ModelID", "TP53_CN"))

# now, add the copy number data to the cell line mutation data
depmap_lines_mutations_cn_tp53 = full_join(depmap_lines_mutations_tp53, depmap_lines_cna_sub, by = "ModelID")

# define cases with a deletion in tp53
# No official answer from the DepMap team as far as I know, but I saw that
# Mina et al, 2020 (https:/pubmed.ncbi.nlm.nih.gov/32989323/ 24) used the following cutoffs for CN:
# Amp: CN > 2^0.75
# Del: CN < 2^-1.2
# where CN = 1 means diploid.

# because we don't care about mutation zygosity at this point, we will simplify the cell line stratification
# based on the presence or absence of a TP53 mutation (SNP or copy number)
depmap_lines_mutations_cn_tp53 =  depmap_lines_mutations_cn_tp53 |>
  mutate(
    tp53_mut_del = case_when(
      CCLEDeleterious == "True" |
        CosmicHotspot == "True" |
        TP53_CN < 0.5 | HugoSymbol == "TP53" ~ "Mut",
      TRUE ~ "WT"
    ),
    tp53_mut_del_detailed = case_when(
      CCLEDeleterious == "True" ~ "Deleterious",
      CosmicHotspot == "True" ~ "Hotspot",
      TP53_CN < 0.5 ~ "Deletion",
      HugoSymbol == "TP53" &
        CCLEDeleterious != "True" &
        CosmicHotspot != "True" & TP53_CN >= 0.5 ~ "Other",
      TRUE ~ "WT"
    )
  )

tp53_mut_lines = depmap_lines_mutations_cn_tp53 |>
  subset(tp53_mut_del == "Mut")

# plot distribution of tp53 mutations in the cell lines
cell_line_status <- depmap_lines_mutations_cn_tp53 |>
  select(ModelID, tp53_mut_del_detailed) |>
  unique() |>
  mutate(broad_status = case_when(tp53_mut_del_detailed == "WT" ~ "WT",
                                  TRUE ~ "Mut"))

# plot the distribution of mutation type across the cell liens
ggplot(cell_line_status,
       aes(x = broad_status, fill = tp53_mut_del_detailed)) +
  geom_bar() +
  scale_fill_manual(
    values = c(
      "WT" = "#b2182b",
      "Hotspot" = "#2166ac",
      "Deleterious" = "#a0da39",
      "Deletion" = "#1fa187",
      "Other" = "#277f8e"
    )
  ) +
  theme_cowplot()


# =========================== #
# CRISPR gene effect data ----
# =========================== #
# download the CRONOS gene effect scores from the DepMap database
download.file("https:/figshare.com/ndownloader/files/40448555",
              destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/CRISPRGeneEffect.csv")
# Gene Effect data for each gene in each cell line
depmap_crispr_effect <-
  read.csv(
    "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/CRISPRGeneEffect.csv"
  )

# remove everything but the gene name from the column headers
names(depmap_crispr_effect) <-
  gsub("\\...*", "", names(depmap_crispr_effect))

# make some duplicate gene names unique
colnames(depmap_crispr_effect) <-
  make.unique(colnames(depmap_crispr_effect))

# select zmat3 and cdkn1a phenotype columns and add the cell line mutation status
depmap_crispr_effect_sub <- depmap_crispr_effect |>
  select(ModelID, ZMAT3, CDKN1A) |>
  mutate(tp53_status = case_when(ModelID %in% tp53_mut_lines$ModelID ~ "Mut",
                                 TRUE ~ "WT")) |>
  mutate(tp53_status = fct_relevel(tp53_status, c("WT", "Mut"))) |>
  reshape2::melt() # melt the dataframe in order to plot with facet_wrap()

# plot distribution differences between gene effect scores based on cell line mutation status
p2 <- ggboxplot(
  depmap_crispr_effect_sub,
  x = "tp53_status",
  y = "value",
  color = "tp53_status",
  palette = c("#b2182b", "#2166ac"),
  add = "jitter",
  shape = 19
) +
  stat_compare_means(aes(label = ..p.signif..), size = 3, label.x = 1.5) +
  facet_wrap(~ variable) +
  ylab("Gene Effect") +
  xlab(NULL) +
  guides(color = guide_legend(title = "TP53 status"))

ggsave(
  plot = p2,
  filename = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/pan_cancer/effect_difference_boxplot_pan_cancer.pdf",
  dpi = 300,
  width = 4.5,
  height = 3.5,
  units = "in"
)

# now plot the distribution of mutation types across the cell lines
cell_line_status <- depmap_lines_mutations_cn_tp53 |>
  select(ModelID, tp53_mut_del_detailed) |>
  unique() |>
  subset(ModelID %in% depmap_crispr_effect_sub$ModelID) |>
  mutate(broad_status = case_when(tp53_mut_del_detailed == "WT" ~ "WT",
                                  TRUE ~ "Mut"))

p3 <- ggplot(cell_line_status,
             aes(x = broad_status, fill = tp53_mut_del_detailed)) +
  geom_bar() +
  scale_fill_manual(
    values = c(
      "WT" = "#b2182b",
      "Hotspot" = "#2166ac",
      "Deleterious" = "#a0da39",
      "Deletion" = "#1fa187",
      "Other" = "#277f8e"
    )
  ) +
  theme_cowplot() +
  ylim(0, 1500) +
  ylab("Number of cell lines") +
  xlab(NULL) +
  guides(fill = guide_legend(title = "TP53 status")) +
  labs(title = "Lines with KO scores + TP53 status\n(N = 1,095)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  plto = p3,
  filename = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/cell_line_status_filtered.pdf",
  dpi = 300,
  width = 5,
  height = 3.5,
  units = "in"
)


# ============================= # 
# ZMAT3 and CDKN1A analyses ----
# ============================= # 
# What we want to do is to calculate a differential KO score for all genes by comparing tp53 mut vs. wt lines
# then we want to see where zmat3 and cdkn1a rank in this list
# download cell line information for DepMap lines (e.g. tissue of origin, mutations, etc.)
download.file(url = "https:/figshare.com/ndownloader/files/40448834",
              destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/Model.csv")
depmap_lines <-
  read.csv(
    "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/Model.csv"
  ) |>
  select(
    ModelID,
    StrippedCellLineName,
    DepmapModelType,
    OncotreeSubtype,
    OncotreePrimaryDisease,
    OncotreeLineage
  )

# add a different file that contains the ModelID and ProfileID relationships in order to pair the CRIPR and mutation data
download.file(url = "https:/figshare.com/ndownloader/files/40449635",
              destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsProfiles.csv")
depmap_cell_line_names = read_csv(
  "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsProfiles.csv"
)

depmap_lines_mutations_all = left_join(depmap_lines_mutations, depmap_cell_line_names, by = "ProfileID")
depmap_lines_mutations_all <-
  left_join(depmap_lines_mutations_all, depmap_lines, by = "ModelID")

rm(depmap_lines_mutations)

# subset mutation data to only lines that have paired CRISPR data
depmap_lines_mutations_crispr <- depmap_lines_mutations_all |>
  subset(ModelID %in% depmap_crispr_effect$ModelID)


# p53 target genes ----
tp53_target_genes = read_excel(
  "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/230524_Brooks_p53_target_list.xlsx"
)

# ========================================= #
# Differential CRISPR KO effect analyses ----
# ========================================= #
# Here, we write a function to automate the gene effect difference between mutated and wt cell lines 
# based on cancer type of interest (or pan-cancer)
depmap_effect_enrichment_by_tp53_mutation <-
  function(pan_cancer, cancer_types) {
    if (pan_cancer %in% c("yes", "Yes", "YES")) {
      sub_lines = depmap_lines_mutations_crispr |>
        select(ModelID) |>
        unique()
      
      cancer_types <- "pan_cancer"
    }
    
    # loop through different cancer types if more than one are provided
    for (i in length(cancer_types)) {
      # create a directory to put the results
      dir.create(
        paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          sep = ""
        )
      )
      
      if (pan_cancer %ni% c("yes", "Yes", "YES")) {
        sub_lines = depmap_lines_mutations_crispr |>
          subset(DepmapModelType == cancer_types[i]) |>
          select(ModelID) |>
          unique()
      }
      
      num_lines <- nrow(sub_lines)
      
      # annotate mutation status and order factors
      depmap_crispr_effect_sub = depmap_crispr_effect |>
        subset(ModelID %in% sub_lines$ModelID) |>
        mutate(tp53_status = case_when(ModelID %in% tp53_mut_lines$ModelID ~ "Mut",
                                       TRUE ~ "WT")) |>
        dplyr::mutate(tp53_status = fct_relevel(tp53_status, c("WT", "Mut")))
      
      # calculate the difference in mean effect score between genotype groups
      depmap_crispr_effect_sub_results <-
        sapply(depmap_crispr_effect_sub[,-1], function(x) {
          mean(x[depmap_crispr_effect_sub$tp53_status == "WT"]) - mean(x[depmap_crispr_effect_sub$tp53_status == "Mut"])
        }) |>
        as.data.frame()
      
      depmap_crispr_effect_sub_results$Gene = rownames(depmap_crispr_effect_sub_results)
      
      colnames(depmap_crispr_effect_sub_results) <-
        c("mean_difference", "Gene")
      
      # remove everything but the gene name from the gene names
      depmap_crispr_effect_sub_results$Gene <-
        gsub("\\...*", "", depmap_crispr_effect_sub_results$Gene)
      
      # calculate p-value from t-test and add this to the results
      exclude_columns <- c(1, 17933)
      p_value <-
        sapply(depmap_crispr_effect_sub[,-exclude_columns], function(x) {
          t.test(x = x[depmap_crispr_effect_sub$tp53_status == "WT"],
                 y = x[depmap_crispr_effect_sub$tp53_status == "Mut"])$p.value
        })
      
      p_value <- as.data.frame(p_value)
      p_value$Gene = rownames(p_value)
      
      # add the p-value to the mean differences
      depmap_crispr_effect_sub_results <-
        left_join(depmap_crispr_effect_sub_results, p_value, by = "Gene")
      depmap_crispr_effect_sub_results$fdr <-
        p.adjust(depmap_crispr_effect_sub_results$p_value)
      
      # Order rows by decreasing value and add order number as a column
      depmap_crispr_effect_sub_results <-
        depmap_crispr_effect_sub_results |>
        na.omit()
      top_genes <-
        as.numeric((nrow(depmap_crispr_effect_sub_results) - 5))
      
      depmap_crispr_effect_sub_results <-
        depmap_crispr_effect_sub_results |>
        arrange(desc(mean_difference)) |>
        mutate(
          order_number = row_number(),
          Gene_name = case_when(
            Gene %in% c("ZMAT3", "CDKN1A") ~ Gene,
            # annotate top and bottom five genes as well as target genes
            order_number < 6 ~ Gene,
            order_number > top_genes ~ Gene
          )
        )
      
      # waterfall plot of mean effect differences between WT and Mut for each gene
      ggplot(
        depmap_crispr_effect_sub_results,
        aes(
          color = mean_difference,
          reorder(mean_difference,-order_number),
          mean_difference,
          label = Gene_name
        )
      ) +
        geom_point() +
        geom_bar(aes(fill = mean_difference), stat = "identity") +
        scale_color_viridis(name = "Effect\ndifference") +
        scale_fill_viridis() +
        geom_label_repel() +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        xlab("Rank order") +
        ylab("Effect difference (WT-Mut)") +
        annotate(
          "segment",
          x = 8601,
          y = -0.6,
          xend = 16500,
          yend = -0.6,
          size = 5,
          linejoin = "mitre",
          arrow = arrow(type = "closed", length = unit(0.01, "npc")),
          color = "#2166ac"
        ) +
        annotate(
          "text",
          x = 12800.5,
          y = -0.6,
          label = "Gene is more essential in Mut",
          color = "white",
          size = 10 / .pt
        ) +
        annotate(
          "segment",
          x = 8601,
          y = 0.6,
          xend = 500,
          yend = 0.6,
          size = 5,
          linejoin = "mitre",
          arrow = arrow(type = "closed", length = unit(0.01, "npc")),
          color = "#b2182b"
        ) +
        annotate(
          "text",
          x = 4550.5,
          y = 0.6,
          label = "Gene is more essential in WT",
          color = "white",
          size = 10 / .pt
        ) +
        ggtitle(paste(cancer_types[i], "\n(", num_lines, " lines)", sep = "")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        guides(fill = FALSE)
      
      ggsave(
        filename = paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          "/",
          cancer_types[i],
          "_effect_difference_waterfall.pdf",
          sep = ""
        ),
        dpi = 300,
        width = 6,
        height = 3.5,
        units = "in"
      )
      write_csv(
        x = depmap_crispr_effect_sub_results,
        file = paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          "/",
          cancer_types[i],
          "_effect_difference_data.csv",
          sep = ""
        )
      )
      
      
      # volcano plot of mean vs p-val
      # plot the results as a volcano plot
      # Calculate quantiles
      quantiles <-
        quantile(
          depmap_crispr_effect_sub_results$mean_difference,
          probs = c(0.1, 0.9)
        )
      
      # Create a new column 'annotation' based on quantiles
      depmap_crispr_effect_sub_results$annotation <-
        ifelse(
          depmap_crispr_effect_sub_results$mean_difference < quantiles[1],
          "10%",
          ifelse(
            depmap_crispr_effect_sub_results$mean_difference > quantiles[2],
            "10%",
            "Middle 80%"
          )
        )
      
      # add column for alpha
      depmap_crispr_effect_sub_results <-
        depmap_crispr_effect_sub_results |>
        mutate(
          Significance = case_when(
            mean_difference < .1 &
              mean_difference > -.1 ~ "Middle 80%",
            annotation == "10%" &
              fdr < 0.05 ~ "10th percentile",
            annotation == "10%" &
              fdr < 0.05 ~ "10th percentile",
            TRUE ~ "Middle 80%"
          )
        )

   # plot the results as a volcano plot
   depmap_crispr_effect_sub_results$Significance <- factor(depmap_crispr_effect_sub_results$Significance, levels = c("10th percentile", "Middle 80%"))
   
   ggplot(depmap_crispr_effect_sub_results, aes(mean_difference,-log(fdr))) +
     geom_vline(
       xintercept = -.1,
       slope = (0),
       color = "lightgrey",
       linetype = "dashed",
       size = .5
     ) +
     geom_vline(
       xintercept = .1,
       slope = (0),
       color = "lightgrey",
       linetype = "dashed",
       size = .5
     ) +
     geom_abline(
       intercept = -log(0.05),
       slope = (0),
       color = "lightgrey",
       linetype = "dashed",
       size = .5
     ) +
     geom_point(aes(fill = mean_difference, alpha = Significance), shape = 21) +
     scale_alpha_discrete(range = c(1, .05), guide = 'none') +
     xlab(label = "Mean effect difference (WT-Mut)") +
     ylab("-log(p-adj)") +
     theme_cowplot() +
     scale_fill_viridis(name = "Effect\ndifference") +
     guides(alpha = guide_legend(override.aes = list(size = 3), title = NULL))
   
   ggsave(
     filename = paste(
       "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
       cancer_types[i],
       "/",
       cancer_types[i],
       "_effect_difference_volcano.pdf",
       sep = ""
     ),
     dpi = 300,
     width = 6,
     height = 3.5,
     units = "in"
   )
   
   
   # annotate tp53 target genes ----
   depmap_crispr_effect_sub_results <- 
     depmap_crispr_effect_sub_results |>
     mutate(label = case_when(
       Gene %in% tp53_target_genes$`Gene Symbol` ~ "Yes",
       TRUE ~ "No"
     ))
   
   depmap_crispr_effect_sub_results$label <- factor(depmap_crispr_effect_sub_results$label, levels = c("Yes", "No"))
   
   depmap_crispr_effect_sub_results <- depmap_crispr_effect_sub_results %>%
     arrange(desc(label == "No"))
   
   # make the volcano plot
   ggplot(depmap_crispr_effect_sub_results, aes(mean_difference,-log(fdr))) +
     geom_vline(
       xintercept = -.1,
       slope = (0),
       color = "lightgrey",
       linetype = "dashed",
       size = .5
     ) +
     geom_vline(
       xintercept = .1,
       slope = (0),
       color = "lightgrey",
       linetype = "dashed",
       size = .5
     ) +
     geom_abline(
       intercept = -log(0.05),
       slope = (0),
       color = "lightgrey",
       linetype = "dashed",
       size = .5
     ) +
     geom_point(aes(fill = mean_difference, alpha = label), shape = 21) +
     scale_alpha_discrete(range = c(1, .05), guide = 'none') +
     xlab(label = "Mean difference") +
     ylab("-log(p-adj)") +
     theme_cowplot() +
     scale_fill_viridis(name = "Effect\ndifference") +
     guides(alpha = guide_legend(order = 1, override.aes = list(size = 3), title = "p53 target gene"))
   
   ggsave(
     filename = paste(
       "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
       cancer_types[i],
       "/",
       cancer_types[i],
       "_effect_difference_volcano_tp53_targets.pdf",
       sep = ""
     ),
     dpi = 300,
     width = 6,
     height = 3.5,
     units = "in"
   )
    
      
      # plot cases where the effect difference between WT and Mut fall in ++ (tumor suppressor activity) and -- (synthetic lethal like) cases
      # Create a new dataframe with mean values based on mutation status
      mean_values <-
        aggregate(. ~ tp53_status, data = depmap_crispr_effect_sub[, -1], FUN = mean)
      
      # Extract the relevant columns from the mean_values dataframe
      mean_values_t <- mean_values |>
        t() |>
        as.data.frame() |>
        row_to_names(1) |>
        rownames_to_column(var = "Gene")
      
      # Add the mean differences to the mean_values dataframe
      depmap_crispr_effect_sub_results <-
        left_join(depmap_crispr_effect_sub_results, mean_values_t, by = "Gene")
      
      # now, plot the effect differences faceted by the type of effect difference group
      depmap_crispr_effect_sub_results <-
        depmap_crispr_effect_sub_results |>
        mutate(
          WT = as.numeric(WT),
          Mut = as.numeric(Mut),
          group = case_when(WT > 0 & Mut > 0 ~ "+ means",
                            WT < 0 & Mut < 0 ~ "- means",
                            TRUE ~ "Other")
        ) |>
        group_by(group) |>
        arrange(desc(mean_difference)) |>
        mutate(
          order_number = row_number(),
          Gene_name = case_when(order_number <= 5 ~ Gene,
                                order_number > n() - 5 ~ Gene)
        ) |>
        ungroup()
      
      ggplot(
        depmap_crispr_effect_sub_results,
        aes(
          color = mean_difference,
          reorder(
            mean_difference,
            mean_difference
          ),
          mean_difference,
          label = Gene_name
        )
      ) +
        geom_point(size = .5) +
        geom_bar(aes(fill = mean_difference), stat = "identity") +
        scale_color_viridis(name = "Effect\ndifference") +
        scale_fill_viridis() +
        geom_label_repel() +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        xlab("Gene rank order\n<- More essential in WT | More essential in Mut ->") +
        ylab("Effect difference (WT-Mut)") +
        ggtitle(paste(cancer_types, "\n(", num_lines, " lines)", sep = "")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        guides(fill = FALSE) +
        facet_wrap(~ group, nrow = 3)
      
      ggsave(
        filename = paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          "/",
          cancer_types[i],
          "_effect_difference_waterfall_facet_mean_group.pdf",
          sep = ""
        ),
        dpi = 300,
        width = 6,
        height = 8,
        units = "in"
      )
      write_csv(
        x = depmap_crispr_effect_sub_results,
        file = paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          "/",
          cancer_types[i],
          "_effect_difference_facet_mean_group_data.csv",
          sep = ""
        )
      )
      
      # now, perform a targeted analysis only within known TP53 target genes
      target_genes <- tp53_target_genes
      
      depmap_crispr_effect_sub_tp53_targets <-
        depmap_crispr_effect_sub_results |>
        subset(Gene %in% target_genes$`Gene Symbol`) |>
        arrange(desc(mean_difference)) |>
        mutate(
          order_number = row_number(),
          Gene_name = case_when(
            Gene %in% c("ZMAT3", "CDKN1A") ~ Gene,
            order_number <= 5 ~ Gene,
            order_number > n() - 5 ~ Gene
          )
        )
      
      ggplot(
        depmap_crispr_effect_sub_tp53_targets,
        aes(
          color = mean_difference,
          reorder(
            mean_difference,
            mean_difference
          ),
          mean_difference,
          label = Gene_name
        )
      ) +
        geom_point(size = 0.5) +
        geom_bar(aes(fill = mean_difference), stat = "identity") +
        scale_color_viridis(name = "Effect\ndifference") +
        scale_fill_viridis() +
        geom_label_repel() +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        xlab("Rank order") +
        ylab("Effect difference (WT-Mut)") +
        ggtitle(paste(cancer_types[i], "\nTP53 target genes")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        guides(fill = FALSE)
      
      ggsave(
        filename = paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          "/",
          cancer_types[i],
          "_effect_difference_waterfall_tp53_target_genes.pdf",
          sep = ""
        ),
        dpi = 300,
        width = 6,
        height = 3.5,
        units = "in"
      )
      write_csv(
        depmap_crispr_effect_sub_tp53_targets,
        file = paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          "/",
          cancer_types[i],
          "_tp53_target_genes_effect_difference_data.csv",
          sep = ""
        )
      )
      
      # plot cases where the effect difference between WT and Mut fall in ++ (tumor suppressor activity) and -- (synthetic lethal like) case
      depmap_crispr_effect_sub_tp53_targets <-
        depmap_crispr_effect_sub_results |>
        # select(depmap_crispr_effect_sub_results, Gene, order_number) |>
        subset(Gene %in% target_genes$`Gene Symbol`) |>
        group_by(group) |>
        arrange(desc(depmap_crispr_effect_sub_tp53_targets)) |>
        mutate(
          order_number = row_number(),
          Gene_name = case_when(
            Gene %in% c("ZMAT3", "CDKN1A") ~ Gene,
            order_number <= 5 ~ Gene,
            order_number > n() - 5 ~ Gene
          )
        ) |>
        ungroup()
      
      ggplot(
        depmap_crispr_effect_sub_tp53_targets,
        aes(
          color = mean_difference,
          reorder(
            mean_difference,
            mean_difference
          ),
          mean_difference,
          label = Gene_name
        )
      ) +
        geom_point(size = 0.5) +
        geom_bar(aes(fill = mean_difference), stat = "identity") +
        scale_color_viridis(name = "Effect\ndifference") +
        scale_fill_viridis() +
        geom_label_repel() +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        xlab("Rank order") +
        ylab("Effect difference (WT-Mut)") +
        ggtitle(paste(cancer_types[i], "\nTP53 target genes")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        guides(fill = FALSE) +
        facet_wrap(~ group, nrow = 3)
      
      ggsave(
        filename = paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          "/",
          cancer_types[i],
          "_effect_difference_waterfall_facet_mean_group_tp53_targets.pdf",
          sep = ""
        ),
        dpi = 300,
        width = 6,
        height = 10,
        units = "in"
      )
      write_csv(
        x = depmap_crispr_effect_sub_tp53_targets,
        file = paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_types[i],
          "/",
          cancer_types[i],
          "_effect_difference_facet_mean_group_tp53_targets_data.csv",
          sep = ""
        )
      )
      
      
    }
  }

depmap_effect_enrichment_by_tp53_mutation(pan_cancer = "Yes", cancer_types = "HCC")



# ============================================ #
# Differential Expression and GSEA analyses ####
# ============================================ #
# We want to perform differential gene expression analysis between TP53 mutated and wild type groups in order to determine how important ZMAT3 and CDKN1A are compared to other genes
# We want to perform the analysis in a pan-cancer and individual cancer types
# load required libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#### libraries and data for DESeq2 ####
BiocManager::install(c("KEGGREST", "DESeq2"))
library("KEGGREST")
library("DESeq2")

# download and read in the gene expression data for the cell lines
download.file(url = "https:/figshare.com/ndownloader/files/40449107",
              destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsExpressionGenesExpectedCountProfile.csv")
depmap_rna_count = read_csv(
  "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsExpressionGenesExpectedCountProfile.csv"
)

# remove the white space between HUGO and ENSG from the header
colnames(depmap_rna_count) <-
  sub(" .*", "", colnames(depmap_rna_count))

# add column name for cell line name column
colnames(depmap_rna_count)[1] <- "ProfileID"

# Need to reformat the data for use in DESeq2
depmap_rna_count_t <- as.data.frame(t(depmap_rna_count)) |>
  row_to_names(row_number = 1) |>
  rownames_to_column() |>
  group_by(rowname) |>
  mutate(HugoID = paste(rowname, ".", 1:n(), sep = "")) |>
  ungroup() |>
  mutate(rowname = NULL) |>
  column_to_rownames(var = "HugoID") |>
  mutate_if(is_character, as.integer)

# create metadata for DESeq
metaData = as.data.frame(colnames(depmap_rna_count_t[, 1:ncol(depmap_rna_count_t)]))
names(metaData) = "ProfileID"

# cell line info for labeling
depmap_cell_line_names = read_csv(
  "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/OmicsProfiles.csv"
)

# join meta data with cell line labeling info
metaData_all = left_join(metaData, depmap_cell_line_names, by = "ProfileID")

# annotate TP53 status in meta data file
metaData_all = metaData_all |>
  mutate(tp53_status = case_when(ModelID %in% tp53_mut_lines$ModelID ~ "Mut", # from line 88
                                 TRUE ~ "WT")) |>
  select(ProfileID, tp53_status) |>
  unique() |>
  mutate_if(is.character, as.factor)

# add the cancer type information to the cell lines in order to filter in the function
depmap_lines <-
  read.csv(
    "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/Model.csv"
  ) |>
  select(
    ModelID,
    StrippedCellLineName,
    DepmapModelType,
    OncotreeSubtype,
    OncotreePrimaryDisease,
    OncotreeLineage
  )

cell_line_info <-
  left_join(depmap_lines,  depmap_cell_line_names, by = "ModelID")

#### libraries and data for GSEA ####
BiocManager::install(c("fgsea", "qusage"), force = TRUE)
library(fgsea)
library(qusage)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")

# download and read in files for GSEA analysis
# annotated hallmark pathways from MSigDB
# download.file(url = "https:/www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.1.Hs/h.all.v2023.1.Hs.symbols.gmt", destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/h.all.v2023.1.Hs.symbols.gmt")
pathways.hallmark <-
  gmtPathways(
    "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/h.all.v2023.1.Hs.symbols.gmt.txt"
  )
# transcription factor targets
# download.file(url = "https:/www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.1.Hs/c3.tft.v2023.1.Hs.symbols.gmt", destfile = "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/c3.tft.v2023.1.Hs.symbols.gmt")
pathways.tfs <-
  gmtPathways(
    "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/raw_data/c3.tft.v2023.1.Hs.symbols.gmt.txt"
  )


# ================================= #
# DESeq2 and GSEA by TP53 status ####
# ================================= #
# now that the metadata file is created for all cell lines, write a function to automate subsetting the cell lines and gene
# expression data by cancer type of interest and perform differential expression analysis

depmap_deseq2_gsea_function <-
  function(pan_cancer,
           cancer_type,
           genes_of_interest) {
    # prompt the user to input the abbreviated cancer type symbol
    if (pan_cancer %ni% c("yes", "Yes", "YES")) {
      print(unique(cell_line_info$DepmapModelType))
      cancer_type <- readline(prompt = "Select cancer type: ")
      
      # create a directory for saving results to
      dir.create(
        paste(
          "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
          cancer_type,
          sep = ""
        )
      )
      
      # select the cell lines in the cancer of interest
      metaData_sub <-
        subset(metaData_all,
               cell_line_info$DepmapModelType == cancer_type)
      
      # subset the gene expression dataframe to only lines present in the cancer of interest
      depmap_rna_count_sub <-
        subset(depmap_rna_count, ProfileID %in% metaData_sub$ProfileID)
      
      # some lines in the MetaData file do not have gene expression so filter to only thoes
      metaData_sub <- metaData_sub |>
        subset(ProfileID %in% depmap_rna_count_sub$ProfileID) |>
        droplevels()
    }
    
    if (pan_cancer %in% c("yes", "Yes", "YES")) {
      # make sure the object name is uniform moving forward
      cancer_type <- "pan_cancer"
      metaData_sub <- metaData_all
      depmap_rna_count_sub <- depmap_rna_count
      
      # subset the gene expression dataframe to only lines present in the cancer of interest
      depmap_rna_count_sub <-
        subset(depmap_rna_count, ProfileID %in% metaData_sub$ProfileID)
      
      # some lines in the MetaData file do not have gene expression so filter to only thoes
      metaData_sub <- metaData_sub |>
        subset(ProfileID %in% depmap_rna_count_sub$ProfileID) |>
        droplevels()
    }
    
    # need to make sure the rownames are the profile ID for DESeq2 analysis
    rownames(metaData_sub) <- metaData_sub$ProfileID
    
    # first, plot the distribution of gene expression for the genes of interest based on tp53 status.
    # then, go back to all genes and perform differential expression
    depmap_rna_count_sub = depmap_rna_count_sub |>
      select(ProfileID, ZMAT3, CDKN1A) |>
      mutate(
        tp53_status = case_when(
          ProfileID %in% metaData_sub$ProfileID &
            metaData_sub$tp53_status == "Mut" ~ "Mut",
          # from line 88
          TRUE ~ "WT"
        )
      ) |>
      reshape2::melt() |>
      mutate(tp53_status = fct_relevel(tp53_status, c("WT", "Mut")),
             value = log(value))
    
    # plot the raw difference in expression for the genes of interest
    ggboxplot(
      depmap_rna_count_sub,
      x = "tp53_status",
      y = "value",
      color = "tp53_status",
      palette = c("#b2182b", "#2166ac"),
      add = "jitter"
    ) +
      stat_compare_means(aes(label = ..p.signif..),
                         size = 3,
                         label.x = 1.5) +
      facet_wrap(~ variable) +
      ylab("log(Gene Expression)") +
      xlab(NULL) +
      ggtitle(paste(
        cancer_type,
        " analysis\n",
        "(",
        nrow(metaData_sub),
        " cell lines)",
        sep = ""
      )) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(
      filename = paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/expression_difference_boxplot.pdf",
        sep = ""
      ),
      dpi = 300,
      width = 4.5,
      height = 4,
      units = "in"
    )
    
    # DESeq2 ----
    # now perform the differential expression
    # need to subset the formatted gene expression file
    depmap_rna_count_deseq <- depmap_rna_count_t |>
      select(which(colnames(depmap_rna_count_t) %in% metaData_sub$ProfileID))
    
    # make the DESeq2 object
    dds_tp53 <-
      DESeq2::DESeqDataSetFromMatrix(
        countData = depmap_rna_count_deseq,
        colData = metaData_sub ,
        design =  ~ tp53_status,
        tidy = F
      )
    
    keep <- rowSums(DESeq2::counts(dds_tp53) >= 10) >= 5
    dds_tp53 <- dds_tp53[keep,]
    dds_tp53$tp53_status <-
      relevel(dds_tp53$tp53_status, ref = "WT")
    dds_tp53 <- DESeq2::DESeq(dds_tp53)
    
    res_tp53 <- DESeq2::results(dds_tp53, tidy = T)
    
    # remove the ".1" annotation
    res_tp53$row = str_sub(res_tp53$row, end = -3)
    
    res_tp53 = res_tp53 |>
      mutate(
        dir_effect = case_when(
          log2FoldChange > 1.5 & padj < 0.01 ~ "upregulated",
          log2FoldChange < -1.5 & padj < 0.01 ~ "downregulated",
          TRUE ~ "NS"
        ),
        fold_change = case_when(
          log2FoldChange < 1.5 & log2FoldChange > -1.5 ~ "small_fc_ns",
          log2FoldChange > 0 & padj < 0.01 ~ "large_fc_sig",
          log2FoldChange < 0 & padj < 0.01 ~ "large_fc_sig",
          TRUE ~ "large_fc_ns"
        ),
        labels = case_when(row %in% c(genes_of_interest) ~ row, # specify the genes of interest
                           # row %in% c("ZMAT3", "CDKN1A", "MDM2") ~ row, # specify the genes of interest
                           TRUE ~ "")
      )
    
    # # select top 10 differentially expressed genes
    # res_tp53$padj <- as.numeric(res_tp53$padj)
    #
    # genes_to_plot <- res_tp53 |>
    #   subset(dir_effect != "NS") |>
    #   group_by(dir_effect) |>
    #   arrange(-padj) |>
    #   subset(padj < 0.01) |>
    #   top_n(10, row) |>
    #   ungroup()
    #
    #   res_tp53 <- res_tp53 |>
    #     mutate(labels = case_when(
    #       row %in% genes_to_plot$row ~ row,
    #       TRUE ~ labels
    #     ))
    
    res_tp53$fold_change = factor(res_tp53$fold_change,
                                  levels = c("large_fc_sig", "large_fc_ns", "small_fc_ns"))
    
    # plot the results as a volcano plot
    ggplot(res_tp53, aes(log2FoldChange,-log(padj))) +
      geom_vline(
        xintercept = -1.5,
        slope = (0),
        color = "#969696",
        linetype = "dashed",
        size = .5
      ) +
      geom_vline(
        xintercept = 1.5,
        slope = (0),
        color = "#969696",
        linetype = "dashed",
        size = .5
      ) +
      geom_abline(
        intercept = -log(0.01),
        slope = (0),
        color = "#969696",
        linetype = "dashed",
        size = .5
      ) +
      geom_point(aes(fill = dir_effect, alpha = fold_change), shape = 21) +
      scale_fill_manual(values = c(
        "downregulated" = "#b2182b",
        "NS" = "#f5f5f5",
        "upregulated" = "#2166ac"
      )) +
      scale_alpha_discrete(range = c(1, .25), guide = 'none') +
      xlab(label = "log2 FC") +
      ylab("-log(p-adj)") +
      theme_cowplot() +
      ggtitle(paste(cancer_type, "\nTP53 WT vs Mut")) +
      theme(
        legend.title.align = 0.5,
        legend.position = "right",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5)
      ) +
      geom_label_repel(
        aes(label = labels),
        # inherit.aes = T,
        force        = 0.5,
        nudge_y      = 10,
        direction    = "y",
        # hjust        = 1,
        segment.size = 0.5
      ) +
      scale_color_manual(values = c(
        "downregulated" = "#b2182b",
        "NS" = "#f5f5f5",
        "upregulated" = "#2166ac"
      )) +
      guides(fill = guide_legend(override.aes = list(size = 3), title = NULL))
    
    # save plot
    ggsave(
      paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/tp53_differential_expression_custom_volcano.pdf",
        sep = ""
      ),
      width = 7.5,
      height = 5
    )
    # save results file
    write_csv(
      res_tp53,
      file = paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/differential_expression_results.csv",
        sep = ""
      )
    )
    
    # highlight the TP53 target genes
    # plot the results as a volcano plot with tp53 target genes highlighted
    res_tp53 <- 
      res_tp53 |>
      mutate(label = case_when(
        row %in% tp53_target_genes$`Gene Symbol` ~ "Yes",
        TRUE ~ "No"
      )) |>
      arrange(label)
    
    
    ggplot(res_tp53, aes(log2FoldChange,-log(padj))) +
      geom_vline(
        xintercept = -1.5,
        slope = (0),
        color = "#969696",
        linetype = "dashed",
        size = .5
      ) +
      geom_vline(
        xintercept = 1.5,
        slope = (0),
        color = "#969696",
        linetype = "dashed",
        size = .5
      ) +
      geom_abline(
        intercept = -log(0.01),
        slope = (0),
        color = "#969696",
        linetype = "dashed",
        size = .5
      ) +
      geom_point(aes(fill = dir_effect, alpha = label), shape = 21, position = "identity") +
      scale_fill_manual(values = c(
        "Other" = "#b2182b",
        "TP53 target gene" = "#2166ac"
      )) +
      scale_alpha_discrete(range = c(0.05, 1)) +
      xlab(label = "log2 FC") +
      ylab("-log(p-adj)") +
      theme_cowplot() +
      ggtitle(paste(cancer_type, "\nTP53 WT vs Mut")) +
      theme(
        legend.title.align = 0.5,
        legend.position = "right",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_fill_manual(values = c(
        "downregulated" = "#b2182b",
        "NS" = "#f5f5f5",
        "upregulated" = "#2166ac"
      )) +
      guides(fill = guide_legend(override.aes = list(size = 3), title = NULL),
             alpha = guide_legend(override.aes = list(size = 3), title = "TP53 target gene"))
    
    ggsave(
      paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/tp53_differential_expression_volcano_tp53_target_genes.pdf",
        sep = ""
      ),
      width = 7.5,
      height = 5
    )
    
    
    # PCA plot
    #First we need to transform the raw count data
    vsdata_tp53 <-
      DESeq2::vst(dds_tp53, blind = FALSE) #vst function will perform variance stabilizing transformation
    
    pcaData_tp53 <-
      DESeq2::plotPCA(vsdata_tp53, intgroup = "tp53_status", returnData = TRUE)
    percentVar <- round(100 * attr(pcaData_tp53, "percentVar"))
    
    pcaData_tp53 = pcaData_tp53[order(pcaData_tp53$group),]
    
    ggplot(pcaData_tp53, aes(PC1, PC2, fill = tp53_status)) +
      geom_point(size = 3, shape = 21) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      scale_fill_manual(values = c("WT" = "#b2182b", "Mut" = "#2166ac")) +
      theme_cowplot() +
      ggtitle(paste(cancer_type, "\nTP53 WT vs Mut PCA")) +
      theme(
        legend.title.align = 0.5,
        legend.title = element_blank(),
        legend.position = "top",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5)
      )
    ggsave(
      paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/tp53_differential_expression_PCA.pdf",
        sep = ""
      ),
      width = 5,
      height = 5
    )
    
    
    # GSEA ----
    # annotated hallmark pathways from MSigDB
    # randomize the ties
    genesTables_tp53 <- res_tp53 |>
      na.omit() |>
      mutate(rank = rank(log2FoldChange,  ties.method = "random")) |>
      dplyr::arrange(desc(rank))
    
    tp53_Ranks = genesTables_tp53$log2FoldChange
    names(tp53_Ranks) = c(genesTables_tp53$row)
    
    fgseaResult_tp53 <- fgsea(
      pathways = pathways.hallmark,
      stats = tp53_Ranks,
      minSize  = 15,
      maxSize  = 500
    )
    
    topPathwaysUp <-
      fgseaResult_tp53[ES > 0][head(order(pval), n = 10), pathway]
    topPathwaysDown <-
      fgseaResult_tp53[ES < 0][head(order(pval), n = 10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    collapsedPathways <-
      collapsePathways(fgseaResult_tp53[order(pval)][padj < 0.01],
                       pathways.hallmark, tp53_Ranks)
    mainPathways_tp53 <-
      fgseaResult_tp53[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
    
    # plot my own
    plot_data_tp53 = as.data.frame(fgseaResult_tp53) |>
      subset(pathway %in% topPathways)
    
    # options(scipen = 999)
    plot_data_tp53 = plot_data_tp53 |>
      mutate(
        `Enriched in` = case_when(NES < 0 ~ "WT",
                                  NES > 0 ~ "Mut"),
        significance = case_when(padj < 0.1 ~ "padj < 0.1",
                                 padj > 0.1 ~ "padj > 0.1", )
      )
    
    plot_data_tp53$significance = factor(plot_data_tp53$significance,
                                         levels = c("padj > 0.1", "padj < 0.1"))
    
    plot_data_tp53$pathway = str_replace_all(plot_data_tp53$pathway, "HALLMARK_", '')
    plot_data_tp53$pathway = str_replace_all(plot_data_tp53$pathway, "_", ' ')
    
    ggplot(
      plot_data_tp53 |>
        arrange(pathway, NES),
      aes(
        x = NES,
        y = reorder(pathway, NES),
        fill = `Enriched in`,
        alpha = sort(significance, increasing = T)
      )
    ) +
      geom_point(color = "black", shape = 21, aes(size = size)) +
      scale_fill_manual(values = c("WT" = "#b2182b", "Mut" = "#2166ac")) +
      theme_cowplot() +
      ylab(NULL) +
      coord_cartesian(ylim = c(0, 21), clip = "off") +
      geom_vline(
        xintercept = 0,
        color = "#969696",
        linetype = "dashed",
        size = .5
      ) +
      scale_alpha_discrete(range = c(0.5, 1), guide = FALSE) +
      annotate(
        "segment",
        x = 0.1,
        y = 21,
        xend = 1.5,
        yend = 21,
        size = 5,
        linejoin = "mitre",
        arrow = arrow(type = "closed", length = unit(0.01, "npc")),
        color = "#2166ac"
      ) +
      annotate(
        "segment",
        x = -0.1,
        y = 21,
        xend = -1.25,
        yend = 21,
        size = 5,
        linejoin = "mitre",
        arrow = arrow(type = "closed", length = unit(0.01, "npc")),
        color = "#b2182b"
      ) +
      annotate(
        "text",
        x = .75,
        y = 21,
        label = "Mut",
        color = "white"
      ) +
      annotate(
        "text",
        x = -.75,
        y = 21,
        label = "WT",
        color = "white"
      ) +
      guides(
        fill = guide_legend(
          override.aes = list(size = 3),
          title = "Enriched in:",
          order = 2
        ),
        size = guide_legend(order = 1, title = "Set size"),
        alpha = guide_legend(override.aes = list(size = 3), title = "Significance")
      ) +
      ggtitle(paste(cancer_type, "\nGSEA for hallmark pathways")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(
      paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/tp53_gsea2.pdf",
        sep = ""
      ),
      width = 9,
      height = 5
    )
    write_csv(
      plot_data_tp53,
      file = paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/gsea_hallmark_pathways_results.csv",
        sep = ""
      )
    )
    
    
    # transcription factor pathways
    # randomize the ties
    genesTables_tp53 <- res_tp53 |>
      na.omit() |>
      mutate(rank = rank(log2FoldChange,  ties.method = "random")) |>
      dplyr::arrange(desc(rank))
    
    tp53_Ranks = genesTables_tp53$log2FoldChange
    names(tp53_Ranks) = c(genesTables_tp53$row)
    
    fgsea_result_tp53 <- fgsea(
      pathways = pathways.tfs,
      stats = tp53_Ranks,
      minSize  = 15,
      maxSize  = 500,
      nPermSimple = 10000
    )
    
    topPathwaysUp <-
      fgsea_result_tp53[ES > 0][head(order(pval), n = 10), pathway]
    topPathwaysDown <-
      fgsea_result_tp53[ES < 0][head(order(pval), n = 10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    collapsedPathways <-
      collapsePathways(fgsea_result_tp53[order(pval)][padj < 0.01],
                       pathways.tfs, tp53_Ranks)
    mainPathways_tp53 <-
      fgsea_result_tp53[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
    
    # plot my own
    plot_data_tp53_tft = as.data.frame(fgsea_result_tp53) |>
      subset(pathway %in% topPathways)
    
    # options(scipen = 999)
    plot_data_tp53_tft = plot_data_tp53_tft |>
      mutate(
        `Enriched in` = case_when(NES < 0 ~ "WT",
                                  NES > 0 ~ "Mut"),
        significance = case_when(padj < 0.1 ~ "padj < 0.1",
                                 padj > 0.1 ~ "padj > 0.1", )
      )
    
    plot_data_tp53_tft$significance = factor(plot_data_tp53_tft$significance,
                                             levels = c("padj > 0.1", "padj < 0.1"))
    
    plot_data_tp53_tft$pathway = str_replace_all(plot_data_tp53_tft$pathway, "_", ' ')
    
    ggplot(
      plot_data_tp53_tft |>
        arrange(pathway, NES),
      aes(
        x = NES,
        y = reorder(pathway, NES),
        fill = `Enriched in`,
        alpha = sort(significance, increasing = T)
      )
    ) +
      geom_point(color = "black",
                 shape = 21,
                 aes(size = size, alpha = significance)) +
      scale_fill_manual(values = c("WT" = "#b2182b", "Mut" = "#2166ac")) +
      theme_cowplot() +
      ylab(NULL) +
      coord_cartesian(ylim = c(0, 21), clip = "off") +
      geom_vline(
        xintercept = 0,
        color = "#969696",
        linetype = "dashed",
        size = .5
      ) +
      scale_alpha_discrete(range = c(0.5, 1), guide = FALSE) +
      annotate(
        "segment",
        x = 0.1,
        y = 21,
        xend = 1.5,
        yend = 21,
        size = 5,
        linejoin = "mitre",
        arrow = arrow(type = "closed", length = unit(0.01, "npc")),
        color = "#2166ac"
      ) +
      annotate(
        "segment",
        x = -0.1,
        y = 21,
        xend = -1.25,
        yend = 21,
        size = 5,
        linejoin = "mitre",
        arrow = arrow(type = "closed", length = unit(0.01, "npc")),
        color = "#b2182b"
      ) +
      annotate(
        "text",
        x = .75,
        y = 21,
        label = "Mut",
        color = "white"
      ) +
      annotate(
        "text",
        x = -.75,
        y = 21,
        label = "WT",
        color = "white"
      ) +
      guides(
        fill = guide_legend(
          override.aes = list(size = 3),
          title = "Enriched in:",
          order = 2
        ),
        size = guide_legend(order = 1, title = "Set size"),
        alpha = guide_legend(override.aes = list(size = 3), title = "Significance")
      ) +
      ggtitle(paste(cancer_type, "\nGSEA for transcription factor targets")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(
      paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/tp53_gsea2_tf_pathways.pdf",
        sep = ""
      ),
      width = 9,
      height = 5
    )
    write_csv(
      plot_data_tp53_tft,
      file = paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_type,
        "/gsea_transcription_factor_targets_results.csv",
        sep = ""
      )
    )
    
    
    
    
    # CRISPR + DESeq2 ----
    # integrate CRISPR effect difference between Mut and WT lines with the differential expression scores
    # to give one plot showing their relationship
    
    colnames(res_tp53)[colnames(res_tp53) == "row"] <- "Gene"
      depmap_crispr_effect_sub_results
    
    joint_cripsr_deseq2 <- inner_join(res_tp53, depmap_crispr_effect_sub_results, by = "Gene")
    
    # order levels
    joint_cripsr_deseq2$Significance = factor(joint_cripsr_deseq2$Significance,
                                  levels = c("Middle 80%", "10th percentile"))
    
    # highlight significant effect differences
    ggplot(joint_cripsr_deseq2, aes(log2FoldChange, mean_difference)) +
      geom_vline(
        xintercept = -1.5,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_vline(
        xintercept = 1.5,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_abline(
        intercept = 0.1,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_abline(
        intercept = -0.1,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_point(aes(fill = mean_difference, alpha = Significance), shape = 21) +
      scale_alpha_discrete(range = c(0.05, 1)) +
      xlab(label = "log2 Fold Change (WT vs. Mut)") +
      ylab("CRISPR effect difference (WT-Mut)") +
      theme_cowplot() +
      scale_fill_viridis(name = "Effect\ndifference") +
      guides(alpha = guide_legend(override.aes = list(size = 3), title = "Effect score difference")) 
    
    ggsave(
      filename = paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_types[i],
        "/",
        cancer_types[i],
        "_crispr_deseq2_scatterplot_effect_differences_highlighted.pdf",
        sep = ""
      ),
      dpi = 300,
      width = 6,
      height = 5,
      units = "in"
    )

    
    # tp53 targets highlighted
    ggplot(joint_cripsr_deseq2, aes(log2FoldChange, mean_difference)) +
      geom_vline(
        xintercept = -1.5,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_vline(
        xintercept = 1.5,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_abline(
        intercept = 0.1,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_abline(
        intercept = -0.1,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_point(aes(fill = mean_difference, alpha = label.x), shape = 21) +
      scale_alpha_discrete(range = c(0.05, 1)) +
      xlab(label = "log2 Fold Change (WT vs. Mut)") +
      ylab("CRISPR effect difference (WT-Mut)") +
      theme_cowplot() +
      scale_fill_viridis(name = "Effect\ndifference") +
      guides(alpha = guide_legend(override.aes = list(size = 3), title = "TP53 target gene")) 
    
    ggsave(
      filename = paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_types[i],
        "/",
        cancer_types[i],
        "_crispr_deseq2_scatterplot_tp53_targets.pdf",
        sep = ""
      ),
      dpi = 300,
      width = 6,
      height = 5,
      units = "in"
    )
    
    
    mean_groups <- 
      depmap_crispr_effect_sub_results |>
      select(Gene, group)
    
    joint_cripsr_deseq2 <- left_join(joint_cripsr_deseq2, mean_groups, by = "Gene")
    
    # visiualized by mean groups
    # tp53 targets highlighted
    ggplot(joint_cripsr_deseq2, aes(log2FoldChange, mean_difference)) +
      geom_vline(
        xintercept = -1.5,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_vline(
        xintercept = 1.5,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_abline(
        intercept = 0.1,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_abline(
        intercept = -0.1,
        slope = (0),
        color = "lightgrey",
        linetype = "dashed",
        size = .5
      ) +
      geom_point(aes(fill = mean_difference, alpha = label.x), shape = 21) +
      scale_alpha_discrete(range = c(0.05, 1)) +
      xlab(label = "log2 Fold Change (WT vs. Mut)") +
      ylab("CRISPR effect difference (WT-Mut)") +
      theme_cowplot() +
      scale_fill_viridis(name = "Effect\ndifference") +
      facet_wrap(~ group, nrow = 3) +
      guides(alpha = guide_legend(override.aes = list(size = 3), title = "TP53 target gene"))
    
    ggsave(
      filename = paste(
        "~/Library/CloudStorage/Box-Box/Brooks Benard's Files/ZMAT3_CDKN1A_TP53/results/",
        cancer_types[i],
        "/",
        cancer_types[i],
        "_crispr_deseq2_scatterplot_tp53_targets_mean_groups.pdf",
        sep = ""
      ),
      dpi = 300,
      width = 6,
      height = 5,
      units = "in"
    )
    
   }

depmap_deseq2_gsea_function(
  pan_cancer = "Yes",
  cancer_type = "AML",
  genes_of_interest = c("ZMAT3", "CDKN1A")
)