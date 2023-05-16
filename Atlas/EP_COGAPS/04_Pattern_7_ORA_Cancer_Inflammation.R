# Author: Jacob Mitchell

# Waterfall plot of pattern 7 hallmarks colored by association
# with cancer progression or inflammation

library(ggplot2)
library(dplyr)
library(forcats)

result <- read.csv("results/MsigDB_ORA/Pattern_7_MsigDB_OverRepresentationAnalysis.csv")[,2:10]

# Limit result to q-values <0.05
result <- result[result$padj < 0.05,]

# Reorder dataframe by ascending adjusted p value
result <- mutate(result, MsigDB_Hallmark=fct_reorder(MsigDB_Hallmark, - padj))
# Add column for hallmark groups for colors
result$Hallmark_Group <- c("inflammation","cancer","inflammation","inflammation",
                           "inflammation","cancer","inflammation","inflammation",
                           "cancer","inflammation","inflammation","cancer",
                           "cancer","cancer")

plot <- ggplot(result, aes_string(y = "neg.log.q", x = "MsigDB_Hallmark", fill = "Hallmark_Group")) +
  ## Specifies barplot
  geom_col() +
  ## Manually set colors to highlight cancer hallmarks
  scale_fill_manual(values=list("cancer" = "red",
                                "inflammation" = "orange",
                                "other" = "gray")) +
  ## Rename y axis
  ylab("-10*log10(FDR q-value)") +
  ## Flips the coordinates
  coord_flip() +
  ## Makes the background white
  theme_minimal() +
  ## Add title
  ggtitle("Pattern 7 Molecular Signature Hallmarks") +
  ## This creates the dotted line at .05 value 
  geom_hline(yintercept=c(13.0103), linetype="dotted") + # Add veritcle line to show significances
  ## Adds the q values
  geom_text(aes(label=format(signif(padj, 4))), hjust = -.04) +
  ## Rename legend title
  labs(fill = "Hallmark Group") +
  theme(axis.text.y = element_text(size = 8)) +
  ## specifies limits 
  ylim(0, ceiling(max(result$"neg.log.q")) + (max(result$"neg.log.q")/4))
ggsave(paste0("results/figures/Pattern_MsigDB_ORA/Pattern_7_cancer_inflammation_ORA.pdf"),
       plot = plot,
       width = unit(10, "in"),
       height = unit(6, "in"),
       device = "pdf")
