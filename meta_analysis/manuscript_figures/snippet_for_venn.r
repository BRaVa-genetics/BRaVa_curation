# Create a venn diagram
library(venneuler)

# Define your counts
# For example, A has 10, B has 12, and intersection A∩B has 5
venn_data <- venneuler(c("BRaVa" = 34, "Jurgens" = 6, "BRaVa&Jurgens" = 31))

# Plot the Venn diagram
pdf("Venn.pdf")
plot(venn_data)
dev.off()
