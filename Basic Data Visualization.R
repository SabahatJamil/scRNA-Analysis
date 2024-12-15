# creating functions
Name <- "Immune"

comb_df <- cell_calculations2(GSE1180286_immune,name = "Immune",path = 
                                "D:/Data_Sets/GSE18/New/Immune/test/",
                              conditions = GSE1180286_immune$new_condition,annotations = GSE1180286_immune$Manual_Annotations)

sbar <- ggplot(comb_df, aes(x = Cell_Type, y = Counts, fill = Conditions)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Cell_types", y = "Counts", fill = "Conditions", title = "Stacked Bar Chart") +
  theme(plot.title = element_text(hjust = 0.5))
sbar 

marker_to_heatmap_function(GSE1180286_immune)

ggplot(df, aes(x = "", y = Counts, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +       
  coord_polar("y") +                             
  theme_void() +                                 
  scale_fill_brewer(palette = "Set3") +          
  labs(title = "Beautiful Pie Chart") +          
  theme(plot.title = element_text(hjust = 0.5,   
                                  size = 16,     
                                  face = "bold"),
        legend.title = element_blank())   


ggplot(df, aes(x = "", y = Counts, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = paste(Name, "Pie Chart")) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_blank()) +
  geom_text(aes(label = paste0(Cell_Type, ": ", format(percentage, nsmall = 1), "%")), 
            position = position_stack(vjust = 0.5), size = 5, color = "white")

ggplot(comb_df, aes(x = "", y = Counts, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = paste(Name, "Pie Chart")) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_blank()) +
  geom_text(aes(label = paste0(Cell_Type, ": ", format(percentage, nsmall = 1), "%")), 
            position = position_stack(vjust = 0.5), size = 5, color = "white")


ggplot(comb_df, aes(x = "", y = Counts, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = paste(Name, "Pie Chart")) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_blank()) +
  geom_text_repel(aes(label = paste0(Cell_Type, ": ", format(percentage, nsmall = 1), "%")), 
                  position = position_stack(vjust = 0.5), size = 5, color = "black", 
                  box.padding = 0.5, nudge_x = 1, nudge_y = 1)

comb_df$percentage

library(ggrepel)

install.packages("plotrix")
library(plotrix)
library(RColorBrewer)

slices <- comb_df$Counts
labels <- comb_df$Category
colors <- brewer.pal(length(slices), "Set3")

pie3D(slices, labels = labels, explode = 0.1, main = "3D Pie Chart from DataFrame",
      col = colors, labelcex = 1.2, start = 0, radius = 1.2, height = 0.2, 
      border = "white", shade = 0.8)

pie3D(slices, 
      labels = labels,         
      explode = 0.1,          
      main = "3D Pie Chart from DataFrame",
      col = colors,            
      labelcex = 1.2,          
      start = 0,               
      radius = 1.2,          
      height = 0.2,           
      border = "white",        
      shade = 0.8)             


library(ggplot2)
library(openxlsx)  # Ensure you have openxlsx package for write.xlsx

cell_calculations3 <- function(seurat_integrated, name, path, conditions, annotations){
  # Initialize empty list to store dataframes
  df_list <- list()
  
  # Get unique conditions
  conditions_count <- unique(conditions)
  
  for(i in 1:length(conditions_count)){
    print(i)
    print(conditions_count[i])
    
    # Set identities and subset data
    Idents(seurat_integrated) <- conditions
    rorpos <- subset(seurat_integrated, idents = conditions[i])
    Idents(rorpos) <- annotations
    
    # Fetch and process data
    n_cells_rorpos <- FetchData(rorpos, vars = "ident") %>%
      dplyr::count(ident) %>%
      tidyr::spread(ident, n)
    
    # Transpose and convert to dataframe
    new <- as.data.frame(t(n_cells_rorpos))
    col <- rownames(new)
    col2 <- new[,1]
    
    new[,1] <- col
    new[,2] <- col2
    new[,3] <- conditions_count[i]
    colnames(new) <- c("Cell_Type", "Counts", "Conditions")
    
    # Calculate percentages
    sumz <- sum(new$Counts)
    new$percentage <- round((new$Counts / sumz) * 100, 1)
    
    # Write dataframe to Excel
    write.xlsx(new, file.path(path, paste0(name, "_", conditions_count[i], "_counts.xlsx")))
    
    # Create pie chart
    p1 <- ggplot(new, aes(x = "", y = Counts, fill = Cell_Type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      theme_void() +
      scale_fill_brewer(palette = "Set3") +
      labs(title = paste(name, "Pie Chart")) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            legend.title = element_blank()) +
      geom_text(aes(label = paste0(Cell_Type, ": ", format(percentage, nsmall = 1), "%")), 
                position = position_stack(vjust = 0.5), size = 5, color = "white")
    
    # Save pie chart
    file_name <- file.path(path, paste0("Pie_Char_of",Name,"Cell_Counts_", conditions_count[i], ".jpeg"))
    ggsave(file_name, plot = p1, width = 8, height = 6, dpi = 300)
    
    # Append to list
    df_list[[i]] <- new
    View(df_list[[i]])
    View(new)
  }
  
  # Combine all dataframes
  combined_df <- do.call(rbind, df_list)
  
  # Write combined dataframe to Excel
  write.xlsx(combined_df, file.path(path, paste0(name, "_Combine_counts.xlsx")))
  
  return(combined_df)
}

comb_df <- cell_calculations3(GSE1180286_immune,name = "Immune",path = 
                                "D:/Data_Sets/GSE18/New/Immune/test/",
                              conditions = GSE1180286_immune$new_condition,annotations = GSE1180286_immune$Manual_Annotations)



##############################################
