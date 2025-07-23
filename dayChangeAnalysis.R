T_cells[['RNA']]$counts[c('CD4', 'CD8A'),]

checkCD <- function(seuObj){
  seuObj_CD = seuObj[["RNA"]]$counts[c("CD4", "CD8A"),]
  seuObj_CD = as.data.frame(t(seuObj_CD))
  
  seuObj_CD <- seuObj_CD %>% mutate(CDLabel = case_when(CD4 > 0 & CD8A == 0 ~ "CD4",
                                                 CD4== 0 & CD8A > 0 ~ "CD8",
                                       .default = "other"))
 
  seuObj_CD
}
t_cell_CD <- checkCD(T_cells)
table(t_cell_CD$CDLabel)

T_cells$cdLabel <- t_cell_CD$CDLabel
T_cells_test = T_cells[["RNA"]]$data[c("CD4", "CD8A"),]

T_cells[[]] %>% mutate(hypoxCAR = paste(hypoxia, CAR, sep = '_')) %>% filter(CAR_pred == 1 & hypoxia == 'HH' & cdLabel != 'other') %>% 
  group_by(day, cdLabel) %>% dplyr::summarise(amount = n()) %>% 
  ggplot(aes(x = day, y = amount, fill = cdLabel))+
  geom_bar(position="dodge", stat="identity")+
  ggtitle('HH')


T_cells[[]] %>% filter(CAR_pred == 1) %>% {table(.$day, .$hypoxia)}
