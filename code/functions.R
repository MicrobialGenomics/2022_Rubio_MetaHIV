theme_Publication <- function(base_size=14, base_family="helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1), hjust = 0.5),
           plot.subtitle=element_text(
             size = rel(0.7)),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.margin = unit(0, "cm"),
           legend.title = element_blank(),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold", size = rel(0.7))
   ))

}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}


#### Calculate alpha diversity indexes:

alpha_div<-function(d_rpkm, d_counts, metadata){
  shannon<-diversity(d_rpkm[ ,-1], MARGIN = 2, index="shannon")
  simpson<-diversity(d_rpkm[, -1], MARGIN = 2, index="simpson")
  invsimpson<-diversity(d_rpkm[, -1], MARGIN = 2, index="invsimpson")
  nARG<-specnumber(d_rpkm[, -1], MARGIN = 2)
  fisher<-fisher.alpha(d_counts[ ,-1], MARGIN = 2) ##we obtain fisher index from counts (not possible to calculate for rpkm values)
  metadata%>%
    filter(SampleID %in% names(d_rpkm))%>%
    mutate(Shannon=shannon, Simpson=simpson, InvSimpson= invsimpson, nARG=nARG, Fisher=fisher)}


#### Plot a diversity index for a clinical variable

alpha_div_plot<-function(data, clin_var, alpha_var){
  clin_var<-enquo(clin_var)
  alpha_var<-enquo(alpha_var)

  data%>%
    filter(!is.na(!!clin_var))%>%
    ggplot(aes(x=!!clin_var, y=!!alpha_var))+
    geom_boxplot(aes(colour=!!clin_var))+
    geom_jitter(width = 0.2, aes(colour=!!clin_var))+
    stat_compare_means()+
    scale_colour_Publication()+ theme_Publication()+
    theme(legend.position = "none",  axis.text.x = element_text(angle=45, hjust = 1))}


#### Plot all diversity indexes for a clinical variable

alpha_div_plot_all<-function(data, clin_var){
  clin_var<-enquo(clin_var)
  data%>%
    filter(!is.na(!!clin_var))%>%
    pivot_longer(cols=c(Shannon, Simpson, InvSimpson, nARG, Fisher), names_to = "Alpha_Index", values_to = "Alpha_value")%>%
    ggplot(aes(x=!!clin_var, y=Alpha_value))+
    geom_boxplot(aes(colour=!!clin_var))+
    geom_jitter(width = 0.2, aes(colour=!!clin_var))+
    facet_wrap(vars(Alpha_Index), ncol = 3, scales = "free_y")+
    stat_compare_means(label="p.format", label.x = 1, label.y.npc = 0.9)+
    scale_colour_Publication()+ theme_Publication()+
    labs(y="Alpha diversity", caption = "P values calculated by Wilcoxon test", title = clin_var)+
    theme(legend.position = "none", axis.title.x=element_blank(),
          axis.text.x = element_text(angle=45, hjust = 1))}


#### Plot correlation between diversity indexes and quantitative clinical variables

alpha_div_cor_all<-function(data, clin_var){
  clin_var<-enquo(clin_var)
  data%>%
    pivot_longer(cols=c(Shannon, Simpson, InvSimpson, nARG, Fisher), names_to = "Alpha_Index", values_to = "Alpha_value")%>%
    ggplot(aes(x=!!clin_var, y=Alpha_value))+
    geom_point(color='#2980B9', size = 2) +
    geom_smooth(method=lm, color="black")+
    facet_wrap(vars(Alpha_Index), ncol = 3, scales = "free_y")+
    stat_cor(method="spearman", label.x.npc = "left", label.y.npc = 0.9)+
    scale_colour_Publication()+ theme_Publication()+
    labs(y="Alpha diversity", caption = "Spearman correlation", title = clin_var)+
    theme(legend.position = "none", axis.title.x=element_blank(),
          axis.text.x = element_text(angle=45, hjust = 1))}


#### NMDS function (for variables with max 3 levels)

beta_nmds<-function(dist, metadata, clin_var){

  clin_var<-enquo(clin_var)

  dist_df<-as.data.frame(as.matrix(dist)); dist_df$SampleID<-rownames(dist_df)
  meta_dist<-inner_join(metadata, dist_df, by="SampleID")%>%
    filter(!is.na(!!clin_var))

  dist<-meta_dist%>%
    select(all_of(.[["SampleID"]]))%>%
    as.dist()

  clin_var2<-meta_dist%>%pull(!!clin_var) ##clin_var vector for adonis test

  test<-adonis(dist~clin_var2, permutations = 999)
  ptest<-test$aov.tab$`Pr(>F)`[1]
  r2test<-round(test$aov.tab$R2[1],2)

  set.seed(200889)
  nmds <- metaMDS(dist, trace = 0)

  ##scores(nmds) %>%
  nmds$points %>%
    as_tibble(rownames = "SampleID") %>%
    rename(NMDS1=MDS1, NMDS2=MDS2) %>%
    inner_join(., metadata, by="SampleID") %>%
    filter(!is.na(!!clin_var))%>%
    ggplot(aes(x=NMDS1, y=NMDS2, color=!!clin_var, fill=!!clin_var)) +
    stat_ellipse(geom="polygon", show.legend = FALSE, alpha=0.2)+
    geom_point()+
    coord_fixed(ratio = 0.8)+
    labs(caption = glue("ADONIS p={ptest}, r2={r2test}"))+
    scale_colour_Publication()+ theme_Publication()+scale_fill_Publication()+
    theme(legend.text = element_text(size=10))}


#### NMDS function (for variables with 4-8 levels)

beta_nmds2<-function(dist, metadata, clin_var){

  clin_var<-enquo(clin_var)

  dist_df<-as.data.frame(as.matrix(dist)); dist_df$SampleID<-rownames(dist_df)
  meta_dist<-inner_join(metadata, dist_df, by="SampleID")%>%
    filter(!is.na(!!clin_var))

  dist<-meta_dist%>%
    select(all_of(.[["SampleID"]]))%>%
    as.dist()

  clin_var2<-meta_dist%>%pull(!!clin_var) ##clin_var vector for adonis test

  test<-adonis(dist~clin_var2, permutations = 999)
  ptest<-test$aov.tab$`Pr(>F)`[1]
  r2test<-round(test$aov.tab$R2[1],2)

  set.seed(200889)
  nmds <- metaMDS(dist, trace = 0)

  ##scores(nmds) %>%
  nmds$points %>%
    as_tibble(rownames = "SampleID") %>%
    rename(NMDS1=MDS1, NMDS2=MDS2) %>%
    inner_join(., metadata, by="SampleID") %>%
    ggplot(aes(x=NMDS1, y=NMDS2, color=!!clin_var)) +
    stat_ellipse(show.legend = FALSE)+
    geom_point()+
    coord_fixed(ratio = 0.8)+
    labs(subtitle = glue("ADONIS p={ptest}, r2={r2test}"))+
    scale_color_manual(name=quo_name(clin_var),
                       values = c("blue", "red","green4", "orange", "aquamarine4", "magenta2", "gold", "black"))+
    scale_colour_Publication()+ theme_Publication()+
    theme(legend.text = element_text(size=10))}



##Function to obtain metadata variables correlations with NMDs axes
corr_meta_envfit<-function(metadata, nmds){
  nmds_positions<-nmds$points%>%
    as_tibble(rownames="SampleID")%>%
    rename(NMDS1=MDS1, NMDS2=MDS2)

  metadata<-metadata%>%
    filter(SampleID %in% nmds_positions$SampleID)%>%
    arrange(SampleID)

  meta<-metadata%>%select(where(is.numeric))
  fit <-envfit(nmds, meta, perm = 999, na.rm=TRUE)

  ##Generate NMDS positions for the plot
  ##Transform values adjusted by pvalue with fx scores  and multiply by arrowmult   value
  cor<-scores(fit, "vectors")%>%
    as_tibble(rownames="Variable")%>%
    mutate(NMDS1=NMDS1*ordiArrowMul(fit), NMDS2=NMDS2*ordiArrowMul(fit), r2=fit$vectors$r, p.value=fit$vectors$pvals)%>%
    return(cor)
}

##Function to obtain AMR correlations with NMDs axes:
corr_amr_envfit<-function(data, refdata, nmds){
  ##transpose data
  ref_name<-pull(data, ref_name)
  data<-as_tibble(cbind(SampleID = names(data), t(data)))%>%slice(-1)%>%
    mutate_at(vars(-("SampleID")),as.numeric)
  colnames(data)<-c("SampleID", ref_name)

  ##Perform envfit
  fit <-envfit(nmds, data[ ,-1], perm = 999)

  ##Generate NMDS positions for the plot: Transform values adjusted by pvalue with fx scores  and multiply by arrowmult value
  cor<-scores(fit, "vectors")%>%
    as_tibble(rownames="ref_name")%>%
    mutate(NMDS1=NMDS1*ordiArrowMul(fit), NMDS2=NMDS2*ordiArrowMul(fit), r2=fit$vectors$r, p.value=fit$vectors$pvals)%>%
    left_join(., refdata, by="ref_name")
  return(cor)}

##Function to obtain functional group correlations with NMDs axes:
corr_group_envfit<-function(data_group, refdata_group, nmds){
  ##transpose data
  group_name<-pull(data_group, 1)
  data_group<-as_tibble(cbind(SampleID = names(data_group), t(data_group)))%>%slice(-1)%>%
    mutate_at(vars(-("SampleID")),as.numeric)
  colnames(data_group)<-c("SampleID", group_name)

  ##Perform envfit
  fit <-envfit(nmds, data_group[ ,-1], perm = 999)

  refdata_group<-refdata_group%>%rename(group_name=1)

  ##Generate NMDS positions for the plot: Transform values adjusted by pvalue with fx scores  and multiply by arrowmult value
  cor<-scores(fit, "vectors")%>%
    as_tibble(rownames="group_name")%>%
    mutate(NMDS1=NMDS1*ordiArrowMul(fit), NMDS2=NMDS2*ordiArrowMul(fit), r2=fit$vectors$r, p.value=fit$vectors$pvals)%>%
    left_join(., refdata_group, by="group_name")
  return(cor)}

##Function to obtain biplot showing selected correlations (AMR, clinical variables or grouped variables):

biplot_amr_envfit<-function(high_cor, nmds, metadata, clin_var, label_var){
  clin_var<-enquo(clin_var)
  label_var<-enquo(label_var)
  nmds_positions<-nmds$points%>%
    as_tibble(rownames="SampleID")%>%
    rename(NMDS1=MDS1, NMDS2=MDS2)

  nmds_positions%>%
    inner_join(., metadata, by="SampleID")%>%
    ##ggplot(aes(x=NMDS1, y=NMDS2)) +
    ggplot( aes(x=NMDS1, y=NMDS2, color=!!clin_var)) +
    geom_point(alpha=0.5, size=2)+
    geom_segment(data=high_cor,
                 aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
                 arrow = arrow(length = unit(0.2, "cm")), alpha=1,colour="gray",
                 inherit.aes=FALSE)+
    geom_text_repel(data=high_cor,
                    aes(x=NMDS1, y=NMDS2, label=!!label_var),
                    min.segment.length = 0.15, segment.alpha=1, segment.color="gray",
                    inherit.aes=FALSE) +
    # theme_bw()+
    # scale_color_manual(values = c("chartreuse4", "red", "gray"))
    scale_colour_Publication()+ theme_Publication()+scale_fill_Publication()
}



##Get significance of each AMR to a clinical variable (dicotomic)

sig_AMR_clin_dic<-function(data, metadata, refdata, clin_var){

  ##transpose data
  ref_name<-pull(data, ref_name)
  data<-as_tibble(cbind(SampleID = names(data), t(data)))%>%slice(-1)%>%
    mutate_at(vars(-("SampleID")),as.numeric)
  colnames(data)<-c("SampleID", ref_name)

  ##merge data with metadata and reference data
  data_all<-data%>%
    pivot_longer(-SampleID, names_to = "ref_name", values_to = "value")%>%
    inner_join(., metadata, by="SampleID")%>%
    inner_join(., refdata, by="ref_name")

  sig_amr<-data_all%>%
    nest(data=-ref_name)%>%
    mutate(test=map(.x=data, ~wilcox.test(value~!!ensym(clin_var), data=.x)%>%tidy))%>%
    unnest(test)%>%
    mutate(p.adjust=p.adjust(p.value, method = "BH"))%>%
    filter(p.adjust<0.05)%>%
    select(ref_name, p.adjust)

  ##Retain data only from significant amr
  data_all<-data_all%>%
    inner_join(sig_amr, by="ref_name")

  ##get log2 fold change
  log_data<-data_all%>%
    group_by(ref_name, !!enquo(clin_var))%>%
    summarise(mean_counts=mean(value), .groups = 'drop')%>%
    pivot_wider(names_from = quo_name(enquo(clin_var)), values_from = "mean_counts")%>%
    as.data.frame()
  log_data$log2f<-log2(log_data[ ,2])-log2(log_data[ ,3])
  log_data<-log_data%>%select(ref_name, log2f)

  ## get median counts per significant group and clinical variable
  sig_data<-data_all%>%
    group_by(ref_name, !!enquo(clin_var))%>%
    summarise(median_counts=median(value), .groups = 'drop')%>%
    pivot_wider(names_from = quo_name(enquo(clin_var)), values_from = "median_counts")%>%
    inner_join(sig_amr, by="ref_name")%>%
    inner_join(log_data, by="ref_name")%>%
    inner_join(refdata, by="ref_name")

  res<-list(data_all, sig_data)
  names(res)<-c("data_plot", "sig_data")
  return(res)
}


##Get significance of each AMR to a clinical variable (factor)

sig_AMR_clin_factor<-function(data, metadata, refdata, clin_var){

  ##transpose data
  ref_name<-pull(data, 1)
  data<-as_tibble(cbind(SampleID = names(data), t(data)))%>%slice(-1)%>%
    mutate_at(vars(-("SampleID")),as.numeric)
  colnames(data)<-c("SampleID", ref_name)

  ##merge data with metadata and reference data
  data_all<-data%>%
    pivot_longer(-SampleID, names_to = "ref_name", values_to = "value")%>%
    inner_join(., metadata, by="SampleID")%>%
    inner_join(., refdata, by="ref_name")

  sig_amr<-data_all%>%
    nest(data=-ref_name)%>%
    mutate(test=map(.x=data, ~kruskal.test(value~!!ensym(clin_var), data=.x)%>%tidy))%>%
    unnest(test)%>%
    mutate(p.adjust=p.adjust(p.value, method = "BH"))%>%
    filter(p.adjust<0.05)%>%
    select(ref_name, p.adjust)

  data_all<-data_all%>%
    inner_join(sig_amr, by="ref_name")

  sig_amr<-sig_amr%>%
    inner_join(refdata, by="ref_name")

  res<-list(data_all, sig_amr)
  names(res)<-c("data_plot", "sig_amr")
  return(res)
}

##Get significance of a grouped matrix to a clinical variable (dicotomic)

sig_group_clin_dic<-function(data_group, metadata, refdata_group, clin_var){

  ##transpose data
  group_name<-pull(data_group, 1)
  data<-as_tibble(cbind(SampleID = names(data_group), t(data_group)))%>%slice(-1)%>%
    mutate_at(vars(-("SampleID")),as.numeric)
  colnames(data)<-c("SampleID", group_name)

  ##merge data with metadata
  data_all<-data%>%
    pivot_longer(-SampleID, names_to = "group_name", values_to = "value")%>%
    inner_join(., metadata, by="SampleID")

  ##get significant groups, adjusted p value
  sig_amr<-data_all%>%
    nest(data=-group_name)%>%
    mutate(test=map(.x=data, ~wilcox.test(value~!!ensym(clin_var), data=.x)%>%tidy))%>%
    unnest(test)%>%
    mutate(p.adjust=p.adjust(p.value, method = "BH"))%>%
    filter(p.adjust<0.05)%>%
    select(group_name, p.adjust)

  ##obtain data only from significant groups (for the plot)
  data_all<-data_all%>%
    inner_join(sig_amr, by="group_name")

  ##get log2 fold change
  log_data<-data_all%>%
    group_by(group_name, !!enquo(clin_var))%>%
    summarise(mean_counts=mean(value), .groups = 'drop')%>%
    pivot_wider(names_from = quo_name(enquo(clin_var)), values_from = "mean_counts")%>%
    as.data.frame()
  log_data$log2f<-log2(log_data[ ,2])-log2(log_data[ ,3])
  log_data<-log_data%>%select(group_name, log2f)

  refdata<-refdata_group%>%rename(group_name=1)

  ## get median counts per significant group and clinical variable
  sig_data<-data_all%>%
    group_by(group_name, !!enquo(clin_var))%>%
    summarise(median_counts=median(value), .groups = 'drop')%>%
    pivot_wider(names_from = quo_name(enquo(clin_var)), values_from = "median_counts")%>%
    inner_join(sig_amr, by="group_name")%>%
    inner_join(log_data, by="group_name")%>%
    inner_join(refdata, by="group_name")

  res<-list(data_all, sig_data)
  names(res)<-c("data_plot", "sig_data")
  return(res)
}

##Get significance of a grouped matrix to a clinical variable (factor)

sig_group_clin_factor<-function(data, metadata, clin_var){
  ##transpose data
  group_name<-pull(data, 1)
  data<-as_tibble(cbind(SampleID = names(data), t(data)))%>%slice(-1)%>%
    mutate_at(vars(-("SampleID")),as.numeric)
  colnames(data)<-c("SampleID", group_name)

  ##merge data with metadata
  data_all<-data%>%
    pivot_longer(-SampleID, names_to = "group_name", values_to = "value")%>%
    inner_join(., metadata, by="SampleID")

  ##get significant groups, adjusted p value
  sig_amr<-data_all%>%
    nest(data=-group_name)%>%
    mutate(test=map(.x=data, ~kruskal.test(value~!!ensym(clin_var), data=.x)%>%tidy))%>%
    unnest(test)%>%
    mutate(p.adjust=p.adjust(p.value, method = "BH"))%>%
    filter(p.adjust<0.05)%>%
    select(group_name, p.adjust)

  ##obtain data only from significant groups (for the plot)
  data_all<-data_all%>%
    inner_join(sig_amr, by="group_name")

  ## get mean counts per significant group and clinical variable
  sig_data<-data_all%>%
    group_by(group_name, !!enquo(clin_var))%>%
    summarise(mean_counts=median(value), .groups = 'drop')%>%
    pivot_wider(names_from = quo_name(enquo(clin_var)), values_from = "mean_counts")%>%
    inner_join(sig_amr, by="group_name")

  res<-list(data_all, sig_data)
  names(res)<-c("data_plot", "sig_data")
  return(res)
}


sig_boxplot<-function(data_plot, clin_var, wrap_var){
  clin_var<-enquo(clin_var)
  wrap_var<-enquo(wrap_var)
  data_plot%>%
    rename(wrap_var2=!!wrap_var)%>%
    ggplot(aes(x=!!clin_var, y=value, color=!!clin_var)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
    geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                                jitter.width = 0.5))+
    facet_wrap(~wrap_var2, scales = "free_y")+
    labs(x= NULL, y=NULL) +
    scale_colour_Publication()+
    theme_Publication()+scale_fill_Publication()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank())
}

sig_boxplot_logdata<-function(data_plot, sig_data, clin_var, wrap_var){
  clin_var<-enquo(clin_var)
  wrap_var<-enquo(wrap_var)

  table<-sig_data%>%
    rename(wrap_var2=!!wrap_var)%>%
    mutate(p_adjust=round(p.adjust, 2),
           p_adjust=as.character(p_adjust),
           p_adjust2=ifelse(p_adjust=="0", "<0.005", glue("={p_adjust}")),
           label=glue("p{p_adjust2}\nlog2FC={round(log2f,2)}"))

  data_plot%>%
    rename(wrap_var2=!!wrap_var)%>%
    ggplot(aes(x=!!clin_var, y=value, color=!!clin_var)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
    geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                                jitter.width = 0.5))+
    facet_wrap(~wrap_var2, scales = "free_y")+
    geom_text(data= table,
              mapping = aes(x = Inf, y = Inf, label = label),
              hjust   = 1.1,
              vjust   = 1.1, size= rel(2), inherit.aes = FALSE)+
    labs(x= NULL, y=NULL) +
    scale_colour_Publication()+
    theme_Publication()+scale_fill_Publication()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank())
}

sig_summary_table<-function(genefamily_sigdata, drugclass_sigdata){

  genfam<-genefamily_sigdata%>%
    select(`AMR Gene Family`=group_name, log2f_gene_family=log2f)

  drugclass<-drugclass_sigdata%>%
    select(-`ARO Name`)

  data<-drugclass%>%
    separate_rows(`AMR Gene Family`, sep = ";")%>%
    inner_join(., genfam, by="AMR Gene Family")%>%
    mutate(tendency=ifelse(log2f*log2f_gene_family>0, "same tendency", "different tendency"))%>%
    select(group_name, `AMR Gene Family`, tendency)%>%
    group_by(group_name, tendency)%>%
    summarise(gene_family = paste(rle(`AMR Gene Family`)$values, collapse=";"))%>%
    pivot_wider(names_from = "tendency", values_from = "gene_family")

  data_sum<-drugclass%>%
    select(-`AMR Gene Family`)%>%
    left_join(., data, by="group_name")%>%
    arrange(log2f)
}



