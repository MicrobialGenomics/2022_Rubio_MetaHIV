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
    theme_classic()+
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
    theme_bw()+
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
    theme_bw()+
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
    scale_color_manual(name=quo_name(clin_var),
                       values = c("blue", "red","green4"))+
    scale_fill_manual(values = c("dodgerblue", "pink", "green"))+
    theme_bw()+
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
    labs(caption = glue("ADONIS p={ptest}, r2={r2test}"))+
    scale_color_manual(name=quo_name(clin_var),
                       values = c("blue", "red","green4", "orange", "aquamarine4", "magenta2", "gold", "black"))+
    theme_bw()+
    theme(legend.text = element_text(size=10))}

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

  data_all<-data_all%>%
    inner_join(sig_amr, by="ref_name")

  sig_amr<-sig_amr%>%
    inner_join(refdata, by="ref_name")

  res<-list(data_all, sig_amr)
  names(res)<-c("data_plot", "sig_amr")
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

sig_group_clin_dic<-function(data, metadata, clin_var){
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


  ## get median counts per significant group and clinical variable
  sig_data<-data_all%>%
    group_by(group_name, !!enquo(clin_var))%>%
    summarise(median_counts=median(value), .groups = 'drop')%>%
    pivot_wider(names_from = quo_name(enquo(clin_var)), values_from = "median_counts")%>%
    inner_join(sig_amr, by="group_name")%>%
    inner_join(log_data, by="group_name")

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


