load("C:/Users/kwebb/Documents/GitHub/dp-for-rct-paper-replication/output/liberia/liberia_ITT_compare_methods_data_v0224.Rda")
#`covset`, `treatment`, `response`, `method`, `ci_overlap`, `absdiff`
names(data.list)
table(full.comparison$data,full.comparison$Sim)

save_compare<-full.comparison|>rename(covset=Sim,response=ModelName)|>
  select(covset,response,data,
         paste0("ci.overlap_",c("cashassonly","tpassonly","tpcashass","Control")),
         paste0("ITT_",c("cashassonly","tpassonly","tpcashass","Control")))
save_compare$method=factor(save_compare$data,levels=c("Confidential","DP-Mb Same Privacy",
                                                      "DP-Mb_Eps2_Del0","DP-Mb_Eps9_Del0","Hybrid","MV Histogram"),
                           labels=c("Confidential","GenM-1",
                                    "GenM-2","GenM-9","Hybrid","MV Histogram"))
save_compare_pivs_itts<-save_compare%>%select(-data)%>%
  tidyr::pivot_longer(cols=-all_of(c("covset","response","method")),
                      names_to=c("measure","treatment"),names_sep="_", values_to="vals")%>%
  tidyr::pivot_wider(names_from=measure,values_from=vals)


confitt=save_compare_pivs_itts%>%filter(method=="Confidential")%>%
  select(covset,response, treatment, ITT)%>%rename(conf_ITT=ITT)

merge_df=save_compare_pivs_itts%>%filter(method!="Confidential")%>%
  left_join(confitt, by=c("covset","response", "treatment"))%>%
  mutate(absdiff=abs(ITT-conf_ITT))%>%
  select(covset, method, absdiff, ci.overlap, response, treatment)%>%
  rename(ci_overlap=ci.overlap)



 head(merge_df)

 write.csv(merge_df,file="output/liberia_utility.csv")
