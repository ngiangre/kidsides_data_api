#
# This is a Plumber API. You can run the API by clicking
# the 'Run API' button above.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#


# Setup -------------------------------------------------------------------

library(plumber)
library(tidyverse)
library(RSQLite)
library(data.table)
library(pool)
library(cowplot)


# db connection -----------------------------------------------------------


con <- dbPool(
    drv = RSQLite::SQLite(),
    dbname = "effect_peds_19q2_v0.3_20211119.sqlite"
)

# data and plot setup -----------------------------------------------------


tmp <- 
    tbl(con,"ade_raw") %>% 
    select(safetyreportid,nichd,sex) %>% 
    distinct()

total_stage_reports <- 
    tmp %>% 
    group_by(nichd) %>% 
    count(name = "total") %>% 
    collect() %>% 
    data.table()
total_stage_sex_reports <- 
    tmp %>% 
    group_by(nichd,sex) %>% 
    count(name = "total") %>% 
    collect() %>% 
    data.table()

min_date <- 
    tbl(con,"ade_raw") %>% 
    select(receive_date) %>% 
    arrange(receive_date) %>% 
    head(1) %>% 
    collect() %>% unlist %>% unname
max_date <- 
    tbl(con,"ade_raw") %>% 
    select(receive_date) %>% 
    arrange(desc(receive_date)) %>% 
    head(1) %>% 
    collect() %>% unlist %>% unname

drug_table <-
    tbl(con,"drug") %>%
    collect() %>%
    data.table() %>%
    na.omit() %>%
    .[,.(atc_concept_id,N = ndrugreports,
         code = atc_concept_code,
         ATC5=atc_concept_name,
         ATC4=atc4_concept_name,
         ATC3=atc3_concept_name,
         ATC2=atc2_concept_name,
         ATC1=atc1_concept_name)] %>% 
    .[order(N,decreasing = T)]

drugNames <- 
    drug_table[order(N,decreasing = T),.(atc_concept_id,code,N,ATC5)][,paste0(ATC5," [",code,"] (N=",scales::comma(N,accuracy = 1),")")]
drugIDs <- 
    drug_table[order(N,decreasing = T),.(atc_concept_id,ATC5)][,atc_concept_id]
names(drugIDs) <- drugNames

event_table <-
    tbl(con,"event") %>%
    collect() %>%
    data.table() %>%
    na.omit() %>%
    .[,.(meddra_concept_id,N = neventreports,
         code = meddra_concept_code_1,
         PT=meddra_concept_name_1,
         HLT=meddra_concept_name_2,
         HLGT=meddra_concept_name_3,
         SOC=meddra_concept_name_4)]

null_dist <- 
    tbl(con,"ade_null_distribution") %>% 
    collect() %>% 
    data.table()

null_dist_summary <- 
    null_dist %>% 
    .[,
      .(
          null_lwr = quantile(gam_score,c(0.01)),
          null_mean = mean(gam_score),
          null_upr = quantile(gam_score,c(0.99))
      ),
      nichd
    ]

theme_big <- 
    theme_classic(base_size=16) + 
    theme(
        strip.text = element_text(color="black",face="bold"),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(color="black",face="bold",
                                   angle=45,vjust=1,hjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color="black",face="bold"),
        legend.position = "top",
        legend.text = element_text(color="black",face="bold"),
        legend.key.size = unit(0.4,"cm"),
        legend.box.margin = margin(-0.4,-0.4,0,-0.4,unit = "cm")
    )

theme_small <- 
    theme_classic(base_size=12) + 
    theme(
        strip.text = element_text(color="black",face="bold",size=8),
        axis.title.x = element_text(color="black",face="bold",size=8),
        axis.title.y = element_text(color="black",face="bold",size=8),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1,color="black",face="bold",size=6),
        axis.text.y = element_text(color="black",face="bold",size=8),
        legend.text = element_text(color="black",face="bold",size=8),
        legend.key.size = unit(0.4,"cm"),
        legend.box.margin = margin(-0.4,-0.4,-0.4,-0.4,unit = "cm")
    )

stages <-
    c("term_neonatal","infancy",
      "toddler","early_childhood",
      "middle_childhood","early_adolescence",
      "late_adolescence")
stages_split <- str_replace(stages,"_","\n")

ade_cohort <- function(drugs=c(),events=c(),database=con){
    
    if(length(events)==0 | length(drugs)==0){errorCondition("no drugs or events given")}
    
    tmp <-
        expand.grid(
            atc_concept_id = drugs,
            meddra_concept_id = events
        ) %>%
        merge(
            tbl(con,"ade_nichd") %>%
                filter(
                    atc_concept_id %in% drugs &
                        meddra_concept_id %in% events) %>%
                collect() %>% data.table(),
            by=c("atc_concept_id","meddra_concept_id")
        ) %>%
        data.table() %>%
        merge(
            tbl(con,"ade") %>%
                filter(
                    atc_concept_id %in% drugs &
                        meddra_concept_id %in% events) %>%
                collect() %>% data.table(),
            by=c("ade","atc_concept_id","meddra_concept_id")
        ) %>%
        merge(
            tbl(con,"ade_null") %>% 
                collect() %>% 
                data.table(),
            by="nichd"
        ) %>% 
        merge(
            tbl(con,"event") %>%
                filter(
                    meddra_concept_id %in% events
                ) %>%
                select(meddra_concept_id,meddra_concept_name_1) %>%
                rename(meddra_concept_name = meddra_concept_name_1) %>%
                collect() %>% data.table() %>% na.omit() %>% unique(),
            by=c("meddra_concept_id")
        ) %>%
        merge(
            tbl(con,"drug") %>%
                filter(
                    atc_concept_id %in% drugs
                ) %>%
                select(atc_concept_id,atc_concept_name) %>%
                collect() %>% data.table() %>% na.omit() %>% unique(),
            by=c("atc_concept_id")
        )
    
    tmp$significance <-
        ifelse(tmp$gt_null_statistic==1,"Nominal","NA")
    tmp[gt_null_99==1,"significance"]="Null model"
    
    tmp
    
}

plot_ade_risks_null_shade <- function(x,color="red",theme=theme_big,ts=5){
    x$nichd_split <- str_replace(x$nichd,"_","\n")
    x$NICHD = factor(x$nichd_split,levels=str_replace(stages,"_","\n"))
    x %>% 
        merge(
            null_dist_summary,
            by="nichd"
        ) %>%
        ggplot(aes(factor(nichd_split,levels=stages_split),gam_score,group=ade)) +
        geom_ribbon(aes(ymin=null_lwr,ymax=null_upr),
                    fill="lightgray",alpha=0.3) +
        geom_bar(stat="identity",aes(y=log10(as.numeric(DE)+1)),fill="lightgray",color="black") +
        ggrepel::geom_label_repel(
            aes(
                label=DE,
                y=log10(as.numeric(DE)+1)
            ),
            vjust=1.3,fontface="bold",size=ts) +
        geom_point(color=color,size=2) +
        geom_errorbar(aes(ymin=gam_score_90mse,ymax=gam_score_90pse),
                      width=0.1,color=color,size=1) +
        geom_hline(yintercept = 0,color="black",size=0.5) +
        facet_wrap(~ade_name,labeller = label_wrap_gen(width = 20)) +
        xlab("") +
        ylab("Risk of ADE (GAM log odds)") + 
        theme +
        theme(
            axis.line = element_line(),
            axis.ticks = element_line(),
            axis.text.y = element_text(color="black",face="bold",size=12)
        )
    
}

population_summary_stage <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
    pop <- 
        tbl(con,"ade_raw") %>% 
        filter(
            atc_concept_id %in% drugs &
                meddra_concept_id %in% events
        ) %>% 
        collect() %>% 
        data.table()
    pop$nichd_split <- str_replace(pop$nichd,"_","\n")
    
    tmp <- 
        merge(
            total_stage_reports,
            pop[,
                .(safetyreportid,nichd)
            ] %>% 
                unique() %>% 
                .[,.N,.(nichd)],
            all.x=T
        ) %>%
        .[,.(nichd,nichd_split = str_replace(nichd,"_","\n"),N = ifelse(is.na(N),0,N))]
    
    merge(
        tmp,
        total_stage_reports,
        by="nichd",
        all.y = T
    ) %>% 
        .[,.(nichd_split,N,total,
             prop = (N/total),
             prop_norm = (((N/total) - min((N/total)))/(max((N/total)) - min((N/total)))),
             prop_norm_mm = (((N/total) - min((N/total)))/(max((N/total)) - min((N/total))))*max(N)
        )
        ] %>% 
        ggplot(aes(factor(nichd_split,levels=stages_split),N)) +
        geom_bar(stat="identity",color="black",fill='lightgray') +
        geom_point(aes(y=prop_norm_mm),color="yellowgreen",pch=21) +
        geom_line(aes(y=prop_norm_mm,group=1),color="yellowgreen",linetype="dotted") +
        ggrepel::geom_label_repel(aes(label=paste0(round(prop*100,2),"% of ",scales::comma(total)),y=prop_norm_mm),
                                  color="yellowgreen",size=ts,fontface="bold",
                                  fill = "white") +
        xlab("") +
        coord_cartesian(clip="off") +
        scale_y_continuous(labels=scales::comma) +
        geom_text(aes(label=N,y=N),vjust=-1,fontface="bold",size=ts) +
        ylab("Number of reports") +
        theme(
            legend.position = "right",
            axis.text.x = element_text(angle=45,vjust=1,hjust=1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm")
        ) +
        theme
}

population_summary_sex <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
    pop <- 
        tbl(con,"ade_raw") %>% 
        filter(
            atc_concept_id %in% drugs &
                meddra_concept_id %in% events
        ) %>% 
        collect() %>% 
        data.table()
    pop$nichd_split <- str_replace(pop$nichd,"_","\n")
    
    merge(
        total_stage_sex_reports,
        pop[,
            .(safetyreportid,nichd,sex)
        ] %>% 
            unique() %>% 
            .[,.N,.(sex,nichd)],
        all.x=T
    ) %>%
        .[,.(sex,nichd,nichd_split = str_replace(nichd,"_","\n"),N = ifelse(is.na(N),0,N))] %>%
        ggplot(aes(factor(nichd_split,levels=stages_split),N,fill=sex)) +
        geom_bar(stat="identity",position = position_dodge(width=0.9),color="black") +
        xlab("") +
        coord_cartesian(clip="off") +
        scale_y_continuous(labels=scales::percent) +
        geom_text(aes(label=N,y=N),vjust=-1,position = position_dodge(width=0.9),fontface="bold",size=ts) +
        colorspace::scale_fill_discrete_qualitative(breaks=c("Female","Male"),labels=c("Female","Male")) +
        guides(fill=guide_legend(title = NULL)) +
        ylab("Number of reports") +
        theme
}

population_summary_poly <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
    pop <- 
        tbl(con,"ade_raw") %>% 
        filter(
            atc_concept_id %in% drugs &
                meddra_concept_id %in% events
        ) %>% 
        collect() %>% 
        data.table()
    pop$nichd_split <- str_replace(pop$nichd,"_","\n")
    
    pop[,
        .(safetyreportid,polypharmacy,nichd_split)
    ] %>% 
        merge(
            data.table(nichd_split = stages_split),
            all.y=T
        ) %>%
        ggplot(aes(factor(nichd_split,levels=stages_split),polypharmacy)) +
        ggbeeswarm::geom_quasirandom(groupOnX = T) +
        xlab("") +
        ylab("Number of drugs") +
        theme +
        theme(
            axis.text.y = element_text(color="black",face="bold"),
            axis.line = element_line()
        )
}

population_summary_reporter <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
    pop <- 
        tbl(con,"ade_raw") %>% 
        filter(
            atc_concept_id %in% drugs &
                meddra_concept_id %in% events
        ) %>% 
        collect() %>% 
        data.table()
    pop$nichd_split <- str_replace(pop$nichd,"_","\n")
    pop$reporter_qualification <- str_replace(pop$reporter_qualification," ","\n")
    
    merge(
        expand.grid(
            nichd_split = stages_split,
            reporter_qualification = pop[,unique(reporter_qualification)]
        ) %>% data.table(),
        pop[
            ,.(safetyreportid,reporter_qualification,nichd_split)
        ] %>% 
            unique() %>% 
            .[,
              .(N = .N),
              .(nichd_split,reporter_qualification)
            ],
        all.x=T
    ) %>%
        .[,.(nichd_split,reporter_qualification,N = ifelse(is.na(N),0,N))] %>% 
        ggplot(aes(factor(nichd_split,levels=stages_split),N,fill=reporter_qualification)) +
        geom_bar(stat="identity",position = position_dodge(width=0.9),color="black") +
        xlab("") +
        scale_y_continuous(labels=scales::percent) +
        coord_cartesian(clip="off") +
        geom_text(aes(label=N,y=N),vjust=-1,position = position_dodge(width=0.9),fontface="bold",size=ts) +
        colorspace::scale_fill_discrete_qualitative() +
        guides(fill=guide_legend(title=NULL,ncol=2)) +
        ylab("Number of reports") +
        theme
}

population_summary_atc <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
    pop <- 
        tbl(con,"ade_raw") %>% 
        filter(
            atc_concept_id %in% drugs &
                meddra_concept_id %in% events
        ) %>% 
        collect() %>% 
        data.table()
    pop$nichd_split <- str_replace(pop$nichd,"_","\n")
    reports=pop[,unique(safetyreportid)]
    
    atc_names <- 
        tbl(con,"atc_raw_map") %>% 
        collect() %>% 
        data.table() %>% 
        .[,unique(atc1_concept_name[order(atc1_concept_name)])] %>% 
        str_replace_all(" ","\n")
    tmp <- 
        merge(
            expand.grid(
                nichd_split = stages_split,
                atc1_concept_name = atc_names
            ) %>% data.table() %>% 
                merge(
                    tbl(con,"drug") %>% 
                        select(atc1_concept_name,atc1_concept_code) %>% 
                        distinct() %>% 
                        collect() %>% 
                        data.table() %>% 
                        .[,.(code = atc1_concept_code,
                             atc1_concept_name = str_replace_all(atc1_concept_name," ","\n"))],
                    by="atc1_concept_name"
                ),
            merge(
                tbl(con,"ade_raw") %>% 
                    filter(safetyreportid %in% reports &
                               !(atc_concept_id %in% drugs)) %>% 
                    select(safetyreportid,atc_concept_id,nichd) %>% 
                    collect() %>% 
                    data.table() %>% 
                    .[,.(safetyreportid,atc_concept_id,nichd_split = str_replace_all(nichd,"_","\n"))] %>% 
                    unique(),
                tbl(con,"drug") %>% 
                    collect() %>% 
                    data.table() %>% 
                    .[,.(atc_concept_id,code = atc1_concept_code,
                         atc1_concept_name = str_replace_all(atc1_concept_name," ","\n"))],
                by="atc_concept_id"
            ) %>% 
                .[,.(safetyreportid,code,atc1_concept_name,code,nichd_split)] %>% 
                unique() %>% 
                .[,.N,
                  .(atc1_concept_name,
                    nichd_split,code)
                ] %>% na.omit(),
            all.x=T,
            by=c("atc1_concept_name","nichd_split","code")
        )
    tmp %>% 
        .[atc1_concept_name!="VARIOUS",
          .(nichd_split,
            atc1_concept_name = paste0(atc1_concept_name,"\n(ATC 1st code: ",code,")"),
            N = ifelse(is.na(N),0,N))] %>% 
        ggplot(aes(factor(nichd_split,levels=stages_split),N)) +
        geom_bar(stat="identity",color="black",fill='lightgray') +
        xlab("") +
        coord_cartesian(clip="off") +
        scale_y_continuous(labels=scales::comma) +
        geom_text(aes(label=N,y=N),vjust=-1,fontface="bold",size=ts) +
        facet_wrap(~atc1_concept_name~.,labeller = label_wrap_gen(width=20)) +
        ylab("Number of reports") +
        theme +
        theme(
            strip.text = element_text(face="bold",size=12),
            axis.text.x = element_text(face="bold",size=10)
        )
}

population_summary_date <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
    pop <- 
        tbl(con,"ade_raw") %>% 
        filter(
            atc_concept_id %in% drugs &
                meddra_concept_id %in% events
        ) %>% 
        collect() %>% 
        data.table()
    pop$nichd_split <- str_replace(pop$nichd,"_","\n")
    pop$receive_date <- as.Date(pop$receive_date)
    
    range_ <- seq.Date(as.Date(min_date),as.Date(max_date),by = 1)
    date_range <- data.table(receive_date = range_)
    date_range$interval <- cut(range_,breaks = "quarters")
    
    pop[,
        .(receive_date,nichd_split,safetyreportid)
    ] %>% 
        unique() %>% 
        merge(
            expand.grid(
                nichd_split = stages_split,
                receive_date = range_
            ) %>% data.table(),
            all.y=T
        ) %>%
        merge(
            date_range,
            by="receive_date",
            all.x=T
        ) %>% 
        .[,.(N = length(unique(na.omit(safetyreportid)))),.(interval,nichd_split)] %>% 
        .[,.(interval = as.Date(interval),N,nichd_split)] %>% 
        ggplot(aes(interval,N,group=1)) +
        geom_point() +
        geom_line() +
        scale_x_date(breaks="3 years",date_labels = "%Y") +
        scale_y_continuous(labels=scales::number_format(accuracy=1),lim=c(0,NA)) +
        facet_wrap(~factor(nichd_split,stages_split)) +
        xlab("Time") +
        ylab("Number\nof reports") +
        theme +
        theme(
            axis.text.y = element_text(color="black",face="bold"),
            axis.line = element_line()
        )
}

population_summary <- function(drugs=c(),events=c(),database=con,theme=theme_small){
    
    pop <- 
        tbl(con,"ade_raw") %>% 
        filter(
            atc_concept_id %in% drugs &
                meddra_concept_id %in% events
        ) %>% 
        collect() %>% 
        data.table()
    pop$receive_date <- as.Date(pop$receive_date)
    pop$nichd_split <- str_replace(pop$nichd,"_","\n")
    pop$reporter_qualification <- str_replace(pop$reporter_qualification," ","\n")
    reports=pop[,unique(safetyreportid)]
    bm <- 0.2
    ts <- 2.5
    ##stage
    tmp <- 
        merge(
            total_stage_reports,
            pop[,
                .(safetyreportid,nichd)
            ] %>% 
                unique() %>% 
                .[,.N,.(nichd)],
            all.x=T
        ) %>%
        .[,.(nichd,nichd_split = str_replace(nichd,"_","\n"),N = ifelse(is.na(N),0,N))]
    
    gstage <- 
        merge(
            tmp,
            total_stage_reports,
            by="nichd",
            all.y = T
        ) %>% 
        .[,.(nichd_split,N,total,
             prop = (N/total),
             prop_norm = (((N/total) - min((N/total)))/(max((N/total)) - min((N/total)))),
             prop_norm_mm = (((N/total) - min((N/total)))/(max((N/total)) - min((N/total))))*max(N)
        )
        ] %>% 
        ggplot(aes(factor(nichd_split,levels=stages_split),N)) +
        geom_bar(stat="identity",color="black",fill='lightgray') +
        geom_point(aes(y=prop_norm_mm),color="yellowgreen",pch=21) +
        geom_line(aes(y=prop_norm_mm,group=1),color="yellowgreen",linetype="dashed") +
        geom_line(aes(y=prop_norm_mm,group=1),color="yellowgreen",linetype="dotted") +
        ggrepel::geom_label_repel(aes(label=paste0(round(prop*100,2),"% of ",scales::comma(total)),y=prop_norm_mm),
                                  color="yellowgreen",size=ts,fontface="bold",
                                  fill = "white") +
        xlab("") +
        coord_cartesian(clip="off") +
        scale_y_continuous(labels=scales::comma) +
        geom_text(aes(label=N,y=N),vjust=-1,fontface="bold",size=ts) +
        ylab("Number\nof reports") +
        theme +
        theme(
            legend.position = "right",
            axis.text.x = element_text(angle=45,vjust=1,hjust=1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_blank(),
            plot.margin = unit(c(0.4,0,0,0), "cm")
        )
    ##Sex
    gsex <- 
        merge(
            total_stage_sex_reports,
            pop[,
                .(safetyreportid,nichd,sex)
            ] %>% 
                unique() %>% 
                .[,.N,.(sex,nichd)],
            all.x=T
        ) %>%
        .[,.(sex,nichd,nichd_split = str_replace(nichd,"_","\n"),N = ifelse(is.na(N),0,N))] %>%
        ggplot(aes(factor(nichd_split,levels=stages_split),N,fill=sex)) +
        geom_bar(stat="identity",position = position_dodge(width=0.9),color="black") +
        xlab("") +
        coord_cartesian(clip="off") +
        scale_y_continuous(labels=scales::percent) +
        geom_text(aes(label=N,y=N),vjust=-1,position = position_dodge(width=0.9),fontface="bold",size=ts) +
        colorspace::scale_fill_discrete_qualitative(breaks=c("Female","Male"),labels=c("Female","Male")) +
        guides(fill=guide_legend(title = NULL)) +
        ylab("Number\nof reports") +
        theme +
        theme(
            legend.position = "top",
            axis.text.x = element_text(angle=45,vjust=1,hjust=1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm"),
            legend.box.margin = margin(-0.4,0.6,0,-0.4,unit = "cm")
        ) 
    ##Polypharmacy
    gpoly <- 
        pop[,
            .(safetyreportid,polypharmacy,nichd_split)
        ] %>% 
        merge(
            data.table(nichd_split = stages_split),
            all.y=T
        ) %>%
        mutate(nichd_split = factor(nichd_split,levels=stages_split)) |> 
        na.omit() |> 
        ggplot(aes(nichd_split,polypharmacy)) +
        scale_x_discrete(drop=FALSE) +
        ggbeeswarm::geom_quasirandom(groupOnX = T,size=0.5) +
        xlab("") +
        ylab("Polypharmacy") +
        theme +
        theme(
            axis.text.x = element_text(angle=45,vjust=1,hjust=1),
            axis.line = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm")
        ) 
    
    ##Reporter
    grep <- merge(
        expand.grid(
            nichd_split = stages_split,
            reporter_qualification = pop[,unique(reporter_qualification)]
        ) %>% data.table(),
        pop[
            ,.(safetyreportid,reporter_qualification,nichd_split)
        ] %>% 
            unique() %>% 
            .[,
              .(N = .N),
              .(nichd_split,reporter_qualification)
            ],
        all.x=T
    ) %>%
        .[,.(nichd_split,reporter_qualification,N = ifelse(is.na(N),0,N))] %>% 
        ggplot(aes(factor(nichd_split,levels=stages_split),N,fill=reporter_qualification)) +
        geom_bar(stat="identity",position = position_dodge(width=0.9),color="black") +
        xlab("") +
        scale_y_continuous(labels=scales::percent) +
        coord_cartesian(clip="off") +
        geom_text(aes(label=N,y=N),vjust=-1,position = position_dodge(width=0.9),fontface="bold",size=ts) +
        colorspace::scale_fill_discrete_qualitative() +
        guides(fill=guide_legend(title=NULL,ncol=2)) +
        ylab("Number\nof reports") +
        theme +
        theme(
            legend.position = "top",
            strip.background = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_blank(),
            axis.text.x = element_text(angle=45,vjust=1,hjust=1),
            plot.margin = unit(c(0,0,0,0), "cm"),
            legend.box.margin = margin(-0.4,0.6,0,-0.4,unit = "cm")
        )
    
    ##ATC class proportion
    atc_names <- 
        tbl(con,"atc_raw_map") %>% 
        collect() %>% 
        data.table() %>% 
        .[,unique(atc1_concept_name[order(atc1_concept_name)])] %>% 
        str_replace_all(" ","\n")
    tmp <- 
        merge(
            expand.grid(
                nichd_split = stages_split,
                atc1_concept_name = atc_names
            ) %>% data.table() %>% 
                merge(
                    tbl(con,"drug") %>% 
                        select(atc1_concept_name,atc1_concept_code) %>% 
                        distinct() %>% 
                        collect() %>% 
                        data.table() %>% 
                        .[,.(code = atc1_concept_code,
                             atc1_concept_name = str_replace_all(atc1_concept_name," ","\n"))],
                    by="atc1_concept_name"
                ),
            merge(
                tbl(con,"ade_raw") %>% 
                    filter(safetyreportid %in% reports &
                               !(atc_concept_id %in% drugs)) %>% 
                    select(safetyreportid,atc_concept_id,nichd) %>% 
                    collect() %>% 
                    data.table() %>% 
                    .[,.(safetyreportid,atc_concept_id,nichd_split = str_replace_all(nichd,"_","\n"))] %>% 
                    unique(),
                tbl(con,"drug") %>% 
                    collect() %>% 
                    data.table() %>% 
                    .[,.(atc_concept_id,code = atc1_concept_code,
                         atc1_concept_name = str_replace_all(atc1_concept_name," ","\n"))],
                by="atc_concept_id"
            ) %>% 
                .[,.(safetyreportid,code,atc1_concept_name,code,nichd_split)] %>% 
                unique() %>% 
                .[,.N,
                  .(atc1_concept_name,
                    nichd_split,code)
                ] %>% na.omit(),
            all.x=T,
            by=c("atc1_concept_name","nichd_split","code")
        )
    gatc <- tmp %>% 
        .[atc1_concept_name!="VARIOUS",
          .(nichd_split,
            atc1_concept_name,
            N = ifelse(is.na(N),0,N))] %>% 
        ggplot(aes(factor(nichd_split,levels=stages_split),N)) +
        geom_bar(stat="identity",color="black",fill='lightgray') +
        xlab("") +
        coord_cartesian(clip="off") +
        scale_y_continuous(labels=scales::comma) +
        geom_text(aes(label=N,y=N),vjust=-1,fontface="bold",size=ts) +
        facet_wrap(~atc1_concept_name~.,labeller = label_wrap_gen(width=20)) +
        ylab("Number of reports") +
        theme +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(face="bold",size=12),
            axis.text.x = element_text(face="bold",size=10)
        )
    
    #Date
    range_ <- seq.Date(as.Date(min_date),as.Date(max_date),by = 1)
    date_range <- data.table(receive_date = range_)
    date_range$interval <- cut(range_,breaks = "quarters")
    
    gdate <- 
        pop[,
            .(receive_date,nichd_split,safetyreportid)
        ] %>% 
        unique() %>% 
        merge(
            expand.grid(
                nichd_split = stages_split,
                receive_date = range_
            ) %>% data.table(),
            all.y=T
        ) %>%
        merge(
            date_range,
            by="receive_date",
            all.x=T
        ) %>% 
        .[,.(N = length(unique(na.omit(safetyreportid)))),.(interval,nichd_split)] %>% 
        .[,.(interval = as.Date(interval),N,nichd_split)] %>% 
        ggplot(aes(interval,N,group=1)) +
        geom_point() +
        geom_line() +
        scale_x_date(breaks="3 years",date_labels = "%Y") +
        scale_y_continuous(labels=scales::number_format(accuracy=1),lim=c(0,NA)) +
        facet_wrap(~factor(nichd_split,stages_split)) +
        xlab("Time") +
        ylab("Number\nof reports") +
        theme +
        theme(
            strip.background = element_blank(),
            axis.text.x = element_text(angle=45,vjust=1,hjust=1),
            axis.line = element_blank(),
            legend.position = "none"
        )
    
    #https://wilkelab.org/cowplot/articles/plot_grid.html
    plot_grid(
        plot_grid(
            plot_grid(gstage,gsex,gpoly,ncol=1,nrow=3),
            plot_grid(gatc),
            nrow=1,ncol=2,
            rel_widths = c(2,4.5)
        ),
        plot_grid(grep,gdate,nrow=1,ncol=2,rel_widths = c(2.5,2)),
        nrow=2,ncol=1,
        rel_heights = c(4.5,2)
    )
    
}

get_drug_name <- function(x){
    drug_table %>%
        .[atc_concept_id %in% x,
          unique(ATC5)]
}

get_event_name <- function(x){
    event_table %>%
        .[meddra_concept_id %in% x,
          unique(PT)]
}


# Plumber API -------------------------------------------------------------


#* @apiTitle Plumber API for KidSIDES and the PDSPortal
#* @apiDescription Plumber example description.

#* Return table from KidSIDES database
#* @param table The table to retrieve from the database
#* @get /table_json
function(table) {
    match.arg(table,DBI::dbListTables(conn))
    dplyr::tbl(conn,table) |> 
        dplyr::collect()
}

#* Return table from KidSIDES database
#* @param table The table to retrieve from the database
#* @get /table_csv
#* @serializer csv
function(table) {
    match.arg(table,DBI::dbListTables(conn))
    dplyr::tbl(conn,table) |> 
        dplyr::collect()
}

#* Return named list of drug names and identifiers
#* @get /drug_names_json
function() {
    drugIDs |> 
        tibble::enframe()
}

#* Return named list of drug names and identifiers
#* @get /drug_names_csv
#* @serializer csv
function() {
    drugIDs |> 
        tibble::enframe()
}

#* Return named list of event names and identifiers co-reported with drug(s)
#* @param drugs_ Identifiers for drug(s)
#* @get /event_names_json
function(drugs_) {
    event_count <-
        tbl(con,"ade_raw") %>%
        filter(atc_concept_id %in% drugs_) %>%
        distinct(safetyreportid,meddra_concept_id) %>%
        group_by(meddra_concept_id) %>%
        count() %>%
        collect() %>%
        data.table()
    
    event_table_update <-
        merge(
            event_table[,.(meddra_concept_id,code,PT,SOC)],
            event_count[,.(meddra_concept_id,N=n)],
            by="meddra_concept_id"
        )
    
    eventNames <-
        event_table_update[
            order(N,decreasing = T),
            .(meddra_concept_id,code,N,PT,SOC)
        ][,
          paste0(stringr::str_replace_all(PT,' ','_'),
                 " [",
                 stringr::str_replace_all(SOC,' ','_'),
                 "] (N=",
                 scales::comma(N,accuracy = 1),
                 ")")
        ]
    eventIDs <-
        event_table_update[order(N,decreasing = T),.(meddra_concept_id,code,N,PT)][,meddra_concept_id]
    names(eventIDs) <- eventNames
    
    eventIDs |> 
        tibble::enframe()
    
}

#* Return an infographic of the population summary for the co occurrence of drug(s) and event(s).
#* @param drug_vec Vector of drugs
#* @param event_vec Vector of events
#* @get /population_summary
#* @serializer png list(width = 800, height = 800)
function(drug_vec, event_vec) {
    population_summary(
        drug_vec,
        event_vec
    ) |> print()
}
