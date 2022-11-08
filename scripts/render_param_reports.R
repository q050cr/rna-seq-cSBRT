


rm(list=ls())  # not in rmd doc, otherwise params deleted

# define render function to stratify groups
render_report = function(strat_groups, pval_threshold, continous, LV_only) {
  rmarkdown::render(
    "rna-seq-rahm/scripts/deseq_radiation_parameterized.Rmd", 
    params = list(
      strat_groups = strat_groups, 
      LV_only = LV_only, 
      continous = continous, 
      pval_threshold = pval_threshold, # default 0.1 in DESeq2
      doc_title = paste0("cSBRT - Differential Gene Expression Analysis - ", strat_groups)
    ),
    output_file = paste0(
      format(Sys.time(), "%Y%m%d"), '_',  strat_groups, "_cont_analysis_", continous, 
      "_LV_only_", LV_only, "_pval_", pval_threshold, "_cSBRT_analysis", '.html'
      ),
    output_dir = "rna-seq-rahm/reports-deseq",
    envir = globalenv()
  )
}


# also RV analysis included..
render_report(strat_groups = "radiation3groups", LV_only = FALSE, continous = FALSE, pval_threshold = 0.05)
render_report(strat_groups = "radiation2groups", LV_only = FALSE, continous = FALSE, pval_threshold = 0.1) # no DEgenes with 0.05

## LV only!
render_report(strat_groups = "radiation2groups", LV_only = TRUE, continous = FALSE, pval_threshold = 0.05)  # no DEgenes
render_report(strat_groups = "radiation2groups", LV_only = TRUE, continous = FALSE, pval_threshold = 0.1)
render_report(strat_groups = "radiation3groups", LV_only = TRUE, continous = FALSE, pval_threshold = 0.05)

## analysis continous radiation
render_report(strat_groups = "radiation3groups", LV_only = FALSE, continous = TRUE, pval_threshold = 0.05) 
render_report(strat_groups = "radiation3groups", LV_only = TRUE, continous = TRUE, pval_threshold = 0.05) 

render_report(strat_groups = "radiation3groups", LV_only = TRUE, continous = TRUE, pval_threshold = 0.1)


