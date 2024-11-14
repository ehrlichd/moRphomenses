
input <- list(
  "tar_length" = 28,
  "tar_alignment" = 16,
  "knn_n" = 3,
  "ind_pick1" = 1
)

rv <- list(
  "mm_dataset" = mm_data,
  "assign_id" = "ID",
  "assign_x" = "CYCLEDAY",
  "assign_y" = "E1G",
  "assign_aln" = "MIDPOINT"
)

aln_11 <-  mm_ArrayData(IDs = rv$mm_dataset %>%
                          pull(!!sym(eval(rv$assign_id))),
                        DAYS = rv$mm_dataset %>%
                          pull(!!sym(eval(rv$assign_x))),
                        VALUE = rv$mm_dataset %>%
                          pull(!!sym(eval(rv$assign_y))),
                        MID = rv$mm_dataset %>%
                          pull(!!sym(eval(rv$assign_aln))),
                        targetLENGTH = input$tar_length,
                        targetMID = input$tar_alignment,
                        transformation = "minmax",
                        impute_missing = input$knn_n

)
