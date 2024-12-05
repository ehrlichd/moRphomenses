




library(shiny)
library(dplyr)
library(moRphomenses)

pc_shape_height <- 100
ind_shape_height <- 156

pc_pad <- 50



# Application title ------------------------

fluidPage(
titlePanel("mm Explorer"),

# tabset layout ---------------------------
tabsetPanel(
# # tab Load Data --------------------------
# tabPanel("Load Data", sidebarLayout(
#   sidebarPanel(width = 2, fileInput("file_select", "Choose data file")),
#   ## close SB
#   mainPanel(
#     actionButton("go", "Analyze!"),
#     uiOutput("dyn_var_assign")
#     ) ## close MP
# )), ## close sbpLayout; close tab


tabPanel("Introduction",

         ),

# tab alignment -------------------------
tabPanel("Alignment", sidebarLayout(
  sidebarPanel(
    width = 2,

    ## input$tar_alignment -----------------
    numericInput(
      "tar_alignment",
      label = "Align by day:",
      value = 16,
      min = 1,
      max = 45
    ),

    ## input$tar_length --------------------
    numericInput(
      "tar_length",
      "Resample cycles to l-days",
      value = 28,
      min = 10,
      max = 50
    ),

    ## input$knn_n ---------------------------------
    numericInput(
      "knn_n",
      "Impute missing using n-neighbors:",
      value = 3,
      min = 1,
      max = 10
    ),

    uiOutput("dyn_pick1")

    ## input$scale_type ---------------------
    #### implement in next update

  ),

  ## mainpanel ---------------------------
  mainPanel(
    fluidRow(
      column(
        width = 6,
        ## plot of all individuals
        ## hover to highlight a single line
        ## click to select IND
        #### also have ID select in SB

        ## PlotOutput aln_scl -------------------------

      plotOutput("aln_scl", height = ind_shape_height)),
    column(
      width = 6,
      ## Plot of selected individual showing real and imputed data
      plotOutput("ind_missing",
                 height = ind_shape_height)
      )),

    fluidRow(
      column(width = 6,
             ## Explanotory text of the above plot
             "Selected individual is highlighted in green."
      ),
      column(width = 6,
             ## explanotry text of missing data
             "Selected individual is highlighted in red circles. Empty circles indicate the imputed value for missing data for that indivdual. Data imputation based on k-nearest neighbors."
             )
    ),

    fluidRow(
      column(
        width = 6,
        verbatimTextOutput("aln_summary")),
      column(
        width = 6,
        "n-missing data (across individuals):",
        tableOutput("n_miss_summary")
    )
  )) ## close MP
)), ## close sidebarlayout; close alingnment tab

# tab shapespace ----------------------------------

tabPanel("Shapespace", sidebarLayout(
  sidebarPanel(
    width = 2,

    numericInput(
      "x_pc",
      "PC on X axis",
      value = 1,
      min = 1,
      max = 20
    ),

    numericInput(
      "y_pc",
      "PC on Y axis",
      value = 2,
      min = 2,
      max = 20
    )

  ),
  ## close SB

  mainPanel(## PC shapes and scatter ---------------------------------

            ## switch this entire plot to a renderUI
            ## include an if/else that checks the output window size
            fluidRow(
              ## entire shapespace plot
              column(
                ## open LH column
                width = 8,

                ## Top row
                fluidRow(column(
                  width = 4,
                  offset = 4,
                  ### plotOutput pc_y_max -----------------------
                  plotOutput("pc_y_max", height = pc_shape_height, width = pc_shape_height *
                               2)
                )), ## close center column; ; close top row


                ## mid row

                fluidRow(
                  column(width = 4, fluidRow(plotOutput("pad_l", height = pc_pad)), fluidRow(
                    plotOutput("pc_x_min", height = pc_shape_height, width = pc_shape_height *
                                 2)
                  )),

                  column(
                    width = 4,
                    plotOutput(
                      "pc_scatter",
                      click = "pc_click",
                      height = pc_shape_height * 2,
                      width = pc_shape_height * 2
                    )
                  ),

                  column(width = 4, fluidRow(plotOutput("pad_r", height = pc_pad)), fluidRow(
                    plotOutput("pc_x_max", height = pc_shape_height, width = pc_shape_height *
                                 2)
                  ))

                ),

                ## bottom row
                fluidRow(column(
                  width = 4,
                  offset = 4,
                  plotOutput("pc_y_min", height = pc_shape_height, width = pc_shape_height *
                               2)
                ))

                ),
                ## close LH column

                ## Open RH column

                column(
                  width = 4,
                  fluidRow(plotOutput("target_shp", height = ind_shape_height)),
                  fluidRow(
                    "Individual profiles shown in grey (above); the hypothetical shape of selected coordinates is visualized in green.",
                    "The mean shape of the data is located at coordindates (0,0).",
                    "Selected coordinates:",
                    textOutput("target_coords")
                  )

                ) ## close RH column

              ), ## close shapespace monster

              fluidRow(## avg distance from centroid, x, y, xy
                ))
                ## close mainpanel
              )), ## close sidebarLAYOUT; close shapespace tab

# tab subgroups ------------------------------------------------
tabPanel("Subgroups", sidebarLayout(
  sidebarPanel(
    width = 2,

    uiOutput("dyn_k_grps"),

    # radioButtons(
    #   "cluster_type",
    #   label = "Clustering Algorithm:",
    #   choices = c("Hierarchical", "k-means")
    # ),
    checkboxInput("show_diagnostics", "Show Diagnostic Plots")

  ),
  ## close SB

  mainPanel(fluidRow(
    column(
      width = 6,
      plotOutput("dendro_plot",
                 height = ind_shape_height*2,
                 dblclick = "cut",
                 brush = "inds")),
    column(
      width = 6,
      fluidRow(
      plotOutput("all_shp1",
                 height = ind_shape_height),
      plotOutput("tree_shp",
                 height = ind_shape_height))
    )
  ), ## close row
  fluidRow(
    column(width = 6,
           "Dendrogram using Ward's clustering to visualize multivaraite relationships of PC scores. Individuals reprepsented by points across the bottom. Smaller branches represent more similar individuals. Use the controls at the left to split the dendrogram into k-clusters, represented by the plot color."
           ),
    column(width = 6,
           "Top, individual profiles amd mean for whole mean.",
           "<br>",
           "Bottom, individual profiles and mean for selected sub-set. Click and drag on the dendrogram to select a subset of individuals to visualize above. Selected individuals are highlted in orange above and on the dendrogram."
           )

  ),


  ## diagnostics row
  fluidRow(
    column(width = 6,
    plotOutput("sse_plot",
               height = ind_shape_height*2)),
    column(width = 6,
    plotOutput("sil_plot",
               height = ind_shape_height*2))
    )
  ) ## close MP
)), ## close sbp layout; close clustering tab

# tab phenotypes ---------------------------------------

tabPanel(
  "Phenotypes",
  sidebarLayout(
    sidebarPanel(
      width = 2
      ),
    mainPanel(
      fluidRow(
        column(width = 8,
          plotOutput("phenotypes")),
        column(width = 4,
               verbatimTextOutput("pheno_summary")
               ## display summary info for each group
               #### n sample for each subgroup
               #### avg + sd of x and y vals
        ))
      ) ## close MP
)) ## close sbp layout; close phenotypes tab

) ## Close TABSET layout
) ## close UI
