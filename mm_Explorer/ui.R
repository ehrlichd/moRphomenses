



library(shiny)
library(dplyr)
library(moRphomenses)

pc_shape_height <- 78
ind_shape_height <- 156

pc_pad <- 40



# Application title ------------------------

ui <- fluidPage(
titlePanel("mm Explorer"),

# tabset layout ---------------------------
tabsetPanel(

# tab Load Data --------------------------
tabPanel("Load Data", sidebarLayout(
  sidebarPanel(
    width = 2,
    fileInput("file_select","Choose data file")

  ), ## close SB
  mainPanel(
    uiOutput("dyn_var_assign"),
    actionButton("go", "Analyze!")
  ) ## close MP
  ) ## close sbpLayout
  ), ## close tab load-data

# tab alignment -------------------------
tabPanel("Alignment", sidebarLayout(
  sidebarPanel(
    width = 2,


    ## input$tar_alignment -----------------
    numericInput("tar_alignment",
                 label = "Align by day:",
                 value = 16,
                 min = 1,
                 max = 45
                 ),

    ## input$tar_length --------------------
    numericInput("tar_length",
      "Resample cycles to l-days",
      value = 28,
      min = 10,
      max = 50
      ),

    ## input$knn_n ---------------------------------
    numericInput("knn_n",
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
      column(width = 6,

             ## plot of all individuals
             ## hover to highlight a single line
             ## click to select IND
             #### also have ID select in SB

            ## PlotOutput aln_scl -----------------------------
            plotOutput("aln_scl",
                       height = ind_shape_height
                       )
            ),
      column(width = 6,

             ## Plot of selected individual showing real and imputed data
             plotOutput("ind_missing",
                        height = ind_shape_height)
             )
      ),
    fluidRow(
      column(width = 8,
             verbatimTextOutput("aln_summary")
             ),
      column(width = 8,
             "n-missing data (across individuals):",
             tableOutput("n_miss_summary")
             )
      )
    ) ## close MP
  ) ## close sidebarlayout
),## close alingment tab




# tab shapespace ----------------------------------

tabPanel("Shapespace", sidebarLayout(
  sidebarPanel(
    width = 2,

    numericInput("x_pc", "PC on X axis",
                 value = 1,
                 min = 1,
                 max = 20
                 ),

    numericInput("y_pc", "PC on Y axis",
                 value = 2,
                 min = 2,
                 max = 20
                 )

    ), ## close SB

  mainPanel(

    ## PC shapes and scatter ---------------------------------
    fluidRow( ## entire shapespace plot
      column( ## open LH column
        width = 8,

        ## Top row
        fluidRow(
          column(
            width = 4, offset = 4,
            ### plotOutput pc_y_max -----------------------
            plotOutput("pc_y_max",
                       height = pc_shape_height,
                       width = pc_shape_height*2)
            ) ## close center column
            ), ## close top row

        ## mid row

        fluidRow(
          column(
            width = 4,
            fluidRow(
              plotOutput("pad", height = pc_pad)),
            fluidRow(
                plotOutput("pc_x_min",
                           height = pc_shape_height,
                           width = pc_shape_height*2)
                )),

          column(
            width = 4,
            plotOutput("pc_scatter",
                       click = "pc_click",
                       height = pc_shape_height*2,
                       width = pc_shape_height*2)
            ),

          column(
            width = 4,
            fluidRow(
              plotOutput("pad", height = pc_pad)
              ),
            fluidRow(
              plotOutput("pc_x_max",
                         height = pc_shape_height,
                         width = pc_shape_height*2)
              ))

          ),

        ## bottom row
        fluidRow(
          column(
            width = 4,
            offset = 4,
            plotOutput("pc_y_min",
                       height = pc_shape_height,
                       width = pc_shape_height*2)
            ))

        ), ## close LH column

      ## Open RH column

      column(
        width = 4,
        fluidRow(
          plotOutput("target_shp",
                     height = ind_shape_height)
          ),
        fluidRow(
          "Individual profiles shown in grey (above); the hypothetical shape of selected coordinates is visualized in green.",
          "The mean shape of the data is located at coordindates (0,0).",
          "Selected coordinates:",
          textOutput("target_coords")
          )

        ) ## close RH column

      ), ## close shapespace monster

    fluidRow(
      ## avg distance from centroid, x, y, xy
      )

    ), ## close mainpanel
  ) ## close sidebarLAYOUT
), ## close tab shapespace

# tab subgroups ------------------------------------------------
tabPanel("Subgroups", sidebarLayout(
  sidebarPanel(
    width = 2,

    uiOutput("dyn_k_grps"),

    radioButtons("cluster_type", label = "Clustering Algorithm:",
                 choices = c("Hierarchical", "k-means")),
    checkboxInput("show_diagnostics", "Show Diagnostic Plots")

    ), ## close SB

  mainPanel(
  fluidRow(
    column(
      width = 6,

      plotOutput("dendro_plot", click = "cut", brush = "inds")
      ),

    column(
      width = 6,
      plotOutput("tree_shp", height = ind_shape_height)
      ) ## close RH
    ), ## close row

  ## diagnostics row
  fluidRow(
    uiOutput(
      plotOutput("sse_plot"),
      plotOutput("sil_plot")

    )
    )

) ## close MP
) ## close sbp layout

), ## close clustering

# tab phenotypes ---------------------------------------

tabPanel("Phenotypes", sidebarLayout(
  sidebarPanel(
    width = 2, ## do we need any parameters for this?
    ),


  mainPanel(
    fluidRow(
      column(width = 9,
             plotOutput("phenotypes")),

      column(width = 3,
             ## display summary info for each group
             #### n sample for each subgroup
             #### avg + sd of x and y vals
      )

      )
    ) ## close MP
) ## close sbp layout
) ## close phenotypes tab

) ## Close TABSET layout
) ## close UI
