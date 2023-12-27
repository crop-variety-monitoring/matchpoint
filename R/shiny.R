
# RH based on code by Luis Mijangos

check_files <- function(snpfile, genfile) {
	if (is.null(snpfile)) return("empty SNP filename")
	if (is.null(genfile)) return("empty genotypes filename")
	if (!file.exists(snpfile)) return("SNP file does not exist")
	if (!file.exists(genfile)) return("genotype file does not exist")
	return("")
}

shiny_IBS <- function(...) {

	# globals
	recoded <- FALSE
	crd <- js <- NULL
	fsnp <- system.file("ex/DCas00-0000_SNP.csv", package="matchpoint")
	fvar <- system.file("ex/DCas00-0000_genotype-info.csv", package="matchpoint")

  #### dataset list ###
	datasetListUI <- function(id) {
		ns <- shiny::NS(id)
		shiny::uiOutput(ns('dataList'))
	}
	datasetListServer <- function(id, md_list) {
		shiny::moduleServer(id,
			function(input, output, session) {
				output$dataList = shiny::renderUI({
					shiny::selectizeInput(
						inputId = 'dataList',
						label = 'Select a sample',
						choices = md_list,
						options = list(placeholder = '')
					)
				})
			})
	}

	# maximum size of file to upload
	options(shiny.maxRequestSize = 50 * 1024 ^ 2)
	box_color <- "blue"
	box_height <- '100%'
	input_class <- c("ui small icon input", "ui fluid icon input")

#######################################################################
############### UI ####################################################
#######################################################################


	ui <- semantic.dashboard::dashboardPage(
		title = "Variety Identification",
		semantic.dashboard::dashboardHeader(
			color = "green",
			menu_button_label = "",
			class = "ui top attached header",
			shiny.semantic::button(
				input_id = "close",
				label = shiny::span(shiny::icon("close"), "Exit"),
				class = c("tiny", "ui red button", "compact ui button")
			),
		),

		### Sidebar content ###
		semantic.dashboard::dashboardSidebar(
			size = "thin",
			color = "green",
			semantic.dashboard::sidebarMenu(
				semantic.dashboard::menuItem(
					text = shiny::span(shiny::icon("upload"), "Inputs"), tabName = "inputs_tab"),
				semantic.dashboard::menuItem(
					text = shiny::span(shiny::icon("gear"), "Parameters"), tabName = "parameters_tab"
				),
				semantic.dashboard::menuItem(
					text = shiny::span(shiny::icon("filter"), "Results"), tabName = "results_tab"
				),
				semantic.dashboard::menuItem(
					text = shiny::span(shiny::icon("eye"), "Visualization"), tabName = "visualization_tab"
				),
				semantic.dashboard::menuItem(
					text = shiny::span(shiny::icon("map"), "Map"), tabName = "map_tab"
				)
			)
		),

		## Body content
		semantic.dashboard::dashboardBody(
			shinyjs::useShinyjs(),
			shinyjs::extendShinyjs(text = "shinyjs.closeWindow = function() { window.close(); }",
										functions = c("closeWindow")),
			shinybusy::add_busy_spinner(spin = "fading-circle"),

			shinybusy::add_busy_gif(
			  src = "https://jeroen.github.io/images/banana.gif",
			  height = 70, width = 70
			),

#### INPUT
			semantic.dashboard::tabItems(
				semantic.dashboard::tabItem(
					tabName = "inputs_tab",
					semantic.dashboard::box(
						title =	"Input files",
						color = box_color,
						width = 16,
						collapsible = FALSE,
						title_side = "top left",
						style = box_height,
						shiny.semantic::fileInput(
							inputId = "snp_file",
							label = "SNP file",
							buttonLabel = "Browse",
							type = input_class,
							placeholder = fsnp,
							accept = "csv"
						),
						shiny.semantic::fileInput(
							inputId = "gen_file",
							label = "Genotypes file",
							buttonLabel = "Browse",
							type = input_class,
							placeholder = fvar,
							accept = "csv"
						),
						shiny.semantic::button(
							input_id = "run_ref_check",
							label = shiny::span(shiny::icon("play"), "RUN"),
							class = "ui green button"
						)
					),
					shiny::fluidRow(
						id = "reference",
						shiny::column(3,
							shiny::verbatimTextOutput('rf_read'),
							shiny::verbatimTextOutput('rf_combine'),
							shiny::verbatimTextOutput('rf_recode')
						)
					)					
				),


#### PARAMA
				semantic.dashboard::tabItem(
					tabName = "parameters_tab",
					semantic.dashboard::box(
						title = "Params",
						width = 16,
						color = box_color,
						collapsible = FALSE,
						title_side = "top left",
						shiny::fluidRow(
							shiny::column(2,
								shiny::numericInput(
									inputId = "ibs",
									label = "IBS cut-off",
									value = .8,
									min = 0.1,
									max = 0.9
								),
								shiny::br(),
								shiny.semantic::button(
									input_id = "run_id",
									label = shiny::span(shiny::icon("play"), "RUN"),
									class = "ui green button"
								),
								shiny::br()
							)
						),
					)
				),

#### OUTPUT
				semantic.dashboard::tabItem(
					tabName = "results_tab",
					semantic.dashboard::box(
						title = "Results",
						width = 16,
						color = box_color,
						collapsible = FALSE,
						title_side = "top left",
						shiny::fluidRow(
							shiny::column(1,
								shiny::br(),
								shiny::downloadButton('download', "Save to .csv"),
								shiny::br()
							)
						),
						DT::dataTableOutput("res.ID"),
						style = "height:800px; overflow-y: scroll;overflow-x: scroll;"
					)
				),


############### VISUALIZATION #########################################

				semantic.dashboard::tabItem(
					tabName = "visualization_tab",
					semantic.dashboard::box(
						title = "Map",
						width = 16,
						color = box_color,
						collapsible = FALSE,
						title_side = "top left",
						shiny.semantic::flow_layout(
							shiny.semantic::button(
								input_id = "make_map",
								label = shiny::span(shiny::icon("play"), "RUN"),
								class = "ui green button"
							),
							shiny::selectInput("variety", "variety:", c("none" = "none"))
						),
						shiny::p(),
						leaflet::leafletOutput("mymap"),
						shiny::p()
					)
				)
			)
		)
	)
	
#######################################################################
############### SERVER ################################################
#######################################################################

	server <- function(input, output, session) {

		options(shiny.maxRequestSize=50*1024^2)
		ID_res <- shiny::reactiveVal(NULL)

		# Close button
		shiny::observeEvent(input$close, {
		  lapply(names(shiny::resourcePaths()), shiny::removeResourcePath)
		  js$closeWindow()
		  shiny::stopApp()
		})

#### INPUT

		shiny::observeEvent(input$run_ref_check, {

			snpfile = input$snp_file$datapath
			genfile = input$gen_file$datapath
#			snpfile = fsnp
#			genfile = fvar

			check = check_files(snpfile, genfile)
			if (check != "") {
				output$rf_read <- shiny::renderText({check})
				return(NULL)
			}

			snps <- try(matchpoint::read_dart(snpfile))
			if (inherits(snps, "try-error")) stop("cannot read SNP file")
			genotypes <- try(data.table::fread(genfile))
			if (inherits(genotypes, "try-error")) stop("cannot read genotype file")
			markers <- matchpoint::marker_positions("")

			geo <- c("dart.id", "longitude", "latitude")
			if (geo %in% names(genotypes)) {
				crd <<- genotypes[, geo]
			} else {
				crd <<- NULL
			}

			snpname = input$snp_file$name
			genname = input$geno_file$name		

			output$rf_read <- shiny::renderText({
				paste0("SNP        : ", snpname, ", ", nrow(snps$snp), " records\n", 
					   "genotypes  : ", genname, ", ", nrow(genotypes), " records\n",
					   "coordinates: ", ifelse(is.null(crd), "No", "Yes\n\n"))
			})
			crf <- match_IBS(snps$snp, genotypes, markers)

			output$rf_combine <- shiny::renderText({
				paste(paste0(names(crf), "\n"))
			})		
			recoded <<- TRUE
			output$rf_recode <- shiny::renderText({
				"recoded successfully"
			})

		})
############### REFERENCE IDENTIFICATION ##############################

		shiny::observeEvent(input$run_id, {

			if (!recoded) {
				output$rf_id <- shiny::renderText({
					"input has not been generated (go back one tab)"
				})
			
			} else {
				
				IBS_cutoff = input$ibs
				ibsvar <- paste0("IBS_cutoff_", IBS_cutoff, "_best_match")
				
				res_ID <<- crf$best_match
					
				ID_res(res_ID)

				res_ID$IBS <- round(res_ID$IBS, 3)
				res_ID$Sample_SNP_mr <- round(res_ID$Sample_SNP_mr, 3)
				
				output$res.ID <- DT::renderDataTable({DT::datatable(res_ID)})

				output$download <- shiny::downloadHandler(
					filename = function() {"match_results.csv"},
					content = function(fname) {
						utils::write.csv(res_ID$best_match, fname)
					}
				)
				
			}
		})


    ############### VISUALIZATION #########################################
		
		output$mymap <- leaflet::renderLeaflet({
			leaflet::leaflet() |> leaflet::setView(0, 0, zoom = 3) |>
			  leaflet::addProviderTiles("OpenStreetMap",
				options = leaflet::providerTileOptions(noWrap = TRUE)
			  ) # |> addMarkers(data = points())
		})

		drawPoints <- function(variety) {
			if (is.null(crd)) return(NULL)
			if (!all(c("longitude", "latitude") %in% colnames(crd))) return (NULL)
			points <- crd[ , c("longitude", "latitude")]
			if (variety == "all varieties") {
				leaflet::leafletProxy("mymap", data = points) |>
					leaflet::clearShapes() |>
					leaflet::addCircles(radius = 15, weight = 10, color = "red")				
			} else {
				leaflet::leafletProxy("mymap", data = points) |>
					leaflet::clearShapes() |>
					leaflet::addCircles(radius = 10, weight = 6, color = "gray")				
				points <- crd[crd[,1] == variety, c("longitude", "latitude")]
				leaflet::leafletProxy("mymap", data = points) |>
					leaflet::addCircles(radius = 15, weight = 10, color = "red")	
			}
		}
		
		shiny::observeEvent(input$make_map, {
			uvars <- c("all varieties", unique(crd[,1]))
			shiny::updateSelectInput(session, "variety",
				label = "variety",
				choices = uvars,
				selected = uvars[1]
			)
			#drawPoints(input$variety)
		})

		shiny::observeEvent(input$variety, {
			#drawPoints(input$variety)
		})

	} # server end
	
	shiny::shinyApp(ui, server, ...)
}



shiny_IBS()

