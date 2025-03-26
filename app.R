# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Make sure packages are installed!

required_packages <- c("shiny", "shinyjs", "devtools", "stats", "mclust", "nortest", "mappeR")


install_if_missing <- function(packages) {
  missing_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  }

  # Load all packages
  for (package in packages) {
    library(package, character.only = TRUE)
  }
}

# install/load all required packages
install_if_missing(required_packages)

# stolen from sir jacob miller
generate_spiral <- function(n=1000, noise=0.1){
  t <- runif(n, 0, 4 * pi)
  h <- runif(n, 0, 1)

  hole_mask <- (t > 1.5*pi & t < 1.75*pi) & (h>0.25 & h<0.75)
  t <- t[!hole_mask]
  h <- h[!hole_mask]

  x <- t*cos(t) + rnorm(length(t), 0, noise)
  z <- t*sin(t) + rnorm(length(t), 0, noise)

  return(data.frame(x=x, y=z))
}

generate_barbell <- function(n) {
  r1 = sqrt(runif(n*.45, 0, .5))
  a1 = runif(n*.45, 0, 2*pi)
  disk1 = data.frame(x=r1*cos(a1) - 1, y=r1*sin(a1))

  r2 = sqrt(runif(n*.45, 0, .5))
  a2 = runif(n*.45, 0, 2*pi)
  disk2 = data.frame(x=r2*cos(a2) + 1, y=r2*sin(a2))

  line = data.frame(x = runif(n/20, -.5, .5), y = rep(0, n/20))

  return(rbind(disk1, disk2, line))
}

# Define UI for application that constructs mapper graph

ui <- fluidPage(
    useShinyjs(),
    titlePanel("1D Mapper"),

    # Sidebar with parameter input options
    sidebarLayout(
        sidebarPanel(
            selectInput(
                "data",
                "Dataset",
                choices = c("circle", "figure 8", "spiral", "barbell")
            ),
            sliderInput(
                "points",
                "Number of points",
                value = 1000,
                min = 100,
                max = 2000,
                step = 100
            ),
            selectInput(
                "cover_type",
                "Covering Scheme:",
                choices = c("Width-Balanced", "GMM-based", "G-Mapper")
            ),
            selectInput(
                "lens",
                "Lens Function: ",
                choices = c("project to x", "project to y", "use eccentricity value", "PCA-1")
            ),

            # we only show bins slider for Width-Balanced and GMM-based
            conditionalPanel(
                condition = "input.cover_type != 'G-Mapper'",
                sliderInput(
                    "bins",
                    "Number of bins:",
                    min = 1,
                    max = 50,
                    value = 10
                )
            ),

            # G-Mapper parameters
            conditionalPanel(
                condition = "input.cover_type == 'G-Mapper'",
                sliderInput(
                    "alpha",
                    "Alpha (significance level):",
                    min = 0.001,
                    max = 0.2,
                    value = 0.05,
                    step = 0.001
                ),
                sliderInput(
                    "min_points_percent",
                    "Minimum points per interval (% of total):",
                    min = 1,
                    max = 50,
                    value = 10,
                    step = 0.5
                ),
                sliderInput(
                    "max_iterations",
                    "Maximum iterations:",
                    min = 5,
                    max = 50,
                    value = 20
                ),
                htmlOutput("gmapper_bins_info")
            ),
            
            sliderInput(
                "percent_overlap",
                "Percent overlap (weight):",
                min = 0,
                max = 100,
                value = 25
            ),

            selectInput(
                "method",
                "Clustering method",
                choices = c("single", "complete", "average", "ward.D2", "mcquitty")
            )
        ),

        # plot mapper graph
        mainPanel(plotOutput("inputdata"), plotOutput("mapper"))
    )
)


# Define server logic required to construct mapper graph
server <- function(input, output, session) {
    # generate sample data
    data = reactive({
        switch(
            input$data,
            "circle" = data.frame(
                x = sapply(1:input$points, cos) + runif(input$points, 0, .1),
                y = sapply(1:input$points, sin) + runif(input$points, 0, .1)
            ),
            "figure 8" = data.frame(
                x = sapply(1:input$points, function(x) cos(x) / (1 + sin(x)^2)) + runif(input$points, 0, .1),
                y = sapply(1:input$points, function(x) sin(x)*cos(x) / (1 + sin(x)^2)) + runif(input$points, 0, .1)
            ),
            "spiral" = generate_spiral(input$points),
            "barbell" = generate_barbell(input$points)
        )
    })

    # filter data
    filtered_data = reactive({
      # grab current data
      data = data()

      switch(input$lens,
             "project to x" = data$x,
             "project to y" = data$y,
             "use eccentricity value" = eccentricity_filter(data),
             "PCA-1" = {
                pca_output <- prcomp(data, center = FALSE, scale. = FALSE)
                pca_output$x[,1]
             })
    })

# A reactive expression for the 'G-Mapper' covering scheme
# Should now only compute the model when needed

    # number of bins display for G-Mapper setting
    output$gmapper_bins_info <- renderUI({
        if(input$cover_type == "G-Mapper" && !is.null(gmapper_model())) {
            bin_count <- nrow(gmapper_model()$intervals)
            HTML(paste0(
                "<div style='margin-top: 20px; margin-bottom: 20px; padding: 10px; background-color: #6ae106; border-radius: 10px;'>",
                "<strong>G-Mapper determined ", bin_count, " bins</strong>",
                "</div>"
            ))
        }
    })

    gmapper_model <- reactive({
        if(input$cover_type == "G-Mapper") {

            # G-Means algorithm for cover generation

            g_means_cover <- function(filtered_data, alpha, min_points_percent, max_iterations) {
                # Calculate minimum points based on percentage of total data
                total_points <- length(filtered_data)
                min_points <- max(30, round(total_points * min_points_percent / 100))

                # initialize with a single interval
                data_range <- range(filtered_data)
                intervals <- matrix(data_range, nrow=1, ncol=2)
                means <- mean(filtered_data)
                variances <- var(filtered_data)

                # Keeps track if intervals need testing
                to_test <- TRUE
                iteration <- 0 # will eventually halt

                while(any(to_test) && iteration < max_iterations) {
                    iteration <- iteration + 1
                    current_intervals <- intervals[to_test, , drop=FALSE]
                    new_intervals <- list()
                    new_means <- c()
                    new_variances <- c()
                    new_to_test <- c()

                    for(i in 1:nrow(current_intervals)) {
                        interval <- current_intervals[i, ]
                        # Get data in this interval
                        interval_data <- filtered_data[filtered_data >= interval[1] & filtered_data <= interval[2]]

                        if(length(interval_data) < min_points) { # Not enough points to test
                            new_intervals[[length(new_intervals) + 1]] <- interval
                            new_means <- c(new_means, mean(interval_data))
                            new_variances <- c(new_variances, var(interval_data))
                            new_to_test <- c(new_to_test, FALSE)
                            next
                        }

                        # Normalize data for AD
                        normalized_data <- (interval_data - mean(interval_data)) / sqrt(var(interval_data))

                        # Anderson-Darling test for normality
                        ad_test <- ad.test(normalized_data)

                        if(ad_test$p.value < alpha) { # Not normal => split the interval
                            # 2-component GMM
                            data_matrix <- as.matrix(interval_data)
                            gmm <- Mclust(data_matrix, G=2, modelNames="V")

                            if(!is.null(gmm) && length(gmm$parameters$mean) == 2) {
                                # get the means and sort
                                gmm_means <- gmm$parameters$mean
                                gmm_vars <- gmm$parameters$variance$sigmasq
                                sorted_idx <- order(gmm_means)
                                gmm_means <- gmm_means[sorted_idx]
                                gmm_vars <- gmm_vars[sorted_idx]

                                # Split point (weighted by variances)
                                split_point <- (gmm_means[1] * sqrt(gmm_vars[2]) + gmm_means[2] * sqrt(gmm_vars[1])) / 
                                                (sqrt(gmm_vars[1]) + sqrt(gmm_vars[2]))

                                # split into two new intervals in cover
                                new_intervals[[length(new_intervals) + 1]] <- c(interval[1], split_point)
                                new_intervals[[length(new_intervals) + 1]] <- c(split_point, interval[2])

                                # find means and variances for the new intervals
                                left_data <- interval_data[interval_data <= split_point]
                                right_data <- interval_data[interval_data > split_point]

                                new_means <- c(new_means, mean(left_data), mean(right_data))
                                new_variances <- c(new_variances, var(left_data), var(right_data))

                                # mark new intervals for testing
                                new_to_test <- c(new_to_test, TRUE, TRUE)
                            } else {
                                # if GMM fitting fails, keep the original interval
                                new_intervals[[length(new_intervals) + 1]] <- interval
                                new_means <- c(new_means, mean(interval_data))
                                new_variances <- c(new_variances, var(interval_data))
                                new_to_test <- c(new_to_test, FALSE)
                            }
                        } else {
                            # if AD determines ~ normal, keep the interval
                            new_intervals[[length(new_intervals) + 1]] <- interval
                            new_means <- c(new_means, mean(interval_data))
                            new_variances <- c(new_variances, var(interval_data))
                            new_to_test <- c(new_to_test, FALSE)
                        }
                    }

                    # Update intervals that weren't tested
                    if(any(!to_test)) {
                        untested_intervals <- intervals[!to_test, , drop=FALSE]
                        untested_means <- means[!to_test]
                        untested_vars <- variances[!to_test]

                        # combine with new intervals
                        intervals <- rbind(do.call(rbind, new_intervals), untested_intervals)
                        means <- c(new_means, untested_means)
                        variances <- c(new_variances, untested_vars)
                        to_test <- c(new_to_test, rep(FALSE, sum(!to_test)))
                    } else {
                        intervals <- do.call(rbind, new_intervals)
                        means <- new_means
                        variances <- new_variances
                        to_test <- new_to_test
                    }

                    # intervals sorted by LOWER bound
                    sort_idx <- order(intervals[,1])
                    intervals <- intervals[sort_idx, , drop=FALSE]
                    means <- means[sort_idx]
                    variances <- variances[sort_idx]
                    to_test <- to_test[sort_idx]
                }

                # results returned as a list
                return(list(
                    intervals = intervals,
                    means = means,
                    variances = variances
                ))
            }

            # run G-Means
            g_means_result <- g_means_cover(
                filtered_data(), 
                alpha = input$alpha, 
                min_points_percent = input$min_points_percent, 
                max_iterations = input$max_iterations
            )
            return(g_means_result)
        } else {
            NULL
        }
    })


    # updates the bins slider when G-Mapper model changes
    observe({
        if(input$cover_type == "G-Mapper" && !is.null(gmapper_model())) {
            n_components <- nrow(gmapper_model()$intervals)
            updateSliderInput(session, "bins", value = n_components)
        }
    })



    # generate cover
    cover = reactive({
    # grab current data
    data = data()

    # grab current filter values
    filtered_data = filtered_data()

    # create cover based on user selection
    switch(input$cover_type,
            "Width-Balanced" = create_width_balanced_cover(min(filtered_data),
                                max(filtered_data),
                                input$bins,
                                input$percent_overlap),
            "GMM-based" = {
                if (!require("mclust")) {
                    install.packages("mclust")
                    library(mclust)
                }
                # define the functions:
                fit_gmm <- function(filtered_data, n_components) {
                    data_matrix <- as.matrix(filtered_data)

                    # Fit a GMM to 1D-filtered values
                    # "V" model (variable variance, univariate)
                    gmm_model <- Mclust(data_matrix, G=n_components, modelNames="V")
                    return(gmm_model)
                }

                create_gmm_cover <- function(filtered_data, gmm_model, percent_overlap) {
                    # grab variables
                    means <- gmm_model$parameters$mean
                    variances <- gmm_model$parameters$variance$sigmasq
                    sorted_indices <- order(means)
                    sorted_means <- means[sorted_indices]
                    sorted_variances <- variances[sorted_indices]
                    intervals <- matrix(0, nrow=length(sorted_means), ncol=2)

                    # interval lengths weighted by std_deviations and user specified overlap
                    for (i in 1:length(sorted_means)) {
                        std_dev <- sqrt(sorted_variances[i])

                        intervals[i,] <- c(sorted_means[i] - 2*std_dev, sorted_means[i] + 2*std_dev)
                    }

                    for(j in 2:length(sorted_means)){
                        midpoint <- (sorted_means[j] + sorted_means[j-1])/2
                        overlap_distance <- percent_overlap*(sorted_means[j] - sorted_means[j-1])/100
                        intervals[j-1, 2] <- midpoint + overlap_distance
                        intervals[j, 1] <- midpoint - overlap_distance

                    }

                    intervals[1,1] <- min(filtered_data)
                    intervals[length(sorted_means), 2] <- max(filtered_data)

                    return(intervals)
                }
                # create the actual cover
                gmm_result <- fit_gmm(filtered_data, input$bins)
                create_gmm_cover(filtered_data, gmm_result, input$percent_overlap)
            },
            "G-Mapper" = {
                # pre-computed G-Means model
                model <- gmapper_model()

                # create proper overlap weighted similarly to GMM-based
                create_gmapper_cover <- function(model, percent_overlap) {
                    intervals <- model$intervals
                    means <- model$means
                    variances <- model$variances

                    # Adjust intervals for overlap
                    for(j in 2:length(means)) {

                        midpoint <- (means[j] + means[j-1])/2

                        sd_left <- sqrt(variances[j-1])
                        sd_right <- sqrt(variances[j])

                        r <- sd_right / sd_left
                        overlap_factor <- percent_overlap * (r * sd_left + sd_right) / (r + 1) / 100

                        intervals[j-1, 2] <- midpoint + overlap_factor
                        intervals[j, 1] <- midpoint - overlap_factor
                    }

                    return(intervals)
                }

                create_gmapper_cover(model, input$percent_overlap)
            }
    )
})

    # generate mapper graph
    mapper = reactive({
        # grab current data
        data = data()

        # grab current filter values
        filtered_data = filtered_data()

        # grab current cover
        cover = cover()

        # create mapper graph
        create_1D_mapper_object(
            data,
            dist(data),
            filtered_data,
            cover,
            hierarchical_clusterer(input$method)
        )
    })

    # output data plot
    output$inputdata <- renderPlot({
        data = data()
        filtered_data = filtered_data()
        cover = cover()

        # plot data
        plot(data, pch = 20)

        # define a function that creates a color gradient
        colorRampAlpha <- function(..., n, alpha) {
            colors <- colorRampPalette(...)(n)
            paste(colors, sprintf("%x", ceiling(255 * alpha)), sep = "")
        }

        # create a color gradient from blue to gold to red with (number of bins) colors
        bincolors = colorRampAlpha(c("blue", "gold", "red"),
                                 alpha = .5,
                                 n = input$bins)

        # plot the bins on top of the data
        switch(input$lens,
               "project to x" = rect(cover[, 1], min(data$y), cover[, 2], max(data$y), col = bincolors),

               "project to y" = rect(min(data$x), cover[, 2], max(data$x), cover[, 1], col = bincolors),

               "PCA-1" = {
                  # draw PCA line
                  pca_output <- prcomp(data, center = FALSE, scale. = FALSE)
                  pca_vector <- pca_output$rotation[,1]
                  slope <- pca_vector[2] / pca_vector[1]
                  abline(0, slope, col = "green", lwd = 3, lty = 3)

                  # calculate perpendicular vector
                  perp_vector <- c(-pca_vector[2], pca_vector[1])
                  perp_vector <- perp_vector / sqrt(sum(perp_vector^2)) # normalize

                  # color (super annoying) bins using a loop since polygon
                  for (i in 1:nrow(cover)) {
                          # calculate cut points for bins on the pca line
                          cut1 <- cover[i, 1] * pca_vector
                          cut2 <- cover[i, 2] * pca_vector
                          # the 100 is just to make sure the bins don't get cutoff in the image
                          corner1 <- cut1 + 100 * perp_vector
                          corner2 <- cut1 - 100 * perp_vector
                          corner3 <- cut2 - 100 * perp_vector
                          corner4 <- cut2 + 100 * perp_vector

                          # Draw filled rectangle using polygon since rect didn't work :((
                          polygon(x = c(corner1[1], corner2[1], corner3[1], corner4[1]),
                                  y = c(corner1[2], corner2[2], corner3[2], corner4[2]),
                                  col = bincolors[i])
                  }
                })
    })

    # output mapper graph
    output$mapper <- renderPlot({
        plot(mapper_object_to_igraph(mapper()))
    })
}

# Run the application
shinyApp(ui = ui, server = server)
