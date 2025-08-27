effect_plot <- function(model, pred, pred.values = NULL, centered = "all",
                        plot.points = FALSE, interval = FALSE, data = NULL, at = NULL,
                        int.type = c("confidence","prediction"), int.width = .95,
                        outcome.scale = "response", robust = FALSE, cluster = NULL, vcov = NULL, 
                        set.offset = 1, x.label = NULL, y.label = NULL, pred.labels = NULL,
                        main.title = NULL, colors = "black", line.colors = colors, 
                        line.thickness = 1.1, 
                        point.size = 1.5, point.alpha = 0.6, jitter = 0, rug = FALSE, 
                        rug.sides = "lb", force.cat = FALSE, cat.geom = c("point", "line", "bar"), 
                        cat.interval.geom = c("errorbar", "linerange"), cat.pred.point.size = 3.5, 
                        partial.residuals = FALSE, color.class = colors, facet.by = NULL, ...) {
  pacman::p_load(ggplot2, rlang)
  source("scripts/plotting_theme.R")
  # Evaluate the pred arg
  pred <- as_name(enquo(pred))
  
  # Get the data right now rather than checking for its presence several times
  if (is.null(data)) {
    data <- get_data(model)
  }
  
  # Evaluate the facet.by arg, knowing it may be NULL
  facet.by <- enquo(facet.by) 
  if (!quo_is_null(facet.by)) { # check for NULL
    facet.by <- as_name(facet.by)
    if (is.null(at) || facet.by %nin% names(at)) {
      # If user isn't telling me the levels, then grab them all
      if (facet.by %nin% names(data)) {
        # Assume issue is variable isn't in the model formula for some reason
        the_formula <- stats::formula(model)
        the_formula <- update(the_formula, paste(". ~ . +", facet.by))
        data <- get_data(model, formula = the_formula)
      }
      at[[facet.by]] <- unique(data[[facet.by]])
      if (length(at[[facet.by]]) > 10) {
        msg_wrap(facet.by, " has ", length(at[[facet.by]]), " levels. This may
                 result in a difficult-to-see plot and/or a slow loading time
                 while the plot is generated. If you'd like to see just a 
                 subset of the levels, you can specify them in the `at` 
                 argument (e.g., at = list(", facet.by, " = c(",
                 dput(at[[facet.by]][1]), ", ", dput(at[[facet.by]][2]), ", ",
                 dput(at[[facet.by]][3]), "))", brk = "")
      }
    }
  } else {facet.by <- NULL} # don't make it a quosure anymore
  
  # Have a sensible interval default for categorical predictors
  if ("interval" %nin% names(match.call())[-1] && !(is.numeric(data[[pred]]) &&
                                                    force.cat == FALSE)) {
    interval <- TRUE
  }
  
  if (force.cat == TRUE && is.null(pred.values)) {
    pred.values <- sort(unique(suppressMessages(data[[pred]])))
  }
  
  # Deal with legacy color argument
  if (!all(color.class == colors)) colors <- color.class
  
  colors <- get_colors(colors)
  
  pred_out <- make_predictions(model, pred = pred, pred.values = pred.values,
                               at = at, center = centered,
                               interval = interval, int.type = int.type,
                               int.width = int.width,
                               outcome.scale = outcome.scale, robust = robust,
                               cluster = cluster, vcov = vcov,
                               set.offset = set.offset, return.orig.data = TRUE,
                               partial.residuals = partial.residuals, 
                               data = data, ...)
  
  # Putting these outputs into separate objects
  pm <- pred_out[[1]]
  d <- pred_out[[2]]
  
  # Check for clashing options "dpar" and "plot.points"
  dots <- list(...)
  if (!is.null(dots$dpar) && plot.points == TRUE) {
    plot.points <- FALSE
    warn_wrap("The plot.points argument is not compatible with distributional
              parameters specified in `dpar`.")
  }
  if (!is.null(dots$point.color)) {
    warn_wrap("The 'point.color' argument is deprecated and is now ignored. 
               You can change the color of points with the 'colors' argument.")
  }
  
  if (is.numeric(d[[pred]]) && force.cat == FALSE) {
    plot_effect_continuous(predictions = pm, pred = pred,
                           plot.points = plot.points | partial.residuals,
                           interval = interval,
                           data = d, x.label = x.label, y.label = y.label,
                           pred.labels = pred.labels, main.title = main.title,
                           colors = line.colors,
                           line.thickness = line.thickness, jitter = jitter,
                           resp = get_response_name(model, ...),
                           weights = get_weights(model, d)$weights_name,
                           rug = rug, rug.sides = rug.sides,
                           point.size = point.size, point.alpha = point.alpha,
                           point.color = colors, facet.by = facet.by)
  } else {
    plot_cat(predictions = pm, pred = pred, data = d,  
             geom = cat.geom, pred.values = pred.values,
             interval = interval, plot.points = plot.points | partial.residuals,
             pred.labels = pred.labels, x.label = x.label,
             y.label = y.label, main.title = main.title,
             colors = line.colors, weights = get_weights(model, d)$weights_name,
             resp = get_response_name(model, ...), jitter = jitter, 
             interval.geom = cat.interval.geom, line.thickness = line.thickness,
             point.size = point.size, pred.point.size = cat.pred.point.size,
             point.alpha = point.alpha, point.color = colors, facet.by = facet.by)
  }
  
}

plot_effect_continuous <- 
  function(predictions, pred, plot.points = FALSE, interval = FALSE, 
           data = NULL, x.label = NULL, y.label = NULL, pred.labels = NULL,
           main.title = NULL, colors = NULL, line.thickness = 1.1,
           jitter = 0.1, resp = NULL, weights = NULL, rug = FALSE,
           rug.sides = "b",
           point.size = 1, point.alpha = 0.6, point.color = "black", 
           facet.by = NULL) {
    
    pm <- predictions
    d <- data
    
    if (is.null(x.label)) {
      x.label <- pred
    }
    
    if (is.null(y.label)) {
      y.label <- resp
    }
    
    # If only 1 jitter arg, just duplicate it
    if (length(jitter) == 1) {jitter <- rep(jitter, 2)}
    
    # Starting plot object
    p <- ggplot(pm, aes(x = .data[[pred]], y = .data[[resp]]))
    
    # Plot observed data — do this first to plot the line over the points
    if (plot.points == TRUE) {
      
      constants <- list(alpha = point.alpha, colour = point.color)
      if (is.null(weights)) {
        # Only use constant size if weights are not used
        constants$size <- point.size
      } 
      # Need to use layer function to programmatically define constant aesthetics
      p <- p + layer(geom = "point", data = d, stat = "identity",
                     inherit.aes = FALSE, show.legend = FALSE,
                     mapping = aes(x = .data[[pred]], y = .data[[resp]], 
                                   size = .data[[weights]]),
                     position = position_jitter(width = jitter[1], 
                                                height = jitter[2]),
                     params = constants) +
        scale_size(range = c(1 * point.size, 5 * point.size))
      
    }
    
    # Define line thickness
    p <- p + geom_path(linewidth = line.thickness, colour = colors)
    
    # Plot intervals if requested
    if (interval == TRUE) {
      p <- p + geom_ribbon(data = pm, 
                           aes(ymin = .data$ymin,  ymax = .data$ymax),
                           alpha = 1/5, show.legend = FALSE, fill = colors)
    }
    
    # Rug plot for marginal distributions
    if (rug == TRUE) {
      p <- p + geom_rug(data = d,
                        mapping = aes(x = .data[[pred]], y = .data[[resp]]), alpha = 0.6,
                        position = position_jitter(width = jitter[1]),
                        sides = rug.sides, inherit.aes = TRUE, color = colors)
    }
    
    # Using theme_apa for theming...but using legend title and side positioning
    p <- p + labs(x = x.label, y = y.label) # better labels for axes
    
    # Getting rid of tick marks for two-level predictor
    if (length(unique(d[[pred]])) == 2) { # Predictor has only two unique values
      # Make sure those values are in increasing order
      brks <- sort(unique(d[[pred]]), decreasing = FALSE)
      if (is.null(pred.labels)) {
        p <- p + scale_x_continuous(breaks = brks)
      } else {
        if (length(pred.labels) == 2) { # Make sure pred.labels has right length
          p <- p + scale_x_continuous(breaks = brks, labels = pred.labels)
        } else {
          warning("pred.labels argument has the wrong length. It won't be used")
          p <- p + scale_x_continuous(breaks = brks)
        }
      }
    }
    
    if (!is.null(facet.by)) {
      p <- p + facet_wrap(facet.by)
    }   
    
    # Give the plot the user-specified title if there is one
    if (!is.null(main.title)) {
      p <- p + ggtitle(main.title)
    }
    
    # Return the plot
    return(p)
    
    
  }

plot_cat <- function(predictions, pred, data = NULL, 
                     geom = c("point", "line", "bar", "boxplot"), pred.values = NULL,
                     interval = TRUE, plot.points = FALSE, pred.labels = NULL, x.label = NULL,
                     y.label = NULL, main.title = NULL, colors = "black", weights = NULL,
                     resp = NULL, jitter = 0.1, geom.alpha = NULL, dodge.width = NULL,
                     errorbar.width = NULL, interval.geom = c("errorbar", "linerange"),
                     line.thickness = 1.1, point.size = 1, pred.point.size = 3.5, 
                     point.alpha = 0.6, point.color = "black", facet.by = NULL) {
  
  pm <- predictions
  d <- data
  
  geom <- geom[1]
  if (geom == "dot") {geom <- "point"}
  
  # If only 1 jitter arg, just duplicate it
  if (length(jitter) == 1) {jitter <- rep(jitter, 2)}
  
  if (is.null(x.label)) {
    x.label <- pred
  }
  
  if (is.null(y.label)) {
    y.label <- resp
  }
  
  # Deal with numeric predictors coerced into factors
  if (is.numeric(pm[[pred]])) {
    pred.levels <- if (!is.null(pred.values)) {pred.values} else {
      unique(pm[[pred]])
    }
    pred.labels <- if (!is.null(pred.labels)) {pred.labels} else {
      unique(pm[[pred]])
    }
    pm[[pred]] <- factor(pm[[pred]], levels = pred.levels,
                         labels = pred.labels)
    
    # Make sure only observations of requested levels of predictor are included
    d <- d[d[[pred]] %in% pred.levels,]
    d[[pred]] <- factor(d[[pred]], levels = pred.levels, labels = pred.labels)
  }
  
  # Convert strings to symbols for tidy evaluation
  pred <- sym(pred)
  resp <- sym(resp)
  if (!is.null(weights)) {weights <- sym(weights)}
  
  # Checking if user provided the colors his/herself
  colors <- suppressWarnings(get_colors(colors, 1))
  
  if (is.null(geom.alpha)) {
    a_level <- 1
    if (plot.points == TRUE) {
      a_level <- 0.5
    } else if (interval == TRUE) {
      a_level <- 0.5
    }
  } else {a_level <- geom.alpha}
  
  if (is.null(dodge.width)) {
    dodge.width <- if (geom %in% c("bar", "point")) {0.9} else {0}
  }
  if (is.null(errorbar.width)) {
    errorbar.width <- if (geom == "point") {
      0.75
    } else if (geom == "bar") {
      0.75
    } else {0.5}
  }
  
  
  p <- ggplot(pm, aes(x = .data[[pred]], y = .data[[resp]], group = 1))
  
  if (geom == "bar") {
    p <- p + geom_bar(stat = "identity", position = "dodge", alpha = a_level,
                      show.legend = FALSE, color = colors)
  } else if (geom %in% c("point", "line")) {
    p <- p + geom_point(size = pred.point.size,
                        position = position_dodge(dodge.width),
                        show.legend = FALSE, color = colors)
  }
  
  if (geom == "line") {
    p <- p + geom_path(position = position_dodge(dodge.width),
                       linewidth = line.thickness, show.legend = FALSE, 
                       color = colors)
  }
  
  # Plot intervals if requested
  if (interval == TRUE && interval.geom[1] == "errorbar") {
    p <- p + geom_errorbar(aes(ymin = .data$ymin, ymax = .data$ymax),
                           alpha = 1, show.legend = FALSE,
                           position = position_dodge(dodge.width),
                           width = errorbar.width,
                           linewidth = line.thickness, color = colors)
  } else if (interval == TRUE && interval.geom[1] %in% c("line", "linerange")) {
    p <- p + geom_linerange(aes(ymin = .data$ymin, ymax = .data$ymax),
                            alpha = 0.8, show.legend = FALSE,
                            position = position_dodge(dodge.width),
                            linewidth = line.thickness, color = colors)
  }
  
  # For factor vars, plotting the observed points
  # and coloring them by factor looks great
  if (plot.points == TRUE) {
    
    constants <- list(alpha = point.alpha, colour = point.color)
    if (is.null(weights)) {
      # Only use constant size if weights are not used
      constants$size <- point.size
    } 
    # Need to use layer function to programmatically define constant aesthetics
    p <- p + layer(geom = "point", data = d, stat = "identity",
                   inherit.aes = FALSE, show.legend = FALSE,
                   mapping = aes(x = .data[[pred]], y = .data[[resp]], 
                                 size = .data[[weights]]),
                   position = position_jitter(width = jitter[1], 
                                              height = jitter[2]),
                   params = constants) +
      scale_size(range = c(1 * point.size, 5 * point.size))
    
  }
  
  # Using theme_apa for theming...but using legend title and side positioning
  p <- p + drop_x_gridlines() + 
    labs(x = x.label, y = y.label) # better labels for axes
  
  if (!is.null(facet.by)) {
    p <- p + facet_wrap(facet.by)
  }   
  
  # Give the plot the user-specified title if there is one
  if (!is.null(main.title)) {
    p <- p + ggtitle(main.title)
  }
  
  # Return the plot
  return(p)
  
}