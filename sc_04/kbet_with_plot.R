#truncated normal distribution distribution function
ptnorm <- function(x, mu, sd, a=0, b=1, alpha = 0.05, verbose = FALSE){
  #this is the cumulative density of the truncated normal distribution
  #x ~ N(mu, sd^2), but we condition on a <= x <= b
  if (!is.na(x)){


    if (a > b) {
      warning("Lower and upper bound are interchanged.")
      tmp <- a
      a <- b
      b <- tmp
    }

    if (sd <= 0 || is.na(sd)) {
      if (verbose) {
        warning("Standard deviation must be positive.")
      }
      if (alpha <= 0) {
        stop("False positive rate alpha must be positive.")
      }
      sd <- alpha
    }
    if (x < a || x > b) {
      warning("x out of bounds.")
      cdf <- as.numeric(x > a)
    } else {
      alp <- pnorm((a - mu) / sd)
      bet <- pnorm((b - mu) / sd)
      zet <- pnorm((x - mu) / sd)
      cdf <- (zet - alp) / (bet - alp)
    }
    cdf
  } else {
    return(NA)
  }
}

# residual score function of kBET
residual_score_batch <- function(knn.set, class.freq, batch) {
    # knn.set: indices of nearest neighbors
    # class.freq: global batch proportions
    # batch: batch labels

    # return NA if all values of neighborhood are NA (which may arise from subsampling a knn-graph)
    if (all(is.na(knn.set))) {
        return(NA)
    }
    else{
        # extracts the batch labels of all neighbors (excluding NA)
        # computes local batch frequencies (observed)
        freq.env <- table(batch[knn.set[!is.na(knn.set)]]) / length(!is.na(knn.set))
        
        # create zero vector with length of batches and initialize it with freqs of batches
        full.classes <- rep(0, length(class.freq$class))
        full.classes[class.freq$class %in% names(freq.env)] <- freq.env
        
        # global batch props (expected)
        exp.freqs <- class.freq$freq

        # compute chi-square test statistics
        ## sum((observed - expected)^2 / expected)
        chi_sq_statistic <- sum((full.classes - exp.freqs)^2 / exp.freqs)
        
        return(chi_sq_statistic)
    }
}

# core function of kBET; called for each sample
chi_batch_test <- function(knn.set, class.freq, batch, df) {
    # knn.set: indices of nearest neighbors
    # class.freq: global batch proportions
    # batch: batch labels
    # df: degrees of freedom

    # return NA if all values of neighborhood are NA (which may arise from subsampling a knn-graph)
    if (all(is.na(knn.set))) {
        return(NA)
    }
    else {
        # extracts the batch labels of all neighbors (excluding NA)
        # computes local batch counts (observed)
        freq.env <- table(batch[knn.set[!is.na(knn.set)]])

        # create zero vector with length of batches and initialize it with counts of batches
        full.classes <- rep(0, length(class.freq$class))
        full.classes[class.freq$class %in% names(freq.env)] <- freq.env

        # global batch counts (for this sample)
        exp.freqs <- class.freq$freq * length(knn.set)

        # compute chi-square test statistics
        ## sum((observed - expected)^2 / expected)
        chi.sq.value <- sum((full.classes - exp.freqs)^2 / exp.freqs)

        # calculate p-value
        result <- 1 - pchisq(chi.sq.value, df)

        # I actually would like to know when 'NA' arises.
        if (is.na(result)) {
            return(NA)
        }
        else {
           return(result)
        }
    }
}

library(ggplot2)

# df: dataset (rows: cells, columns: genes)
# batch: batch labels for each cell
# k0 = NULL: # neighbors to test on
# knn = NULL: n_cells x k0 matrix of kNN indices (optional)
# testSize = NULL: # data points to test (default is 10% of entire dataset, but at least 25)
# do.pca = TRUE: perform PCA before kNN search
# dim.pca = 50: if do.pca=TRUE then # dims
# heuristic = TRUE: compute an optimal neighborhood size k
# n_repeat = 100: create statistics on batch estimates running test on n_repeat subsets
# alpha = 0.05: significance level
# addTest = FALSE: perform Likelihood Ratio Test approximation to the multinimial test and
#                a multinomial exact test (if appropriate)
# verbose = FALSE: display stages of current computation
# plot = TRUE: if stats > 10 then a boxplot of the resulting rejection rates is created
# adapt = TRUE: in some cases a number of cells do not contribute to any neighborhood and 
#              this may cause an imbalance in observed and expected batch label frequencies.
#              Frequencies will be adapted if adapt=TRUE.
# scenario = NULL: name of the scenario
# sc_dir = NULL: directory to save plots

kBET <- function(
    df, batch, k0, knn, testSize = NULL, do.pca = FALSE, dim.pca = 50,
    heuristic = FALSE, n_repeat = 100, alpha = 0.05, addTest = FALSE,
    verbose = TRUE, plot = TRUE, adapt = FALSE, scenario = NULL, sc_dir = NULL
) {
    dof <- length(unique(batch)) - 1    # degrees of freedom

    if (is.factor(batch)) {
        batch <- droplevels(batch)
    }

    frequencies <- table(batch) / length(batch)

    # get 3 different permutations of the batch label
    batch.shuff <- replicate(3, batch[sample.int(length(batch))])

    class.frequency <- data.frame(
        class = names(frequencies),
        freq = as.numeric(frequencies)
    )

    dataset <- df
    dim.dataset <- dim(dataset)

    # check dimensions
    if (dim.dataset[1] != length(batch) && dim.dataset[2] != length(batch)) {
        stop("Input matrix and batch information do not match. Execution halted.")
    }

    if (dim.dataset[2] == length(batch) && dim.dataset[1] != length(batch)) {
        if (verbose) {
            cat('Input matrix has samples as columns. kBET needs samples as rows. Transposing...\n')
        }
        dataset <- t(dataset)
        dim.dataset <- dim(dataset)
    }

    # check if the dataset is too small
    if (dim.dataset[1] <= 10) {
        if (verbose) {
            cat("Your dataset has less than 10 samples. Abort and return NA.\n")
        }
        return(NA)
    }

    stopifnot(class(n_repeat) == 'numeric', n_repeat > 0)

    # if k0 was set by the user and is too small & we do not operate on a knn graph, abort

    # the reason is that if we want to test kBET on knn graph data integration methods,
    # we usually face small numbers of nearest neighbours.
    if (k0 < 10 & is.null(knn))
    {
        if (verbose)
        {
            warning(
                "Your dataset has too few samples to run a heuristic.\n",
                "Return NA.\n",
                "Please assign k0 and set heuristic=FALSE."
            )
        }
        return(NA)
    }

    # set # tests
    if (is.null(testSize) || (floor(testSize) < 1 || dim.dataset[1] < testSize))
    {
        test.frac <- 0.10   # 10% of datapoints
        testSize <- ceiling(dim.dataset[1] * test.frac)
        
        if (testSize < 2 && dim.dataset[1] > 25)
        {
            testSize <- 25
        }

        if (verbose)
        {
            cat('Number of kBET tests is set to ')
            cat(paste0(testSize, '.\n'))
        }
    }

    # result list
    rejection <- list()
    rejection$summary <- data.frame(
        kBET.expected = numeric(4),
        kBET.observed = numeric(4),
        kBET.signif = numeric(4)
    )

    rejection$results <- data.frame(
        tested = numeric(dim.dataset[1]),
        kBET.pvalue.test = rep(0, dim.dataset[1]),
        kBET.pvalue.null = rep(0, dim.dataset[1])
    )

    # concat knn matrix so
    ## 1st col comes 1st then 2nd col ... then k0-1 th col then cells itself
    env <- as.vector(cbind(knn[, seq_len(k0 - 1)], seq_len(dim.dataset[1])))

    cf <- if (adapt && is.imbalanced) new.class.frequency else class.frequency

    # get avg p-value, p = 1 - CDF
    rejection$average.pval <- 1 - pchisq(k0 * residual_score_batch(env, cf, batch), dof)

    # initialize intermediates
    kBET.expected <- numeric(n_repeat)
    kBET.observed <- numeric(n_repeat)
    kBET.signif <- numeric(n_repeat)

    # kBET; run n_repeat = 100 times
    for (i in seq_len(n_repeat))
    {
        # choose a random sample from dataset (rows: cells, cols: genes)
        # size of random sample is 10% of dataset by default
        idx.runs <- sample.int(dim.dataset[1], size=testSize)
        
        # get the neighbors of sample cells; also attach sample cells as self-neighbor
        env <- cbind(knn[idx.runs, seq_len(k0 - 1)], idx.runs)
        # print(env)

        # perform test for each cell in sample
        cf <- if (adapt && is.imbalanced) new.class.frequency else class.frequency  # global batch props
        # apply chi_batch_test to each row of env matrix; 1 is rows
        p.val.test <- apply(env, 1, chi_batch_test, cf, batch, dof)
        # print(p.val.test)

        is.rejected <- p.val.test < alpha

        # apply the same test but this time batch labels are permuted. 
        # permutations are columns of batch.shuff; 2 is columns
        p.val.test.null <- apply(
            batch.shuff, 2,
            function(x) apply(env, 1, chi_batch_test, class.frequency, x, dof)
        )

        # expected rejection rate under the null for this run
        # for each shuffle compute the fraction of cells with p < alpha, then average across shuffles
        kBET.expected[i] <- mean(
            apply(
                p.val.test.null, 2,
                function(x) sum(x < alpha, na.rm = TRUE) / sum(!is.na(x))
            )
        )

        # observed rejection rate in this run
        # fraction of tested cells rejecting the null with the real labels
        kBET.observed[i] <- sum(is.rejected, na.rm = TRUE) / sum(!is.na(p.val.test))

        # compute significance
        kBET.signif[i] <- 1 - ptnorm(
            kBET.observed[i],
            mu = kBET.expected[i],
            sd = sqrt(kBET.expected[i] * (1 - kBET.expected[i]) / testSize),
            alpha = alpha
        )

        # mark which cells were sampled in this run and store
        ## their observed per-cell p-values
        ## their mean null per-cell p-values (averaged over shuffles)
        rejection$results$tested[idx.runs] <- 1
        rejection$results$kBET.pvalue.test[idx.runs] <- p.val.test
        rejection$results$kBET.pvalue.null[idx.runs] <- rowMeans(p.val.test.null, na.rm=TRUE)
    }

    # plot
    pvals <- rejection$results$kBET.pvalue.test
    tested <- rejection$results$tested

    df <- data.frame(
        cell = seq_along(pvals),
        pvalue = pvals,
        tested = factor(tested)
    )

    p <- ggplot(df, aes(x = cell, y = pvalue)) +
        # points for cells
        geom_point(aes(color = tested), alpha = 0.7) +
        # horizontal line as separate legend item
        geom_hline(aes(yintercept = alpha, color = "Significance level"), linetype = "dashed") +
        scale_color_manual(
            values = c(
            "0" = "lightgray",
            "1" = "steelblue",
            "Significance level" = "black"
            ),
            labels = c(
            "0" = "Not tested",
            "1" = "Tested",
            "Significance level" = "Significance level (0.05)"
            )
        ) +
        labs(
            title = paste0(
            "kBET p-values per Cell (", scenario, "; kBET=", round(mean(kBET.observed), 2), ")"),
            x = "Cell Index",
            y = "p-value",
            color = NULL
        ) +
        theme_minimal() +
        theme(
            legend.position = "top",
            legend.title = element_blank()
        )
    
    p <- p + theme(
        plot.title = element_text(size = 9)  # size in points
    )

    # save it
    ggsave(file.path(sc_dir, "kbet_pvalues_plot.pdf"), plot = p, width = 8, height = 5)
    ggsave(file.path(sc_dir, "kbet_pvalues_plot.png"), plot = p, width = 8, height = 5, dpi = 300, bg = "white")

    # summarize exptected, observed and run-level significance across runs
    if (n_repeat > 1)
    {
        CI95 <- c(0.025, 0.5, 0.975)

        rejection$summary$kBET.expected <- c(
            mean(kBET.expected, na.rm=TRUE),
            quantile(kBET.expected, CI95, na.rm=TRUE)
        )

        rownames(rejection$summary) <- c('mean', '2.5%', '50%', '97.5%')

        rejection$summary$kBET.observed <- c(
            mean(kBET.observed, na.rm=TRUE),
            quantile(kBET.observed, CI95, na.rm=TRUE)
        )

        rejection$summary$kBET.signif <- c(
            mean(kBET.signif, na.rm=TRUE),
            quantile(kBET.signif, CI95, na.rm=TRUE)
        )

        # also return per run vectors
        rejection$stats$kBET.expected <- kBET.expected
        rejection$stats$kBET.observed <- kBET.observed
        rejection$stats$kBET.signif <- kBET.signif

        if (n_repeat < 10)
        {
            cat('Warning: The quantile computation for ')
            cat(paste0(n_repeat))
            cat(' subset results is not meaningful.')
        }

        if (plot)
        {
            plot.data <- data.frame(
                class = rep(c('observed(kBET)', 'expected(random)'), each=n_repeat),
                data = c(kBET.observed, kBET.expected)
            )

            g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + theme_bw() + labs(title = scenario, x = 'Test', y = 'Rejection rate') + scale_y_continuous(limits = c(0, 1))
            # print(g)
            ggsave(file.path(sc_dir, "kBET_boxplot.pdf"), plot = g, width = 7, height = 5, dpi = 300)
        }
    } else
    {
        rejection$summary$kBET.expected <- kBET.expected[1]
        rejection$summary$kBET.observed <- kBET.observed[1]
        rejection$summary$kBET.signif <- kBET.signif[1]
    }

    # collect parameters
    rejection$params <- list()
    rejection$params$k0 <- k0
    rejection$params$testSize <- testSize
    rejection$params$do.pca <- do.pca
    rejection$params$dim.pca <- dim.pca
    rejection$params$heuristic <- heuristic
    rejection$params$n_repeat <- n_repeat
    rejection$params$alpha <- alpha
    rejection$params$addTest <- addTest
    rejection$params$verbose <- verbose
    rejection$params$plot <- plot

    return(rejection)
}