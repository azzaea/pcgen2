plot.nw <- function(object, title, fillcolor = "lightblue", sub = F, fontsize = 50){
  paroptions <- par()$mar
  if (class(object) == 'pcAlgo') {
    Rgraphviz::plot(object@graph, main = title,
                    attrs = list(node = list(fillcolor = fillcolor, fontsize = fontsize)))

    print(object@call)
    if (isTRUE(sub)) {
      title(sub = paste(strwrap(object@call, 150), collapse = "\n"),  outer = T,
            line = -2, adj = .8)
    }
  } else
    Rgraphviz::plot(object, main = title,
                    attrs = list(node = list(fillcolor = fillcolor, fontsize = fontsize)))
  on.exit(par(mar = paroptions))
}


eval.nw <- function(runs, refgraph, suffStat, K, cov.method = "uni", m.max = Inf,
                    method = "pc") {
  bench <- tibble::tibble(tpr = numeric(runs), fpr = numeric(runs),
                  tdr = numeric(runs), shd = numeric(runs))

  for (i in 1:runs){
    print(i)
    if (method == "pc"){
      nw.res <- pcRes(suffStat, K = K, cov.method = cov.method,
                        m.max = m.max)
      bench[i, c("tpr", "fpr", "tdr")] <- t(pcalg::compareGraphs(nw.res@graph, refgraph))
      bench[i, "shd"] <- pcalg::shd(nw.res@graph, refgraph)
    } else{
      nw.res <- fciRes(suffStat, K = K, cov.method = cov.method)
      bench[i, c("tpr", "fpr", "tdr")] <- t(pcalg::compareGraphs(as(nw.res@amat, "graphNEL"), refgraph))
      bench[i, "shd"] <- pcalg::shd(as(nw.res@amat, "graphNEL"), refgraph)
      }

  }

  print(paste(deparse(substitute(refgraph)), "vs pcRes(K = ",
              deparse(substitute(K)), ") run", runs, "times"))

  g1 <- bench %>% pivot_longer(-shd, names_to = "metric", values_to = "value") %>%
    ggplot() +
    geom_point(aes(x = metric, y = value)) + geom_boxplot(aes(x = metric, y = value)) + theme_bw()

  g2 <- bench %>% pivot_longer(shd, names_to = "metric", values_to = "value") %>%
               ggplot(aes(x = metric, y = value)) +
               geom_point() + geom_boxplot() + theme_bw() + theme(axis.title.y = element_blank())

  plot_grid( g1, g2 , rel_widths = c(2, 1))


}

#' Based on pcalg implementation upon ordering the nodes of the true network, and
#' making the names of these nodes identical (here, by replace space with '.')
compareGraphs.order <- function (gl, gt) {
  graph::nodes(gt) <- gsub(" ", ".", graph::nodes(gt))
  ml <- pcalg::wgtMatrix(graph::ugraph(gl))
  mt <- pcalg::wgtMatrix(graph::ugraph(gt))
  mt <- mt[rownames(ml), colnames(ml)]

  p <- dim(ml)[2]
  mt[mt != 0] <- rep(1, sum(mt != 0))
  ml[ml != 0] <- rep(1, sum(ml != 0))
  diffm <- ml - mt
  nmbTrueGaps <- (sum(mt == 0) - p)/2
  fpr <- ifelse (nmbTrueGaps == 0, 1,  (sum(diffm > 0)/2)/nmbTrueGaps)
  diffm2 <- mt - ml
  nmbTrueEdges <- (sum(mt == 1)/2)
  tpr <- ifelse (nmbTrueEdges == 0, 0, 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges)
  trueEstEdges <- (nmbTrueEdges - sum(diffm2 > 0)/2)
  tdr <- ifelse (sum(ml == 1) == 0,
                 ifelse (trueEstEdges == 0, 1, 0),
                 trueEstEdges/(sum(ml == 1)/2))

  c(tpr = tpr, fpr = fpr, tdr = tdr)
}

#' Similarly, based on pcalg implementation upon rearranging nodes' order in g1
shd.order <- function(g1, g2){
  if (is(g1, "pcAlgo"))
    g1 <- g1@graph
  if (is(g2, "pcAlgo"))
    g2 <- g2@graph
  if (is(g1, "graphNEL")) {
    m1 <- wgtMatrix(g1, transpose = FALSE)
    m1[m1 != 0] <- 1
  }
  graph::nodes(g2) <- gsub(" ", ".", graph::nodes(g2))
  if (is(g2, "graphNEL")) {
    m2 <- wgtMatrix(g2, transpose = FALSE)
    m2[m2 != 0] <- 1
  }
  m2 <- m2[rownames(m1), colnames(m1)]
  shd <- 0
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
  d <- abs(m1 - m2)
  shd + sum((d + t(d)) > 0)/2
}


plot.performance <- function(dscobject, method) {
  dscobject <- dscobject %>%
    group_by(simulate.dagcoeff, simulate.rangeGenVar, simulate.rangeEnvVar, learn) %>%
    summarise(across(learn.dagfit_tpr:learn.genfit_tdr, mean)) %>%
    filter(learn == method)

  nw.traits <-  dscobject %>%
    pivot_longer(cols = starts_with("learn.dagfit")) %>%
    mutate(name = str_replace(name, "learn.", "")) %>%
    mutate(name = fct_relevel(name, "dagfit_shd")) %>% ggplot() +
    geom_contour_filled(aes(x = simulate.rangeGenVar, y = simulate.rangeEnvVar,
                            z = value), bins = 10) +
    facet_grid(cols = vars(name), rows = vars(simulate.dagcoeff)) +
    labs(title = bquote('Learnt network via' ~ bold(.(method)))) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_classic()
  print(nw.traits)

  if (method == "pcgen"){
    nw.genotype <- dscobject %>%
      pivot_longer(cols = starts_with("learn.genfit")) %>%
      mutate(name = str_replace(name, "learn.", "")) %>% ggplot() +
      geom_contour_filled(aes(x = simulate.rangeGenVar, y = simulate.rangeEnvVar,
                              z = value), bins = 10) +
      facet_grid(cols = vars(name), rows = vars(simulate.dagcoeff)) +
      labs(title = bquote('Learnt genetic network via' ~ bold(.(method)))) +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(trans = "log10") +
      theme_classic()
    print(nw.genotype)
  }
}




