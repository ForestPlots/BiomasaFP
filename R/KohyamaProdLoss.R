#'Calculating biomass production and loss from repeated censuses
#'
#' Code provided in SI of Koyhama et al. 2019 for estimating biomass gains and losses. Used by SummaryAGWP_Kohyama wrapper.
#' @param x1 Numeric vector. Biomass of stems in first census (Mg)
#' @param x2 Numeric vector. Biomass of stems in second census (Mg)
#' @param xs1 Numeric vector. Initial biomass of stems that survived from first census to second census (Mg)
#' @param intvl Numeric. Census interval length (years).
#' @param area Numeric. Plot area (ha)
#' @param subpop Vector giving subpopulation membership of each stem
#' @param simple Logical. Only used to obtain the different outputs from Koyhama et al.
#' @param annual Logical. Only used to obtain the different outputs from Koyhama et al.
#' @references Kohyama et al. 2019 Estimating net biomass production and loss from repeated measurements of trees in forests and woodlands: Formulae, biases and recomendations. Forest Ecology and Management 433: 729-740.
#' @author T.S. Kohyama, T.I. Kohyama and D. Sheil

turnover = function (x1, x2, xs1, intvl, area,
                 subpop = NULL, simple = FALSE, annual = FALSE) {
    if (length(subpop) == 0) subpop = rep(1, length(x1))
    B0k = tapply(x1, subpop, sum)
    BTk = tapply(x2, subpop, sum)
    Bs0k = tapply(xs1, subpop, sum)
    Bwk = ifelse(BTk != B0k, (BTk-B0k)/log(BTk/B0k), B0k)  # Eq (10)
    Bw_annk = ifelse(BTk != B0k, (BTk-B0k)/((BTk/B0k)^(1/intvl) - 1)/intvl, B0k)  # Eq (14)

    B0k = B0k/area
    BTk = BTk/area
    Bs0k = Bs0k/area
    Bwk = Bwk/area
    Bw_annk = Bw_annk/area

    L_simple  = sum(B0k - Bs0k)/intvl
    P_simple  = sum(BTk - Bs0k)/intvl
    Bw_simple = sum(B0k + BTk)/2

    L_ann  = sum(Bw_annk * (1 - (Bs0k/B0k)^(1/intvl)))
    P_ann  = sum(Bw_annk * (BTk/B0k)^(1/intvl) * (1 - (Bs0k/BTk)^(1/intvl)))
    Bw_ann = sum(Bw_annk)

    L = sum(Bwk * log(B0k/Bs0k))/intvl
    P = sum(Bwk * log(BTk/Bs0k))/intvl
    Bw = sum(Bwk)

    if (simple == TRUE) {
        return(list('P' = P_simple, 'L' = L_simple, 'Bw' = Bw_simple))
    } else if (annual == TRUE) {
        return(list('P' = P_ann, 'L' = L_ann, 'Bw' = Bw_ann))
    } else {
        return(list('P' = P, 'L' = L, 'Bw' = Bw))
    }
}
