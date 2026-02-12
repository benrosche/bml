#' Coalition Governments in Western Democracies (1944-2014)
#'
#' @description
#' A dataset containing information on coalition governments and their member parties
#' across 30 parliamentary democracies. The data are in long format where the unit of
#' analysis is parties in governments, making it suitable for multiple-membership
#' multilevel models where governments (groups) are composed of multiple parties (members).
#'
#' @format A tibble with 2,077 rows and 18 variables. Each row represents a party's
#'   participation in a specific coalition government. The sample contains 628 governments
#'   formed by 312 unique parties across 29 countries.
#'
#' \strong{Identifiers:}
#' \describe{
#'   \item{gid}{Government identifier (group-level unit in \code{mm()} specification).
#'     Range: [3, 1105]}
#'   \item{pid}{Party identifier (member-level unit in \code{mm()} specification).
#'     Range: [11110, 96955]}
#'   \item{cid}{Country identifier (nesting-level unit in \code{hm()} specification).
#'     Range: [11, 96]}
#'   \item{cname}{Three-letter country code (ISO 3166-1 alpha-3)}
#'   \item{pname}{Full party name}
#' }
#'
#' \strong{Government-level variables:}
#' \describe{
#'   \item{election}{Date of the preceding election that led to the government's formation.
#'     Range: [1939-04-02, 2014-12-14]}
#'   \item{n}{Number of parties in the coalition (group size for weight functions).
#'     Range: [2, 9], mean: 3.31}
#'   \item{dur_wkb}{Government duration in days, measured from investiture to termination
#'     (outcome variable for survival models). Range: [7, 1840], mean: 554.5}
#'   \item{event_wkb}{Early termination indicator: 1 = government terminated due to
#'     political conflict (voluntary resignation, dissension within government, lack of
#'     parliamentary support, or head of state intervention) more than one year before
#'     the official end of term; 0 = censored (regular elections, other reasons, or
#'     termination within one year of scheduled elections). Range: [0, 1], mean: 0.39}
#'   \item{majority}{Majority government indicator: 1 = coalition controls majority of
#'     parliamentary seats, 0 = minority government. Range: [0, 1], mean: 0.80}
#'   \item{mwc}{Minimal winning coalition indicator: 1 = coalition would lose its majority
#'     if any party left, 0 = oversized coalition. Range: [0, 1], mean: 0.35}
#'   \item{rile_SD}{Inter-party ideological heterogeneity. Standard deviation of coalition
#'     parties' left-right positions (from CMP) relative to the ideological distribution
#'     of all parties in parliament. Standardized and inverted so higher values indicate
#'     greater ideological cohesion. Range: [-8.40, 2.12], mean: 0.04}
#' }
#'
#' \strong{Country-level variables:}
#' \describe{
#'   \item{investiture}{Investiture vote requirement (time-constant country characteristic):
#'     1 = country requires formal parliamentary investiture vote, 0 = no formal requirement.
#'     Range: [0, 1], mean: 0.46}
#' }
#'
#' \strong{Party-level variables:}
#' \describe{
#'   \item{pseat}{Party's relative seat share within the coalition, computed as
#'     \code{pseat / sum(pseat)} within each government. Sums to 1 within each coalition.
#'     Range: [0.00, 1.00], mean: 0.33}
#'   \item{prime}{Prime minister party indicator: \code{TRUE} = party holds prime
#'     ministership (n = 628), \code{FALSE} = junior coalition partner (n = 1,449)}
#'   \item{cohesion}{Intra-party ideological cohesion, measured using an adaptation of
#'     the Cowles-Jones ratio. Computed as the ratio of continuous ideological shifts
#'     to reversals in a party's left-right position over time. Higher values indicate
#'     more consistent ideological trajectories (greater cohesion). Standardized.
#'     Range: [-1.13, 3.85], mean: 0.00}
#'   \item{rile}{Party's left-right ideological position (from CMP). Measured on a
#'     continuous scale where higher values indicate more right-wing positions and lower
#'     values indicate more left-wing positions. Standardized.
#'     Range: [-3.21, 3.68], mean: 0.00}
#'   \item{finance}{Party's economic dependence on member contributions (from PPDB).
#'     Measured as the share of party funding from member dues relative to total income.
#'     Standardized; higher values indicate greater dependence on member financing.
#'     Treated as time-constant due to data limitations. Range: [-0.98, 4.40], mean: 0.00}
#'   \item{Nmembers}{Number of party members (from PPDB). Standardized; treated as
#'     time-constant due to data limitations. Range: [-0.33, 15.02], mean: 0.00}
#' }
#'
#' @details
#' This dataset demonstrates multiple-membership multilevel modeling where:
#' \itemize{
#'   \item \strong{Members:} Political parties (identified by \code{pid})
#'   \item \strong{Groups:} Coalition governments (identified by \code{gid})
#'   \item \strong{Nesting:} Governments nested within countries (identified by \code{cid})
#' }
#'
#' Each coalition government comprises multiple parties, and parties can participate in
#' multiple governments over time. This creates a multiple-membership structure where
#' party-level characteristics are aggregated to the government level using weighting
#' functions specified in \code{mm()} blocks.
#'
#' \strong{Sample:} After matching party data across sources and excluding single-party
#' and caretaker governments, the sample comprises 628 governments formed by 312 unique
#' parties across 29 countries: Australia, Austria, Belgium, Bulgaria, Croatia, Czech
#' Republic, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Israel,
#' Italy, Japan, Latvia, Lithuania, Netherlands, Norway, Poland, Portugal, Romania,
#' Slovakia, Spain, Sweden, Switzerland, and United Kingdom.
#'
#' \strong{Measurement notes:}
#' \itemize{
#'   \item Government duration follows the WKB convention: time from investiture to
#'     termination or new elections
#'   \item Early termination events focus on political gridlock (conflict-related endings)
#'     and exclude terminations within one year of scheduled elections
#'   \item Party-level variables (\code{cohesion}, \code{finance}, \code{Nmembers}) are
#'     standardized (mean = 0) for analysis
#' }
#'
#' @source
#' Data compiled from multiple sources:
#' \itemize{
#'   \item \strong{Coalition governments:} Woldendorp, Keman, and Budge (WKB) dataset,
#'     updated by Seki and Williams (2014)
#'   \item \strong{Party ideology:} Comparative Manifesto Project (CMP; Volkens et al. 2016)
#'   \item \strong{Party organization:} Political Party Database (PPDB; Scarrow, Poguntke,
#'     and Webb 2017)
#' }
#'
#' Missing party-level data imputed using multiple imputation by chained equations with
#' predictive mean matching.
#'
#' @references
#' Seki, K., & Williams, L. K. (2014). Updating the Party Government data set.
#' \emph{Electoral Studies}, 34, 270-279.
#'
#' Volkens, A., et al. (2016). The Manifesto Data Collection. Manifesto Project
#' (MRG/CMP/MARPOR). Version 2016a. Berlin: Wissenschaftszentrum Berlin fur Sozialforschung.
#'
#' Scarrow, S. E., Webb, P. D., & Poguntke, T. (Eds.). (2017). \emph{Organizing Political
#' Parties: Representation, Participation, and Power}. Oxford University Press.
#'
#' @seealso \code{\link{bml}} for modeling examples using this dataset
#'
#' @examples
#' data(coalgov)
#'
#' # Explore data structure
#' str(coalgov)
#' table(coalgov$cname)
#'
#' # Number of unique units
#' length(unique(coalgov$gid))   # Governments
#' length(unique(coalgov$pid))   # Parties
#' length(unique(coalgov$cid))   # Countries
#'
#' \dontrun{
#' # Model: government duration as function of majority status and party characteristics
#' m1 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = TRUE) +
#'     hm(id = id(cid), type = "RE"),
#'   family = "Weibull",
#'   data = coalgov
#' )
#' summary(m1)
#' }
#'
#' @docType data
#' @keywords datasets
"coalgov"
