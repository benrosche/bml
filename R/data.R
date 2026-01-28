#' Coalition governments in comparative perspective (1944-2013)
#'
#' @description
#' A dataset containing information about political parties participating in coalition
#' governments across multiple countries and time periods. The data are structured at
#' the party-government level (member level), making it suitable for multiple-membership
#' multilevel models where governments (groups) are composed of multiple parties (members).
#'
#' @format A data frame with 1,286 rows and 27 variables. Each row represents a party's
#'   participation in a specific government:
#' \describe{
#'   \item{pid}{Integer. Party identifier (member-level unit in mm() specification)}
#'   \item{pname}{Character. Party name}
#'   \item{gid}{Integer. Government identifier (group-level unit in mm() specification)}
#'   \item{gname}{Character. Government name}
#'   \item{cid}{Integer. Country identifier (nesting-level unit in hm() specification)}
#'   \item{cname}{Character. Country name}
#'   \item{gstart}{Date. Government start date}
#'   \item{gend}{Date. Government end date}
#'   \item{n}{Integer. Number of parties in the government (group size for weight functions)}
#'   \item{prime}{Integer. Prime minister party indicator}
#'   \item{pfam}{Character. Party family}
#'   \item{rile}{Numeric. Party's right-left ideological position}
#'   \item{ipd}{Numeric. Intra-party democracy score}
#'   \item{fdep}{Numeric. Party's financial dependency}
#'   \item{pseatrel}{Numeric. Party's relative seat share within coalition}
#'   \item{majority}{Integer. Majority government indicator: 1 = government controls majority
#'     of seats, 0 = minority government}
#'   \item{mwc}{Integer. Minimal winning coalition indicator}
#'   \item{hetero}{Numeric. Ideological heterogeneity within government}
#'   \item{investiture}{Numeric. Investiture vote requirement at country level}
#'   \item{pmpower}{Numeric. Prime ministerial powers at country level}
#'   \item{earlyterm}{Integer. Early termination indicator: 1 = government terminated early,
#'     0 = government completed full term (event indicator for survival models)}
#'   \item{govdur}{Numeric. Government duration in days (outcome variable for survival models)}
#'   \item{govmaxdur}{Numeric. Maximum possible government duration in days}
#'   \item{sim.w}{Numeric. Simulated weights for demonstration purposes}
#'   \item{sim.y}{Numeric. Simulated linear outcome}
#'   \item{sim.st}{Numeric. Simulated survival time}
#'   \item{sim.e}{Integer. Simulated event status}
#' }
#'
#' @details
#' This dataset is particularly useful for demonstrating multiple-membership models where:
#' \itemize{
#'   \item \strong{Members:} Political parties (identified by \code{pid})
#'   \item \strong{Groups:} Coalition governments (identified by \code{gid})
#'   \item \strong{Nesting:} Governments nested within countries (identified by \code{cid})
#' }
#'
#' Each government can have multiple parties, and parties can participate in multiple
#' governments over time. This creates a multiple-membership structure where party-level
#' characteristics need to be aggregated to the government level using appropriate
#' weighting functions.
#'
#' @source
#' Coalition government data compiled from multiple comparative politics sources covering
#' parliamentary democracies from 1944-2013.
#'
#' @seealso \code{\link{bml}} for modeling examples using this dataset
#'
#' @examples
#' data(coalgov)
#'
#' # Explore structure
#' head(coalgov)
#' table(coalgov$cname)  # Countries in sample
#'
#' # Summarize government durations
#' summary(coalgov$govdur[coalgov$earlyterm == 1])  # Early terminated
#' summary(coalgov$govdur[coalgov$earlyterm == 0])  # Full term
#'
#' # Party participation patterns
#' table(table(coalgov$pid))  # Distribution of government participations per party
#'
#' \dontrun{
#' # Basic model: government duration as function of majority status and party fragmentation
#' m1 <- bml(
#'   Surv(govdur, earlyterm) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE) +
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

#' School networks data from CILS4EU study
#'
#' @description
#' Network and survey data from the Children of Immigrants Longitudinal Survey in
#' Four European Countries (CILS4EU). The data include friendship networks among
#' students nested within classrooms and schools, along with individual-level
#' characteristics and outcomes. Structured for hierarchical and multiple-membership
#' network models.
#'
#' @format Two data objects are provided:
#'
#' \strong{nodedat}: Node-level (student) attributes
#' \describe{
#'   \item{youthid}{Integer. Student identifier (unique across all schools)}
#'   \item{schoolid}{Integer. School identifier}
#'   \item{classid}{Integer. Classroom identifier}
#'   \item{sex}{Factor. Student gender}
#'   \item{etn}{Factor. Ethnicity or immigrant background}
#'   \item{parent_edu}{Numeric. Parental education level}
#'   \item{parent_inc}{Numeric. Parental income}
#'   \item{test_lang}{Numeric. Language test score}
#'   \item{test_cogn}{Numeric. Cognitive test score}
#' }
#'
#' \strong{edgedat}: Edge-level (friendship) data
#' \describe{
#'   \item{youthid_from}{Integer. Nominator (student reporting friendship)}
#'   \item{youthid_to}{Integer. Nominee (student nominated as friend)}
#'   \item{rank}{Integer. Ranking of friendship nomination}
#' }
#'
#' @details
#' This dataset is designed for social network analysis with multilevel structure:
#' \itemize{
#'   \item \strong{Level 1:} Students (nodes)
#'   \item \strong{Level 2:} Classrooms
#'   \item \strong{Level 3:} Schools
#' }
#'
#' Students can form friendships both within and across classrooms (though primarily
#' within classrooms). This creates opportunities for modeling peer effects, social
#' influence, and network formation processes using multiple-membership or hierarchical
#' specifications.
#'
#' \strong{Potential Applications:}
#' \itemize{
#'   \item Peer influence on academic performance
#'   \item Homophily in friendship formation
#'   \item Cross-level interactions (school/classroom effects on peer networks)
#'   \item Multiple-membership models where students are influenced by multiple peers
#' }
#'
#' @source
#' CILS4EU - Children of Immigrants Longitudinal Survey in Four European Countries.
#' A cross-national longitudinal survey of children of immigrants and native-born
#' children in England, Germany, the Netherlands, and Sweden.
#'
#' @seealso \code{\link{bml}} for modeling approaches
#'
#' @examples
#' data(schoolnets)
#'
#' # Explore node data
#' head(nodedat)
#' table(nodedat$schoolid)  # Students per school
#'
#' # Explore network structure
#' head(edgedat)
#' nrow(edgedat)  # Total friendships
#'
#' \dontrun{
#' # Example: Model test scores as function of peer characteristics
#' # This would require preparing data with peer-level aggregation
#' }
#'
#' @rdname schoolnets
#' @docType data
#' @keywords datasets
"nodedat"

#' @rdname schoolnets
"edgedat"
