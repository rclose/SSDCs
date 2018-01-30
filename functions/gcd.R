#' Function For Calculating Great Circle Distances
#'
#' Calculate Great Circle Distance
#' @param arg Description
#' @keywords key 
#' @export
#' @examples
#' example()

gcd <- function(x,y,z) {

	if (x == y && z == 0) return(0) else 6371 * acos( (sin(x * pi / 180) * sin(y * pi / 180) ) + (cos(x * pi / 180) * cos(y * pi/180) * cos(z * pi / 180) ) )

}
