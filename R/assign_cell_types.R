#' assign_cell_types
#' 
#' The function chooses the cell type for each cell.
#' @param ttype tissue types of the cells
#' @param nttype total number of possible tissue types. should be no less than the max value of ttype and the 
#' number of unique values in ttype. some tissue types may not exist in ttype
#' @param cell_type_proportion matrix of cell type proportion, cell type (row) by tissue type (column)
#' @return A vector of the assigned cell type
#' @details This function choose the cell type for each cell. Each tissue type could contains
#' multiple types of cells. Based on the cell type proportion, this function randomly samples 
#' the cell type for each cell based on their corresponding tissue type.
#' @export
#' 
assign_cell_types <- function(ttype, nttype, cell_type_proportion = NULL){
  if (is.null(cell_type_proportion)){
    cell_type_proportion <- matrix(0, nttype, nttype)
    diag(cell_type_proportion) <- 0.8
    cell_type_proportion[cell_type_proportion == 0] <- 0.2 / (nttype - 1)
  }
  if(ncol(cell_type_proportion) != nttype){
    stop("Number of columns of cell_type_proportion matrix does not equal to nttype.")
  }
  true_cell_type <- sapply(ttype, function(x){which(rmultinom(1, 1, prob = cell_type_proportion[, x]) == 1)})
  return(true_cell_type)
}