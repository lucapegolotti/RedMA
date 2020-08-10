// namespace RedMA
// {
//
// template <class InMatrixType>
// aMatrix<InMatrixType>::
// aMatrix() :
//   M_matrix(nullptr),
//   M_nRows(0),
//   M_nCols(0)
// {
// }
//
// template <class InMatrixType>
// void
// aMatrix<InMatrixType>::
// softCopy(std::shared_ptr<aMatrix> other)
// {
//     if (other)
//         M_matrix = other->M_matrix;
//
//     M_nRows = other->M_nRows;
//     M_nCols = other->M_nCols;
// }
//
// template <class InMatrixType>
// void
// aMatrix<InMatrixType>::
// hardCopy(std::shared_ptr<aMatrix> other)
// {
//     if (other)
//         M_matrix.reset(new InMatrixType(*other->M_matrix));
//
//     M_nRows = other->M_nRows;
//     M_nCols = other->M_nCols;
// }
//
//
// }
