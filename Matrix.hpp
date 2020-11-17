#include <stdlib.h>
#include <iostream>
#include <vector>
#include "helper_funcs.cpp"

#define STORAGE_TYPE long double
#define OVERFLOW_PROTECTION 13

class Matrix {
private:
    bool SwapRowsWith0Pivot(); // self explanatory name
    Matrix& swap_row(u_int32_t, u_int32_t); // Swap row[ind1] with row[ind2]
    Matrix& swap_col(u_int32_t, u_int32_t); // Swap col[ind1] with col[ind2]
    Matrix& row_mult(u_int32_t, STORAGE_TYPE); // multiply data[ind] by scale
    Matrix& row_divide(u_int32_t, STORAGE_TYPE); // divide data[ind] by scale
    Matrix& row_sub(u_int32_t, u_int32_t); // substract data[ind2] from data[ind1] 
public:
    u_int32_t rows;
    u_int32_t cols;
    std::vector<std::vector<STORAGE_TYPE>> data; // 2d array
    
    Matrix (u_int32_t, u_int32_t); // Initialize the data array
    Matrix* copy(); // copy this matrix into a new one
    void random(u_int32_t); // randomize data from 1 to the value passed into Matrix.random
    Matrix& operator= (Matrix&); // assignment operator implementation
    static Matrix* from_array (STORAGE_TYPE*[], u_int32_t, u_int32_t); // initialize the data array from an array 
    ~Matrix (); // Free memory holding data array
    Matrix& operator+ (Matrix&); // add a matrix to this one
    Matrix& operator- (Matrix&); // subtract a matrix from this one
    Matrix& operator* (STORAGE_TYPE); // scale the matrix
    Matrix& operator* (Matrix&); // dot product stored in new matrix
    Matrix& operator* (Matrix*); // dot product stored in new matrix
    static Matrix* transpose (Matrix); // transpose a matrix into a new one
    std::vector<STORAGE_TYPE>& operator[] (size_t); // return the row at a given index
    void print(); // display the contents of data in stdout
    STORAGE_TYPE determinant(); // Take the determinant of this using the Bareiss algorithm
    Matrix& inverse(); // Take the matrix inverse of this
    Matrix& to_identity(); // turn this.data into an identity matrix
}; // end of class

// Initialize the data array
Matrix::Matrix (u_int32_t rows_, u_int32_t cols_) 
: rows(rows_), cols(cols_)
{
    for (int row = 0; row < rows; row++) 
    {
        std::vector<STORAGE_TYPE> row2add;
        for (int col = 0; col < cols; col++) 
        {
            row2add.push_back(0.0);
        }
        data.push_back(row2add);
    }
}

// copy this matrix into a new one
Matrix* Matrix::copy()
{
    Matrix* mat = new Matrix(rows, cols);
    for (int row = 0; row < rows; row++) 
    {
        for (int col = 0; col < cols; col++) 
        {
            mat->data[row][col] = data[row][col];
        }
    }
    return mat;
}

// randomize data from 1 to the value passed into Matrix.random
void Matrix::random(u_int32_t max) 
{
    for (int row = 0; row < rows; row++) 
    {
        for (int col = 0; col < cols; col++) 
        {
            data[row][col] = rand() % max + 1;  
        }
    }
}

// assignment operator implementation
Matrix& Matrix::operator= (Matrix& mat)
{
    if (this == &mat) 
        return *this;

    rows = mat.rows;
    cols = mat.cols;
    data = mat.data;
    return *this;
}

// initialize the data array from an array
Matrix* Matrix::from_array (STORAGE_TYPE* arr[], u_int32_t rows, u_int32_t cols) 
{
    Matrix* mat = new Matrix(rows, cols);
    for (int row = 0; row < mat->rows; row++)
    {
        for (int col = 0; col < mat->cols; col++)
        {
            (*mat)[row][col] = arr[row][col];
        }
    }
    return mat;
}  

// Free memory holding data array
Matrix::~Matrix ()
{
    // Free the data array
    // JK! got ya there XD
    // Vectors handle themselves and at the moment
    // there are no pointers to free in this class
}

// add a matrix to this one
Matrix& Matrix::operator+ (Matrix& mat) 
{
    // Catch any errors
    if (! (mat.cols == cols && mat.rows == rows) )
        throw "Different sized matrices cannot be added";
    
    for (int row = 0; row < rows; row++) 
    {
        for (int col = 0; col < cols; col++) 
        {
            data[row][col] += mat[row][col];
        }
    }

    return *this; // allows chaining of operations
}

 // subtract a matrix from this one
Matrix& Matrix::operator- (Matrix& mat)
{
    // Catch any errors
    if (! (mat.cols == cols && mat.rows == rows) )
        throw "Different sized matrices cannot be added";
    
    for (int row = 0; row < rows; row++) 
    {
        for (int col = 0; col < cols; col++) 
        {
            data[row][col] -= mat[row][col];
        }
    }

    return *this; // allows chaining of operations
}

// scale the matrix
Matrix& Matrix::operator* (STORAGE_TYPE scalar)
{    
    for (int row = 0; row < rows; row++) 
    {
        for (int col = 0; col < cols; col++) 
        {
            data[row][col] *= scalar;
        }
    }
    // Correct any computer bit errors (bc working with big/small numbers)
    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++) 
        {
            if ((*this).data[row][col] <= (1.0/pow(10, OVERFLOW_PROTECTION)))
            {
                (*this).data[row][col] = 0;
            }
        }
    }

    return *this; // allows chaining of operations
} 

// dot product stored in new matrix
Matrix& Matrix::operator* (Matrix& mat)
{
    if (! (rows == mat.cols && cols == mat.rows) )
        throw "Improper sizes for dot product";
    
    Matrix* out = new Matrix(rows, mat.cols);
    for (int out_row = 0; out_row < out->rows; out_row++) 
    {
        for (int out_col = 0; out_col < out->cols; out_col++) 
        {
            STORAGE_TYPE sum = 0;
            for (int col = 0; col < cols; col++)
            {
                sum += data[out_row][col] * mat[col][out_col];
            }
            (*out)[out_row][out_col] = sum;
        }
    }

    // Correct any computer bit errors (bc working with big/small numbers)
    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++) 
        {
            if ((*out).data[row][col] <= (1.0/pow(10, OVERFLOW_PROTECTION)))
            {
                (*out).data[row][col] = 0;
            }
        }
    }

    return *out;
} 

// dot product stored in new matrix
Matrix& Matrix::operator* (Matrix* mat)
{
    if (! (rows == mat->cols && cols == mat->rows) )
        throw "Improper sizes for dot product";
    
    Matrix* out = new Matrix(rows, mat->cols);
    for (int out_row = 0; out_row < out->rows; out_row++) 
    {
        for (int out_col = 0; out_col < out->cols; out_col++) 
        {
            STORAGE_TYPE sum = 0;
            for (int col = 0; col < cols; col++)
            {
                sum += data[out_row][col] * (*mat)[col][out_col];
            }
            (*out)[out_row][out_col] = sum;
        }
    }
    // Correct any computer bit errors (bc working with big/small numbers)
    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++) 
        {
            if ((*out).data[row][col] <= (1.0/pow(10, OVERFLOW_PROTECTION)))
            {
                (*out).data[row][col] = 0;
            }
        }
    }

    return *out;
} 

// transpose a matrix into a new one
Matrix* Matrix::transpose (Matrix mat)
{
    Matrix* out = new Matrix(mat.cols, mat.rows);

    for (int row = 0; row < out->rows; row++) 
    {
        for (int col = 0; col < out->cols; col++) 
        {
            (*out)[row][col] = mat[col][row];
        }
    }

    return out;
} 

// return the row at a given index
std::vector<STORAGE_TYPE>&  Matrix::operator[] (size_t row)
{
    if (row >= rows) 
        throw "Index out of bounds";

    return data[row];
} 

// display the contents of data in stdout
void Matrix::print() 
{
    for (std::vector<STORAGE_TYPE> row : data) 
    {
        for (STORAGE_TYPE value : row)
        {
            std::cout << value << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Swap row[ind1] with row[ind2]
Matrix& Matrix::swap_row(u_int32_t ind1, u_int32_t ind2)
{
    std::vector<STORAGE_TYPE> temp = data[ind1];
    data[ind1] = data[ind2];
    data[ind2] = temp;

    (*this) = (*this) * -1;

    return *this;
} 

// Swap col[ind1] with col[ind2]
Matrix& Matrix::swap_col(u_int32_t ind1, u_int32_t ind2)
{
    for (int row = 0; row < rows; row++) 
    {
        STORAGE_TYPE temp = data[row][ind1];
        data[row][ind1] = data[row][ind2];
        data[row][ind2] = temp;
    }

    (*this) = (*this) * -1;

    return *this;
} 

bool Matrix::SwapRowsWith0Pivot()
{
    bool singular = false;

    int row = -1;
    
    while ( row++ < rows and !singular )
    {
        if ( data[row][row] == 0.0 ) // check to see if matrix is singular
        {
            int row2swap = 0;
            while ( row2swap++ < rows && data[row2swap][row] == 0.0 ) {}; // find which row to swap
            if (row2swap < rows)
            {
                swap_row(row, row2swap);
            } else
            {
                singular = true;
            }
            
        }
    }
    
    return singular;
}

// Take the determinant of this using the Bareiss algorithm
STORAGE_TYPE Matrix::determinant()
{
    if (SwapRowsWith0Pivot()) // if we are singular
        return 0; // determinant is 0
    
    Matrix copy = *(*this).copy();

    STORAGE_TYPE pivot = 1.0;
    // Traverse pivots 
    for (int k = 0; k < rows-1; k++)
    {
        // Traverse rows
        for (int row = k + 1; row < rows; row++)
        {// Traverse columns
            for (int col = k + 1; col < cols; col++)
            {
                copy[row][col] = copy[k][k] * copy[row][col] - copy[row][k] * copy[k][col];
                copy[row][col] = copy[row][col] / pivot;
            }
        }
        pivot = copy[k][k];
    }

    return copy[rows-1][cols-1];
}

// Take the matrix inverse of this
Matrix& Matrix::inverse()
{
    Matrix* backup = (*this).copy();
    Matrix* out = new Matrix(rows, cols);
    out->to_identity();

    for (int col = 0; col < cols; col++)
    {
        int target_row_ind = col;
        STORAGE_TYPE scale = data[target_row_ind][col];
        (*this).row_divide(target_row_ind, scale);
        (*out).row_divide(target_row_ind, scale);
        for (int row = 0; row < rows; row++) 
        {
            STORAGE_TYPE target = (STORAGE_TYPE)( col == row );
            if (target == 1)
            {
                // already did this above (ie outer loop)
            } else if (target == 0)
            {
                STORAGE_TYPE scale = data[row][col];
                (*this).row_mult(target_row_ind, scale);
                (*this).row_sub(row, target_row_ind);
                (*this).row_divide(target_row_ind, scale);
                (*out).row_mult(target_row_ind, scale);
                (*out).row_sub(row, target_row_ind);
                (*out).row_divide(target_row_ind, scale);

            }
        }
    }

    // restore the original data;
    (*this).data = (*backup).data;

    return *out;
} 

// multiply data[ind] by scale
Matrix& Matrix::row_mult(u_int32_t ind, STORAGE_TYPE scale)
{
    for (int col = 0; col < cols; col++)
    {
        data[ind][col] *= scale;
    }
    return *this;
} 

// divide data[ind] by scale
Matrix& Matrix::row_divide(u_int32_t ind, STORAGE_TYPE scale)
{
    for (int col = 0; col < cols; col++)
    {
        data[ind][col] /= scale;
    }
    return *this;
}

// substract data[ind2] from data[ind1] 
Matrix& Matrix::row_sub(u_int32_t ind1, u_int32_t ind2)
{
    for (int col = 0; col < cols; col++)
    {
        data[ind1][col] -= data[ind2][col];
    }
    return *this;
}

// turn this.data into an identity matrix
Matrix& Matrix::to_identity()
{
    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            data[row][col] = (STORAGE_TYPE)(row == col);
        }
    }
    return *this;
}