#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using std::vector;
using std::string;
using std::cin;
using std::cout;
using std::endl;

template <typename T>
using matrix = vector<vector<T>>;

namespace io_util {
  template <typename Uint>
  matrix<Uint> read_matrix(std::ifstream& in) {
    matrix<Uint> mat;
    vector<Uint> row;
    
    while (in.good()) {
      char c = in.get();

      switch (c) {
      case '0':
        row.push_back(0);
        break;
      case '1':
        row.push_back(1);
        break;
      case '\n':
        mat.push_back(row);
        row.clear();
        break;
      default:
        break;
      }
    }

    return mat;
  }

  template<typename T>
  void operator<<(std::ostream& out, const matrix<T>& mat) {
    out << mat.size() << "x" << mat.at(0).size() << " matrix:" << endl;
    
    for (const auto& row : mat) {
      for (const T& e : row) {
        out << e << " ";
      }

      out << endl;
    }
  }
}

/*

Sobre el tipo Uint:
- debe ser un tipo de número entero sin signo (unsigned short, unsigned int, unsigned long).
- determina el tamaño máximo del cuadrado maximal: por ejemplo, usar unsigned short será válido si el lado del cuadrado maximal mide 65535 o menos.

*/
template<typename Uint>
class LocalSolution: public matrix<Uint> {
  bool locally_solved;
  Uint max_square_size;
  
  inline Uint triangle(Uint left, Uint right, Uint up, Uint down) {
    return (std::min(std::min(left, right), up) + 1) * down;
  }
  
public:
  LocalSolution(matrix<Uint> m) :matrix<Uint>(m), locally_solved(false), max_square_size(0) {}
  
  Uint maximal_square() {
    if (!locally_solved) {
      Uint max = 0;
      size_t max_col = (*this)[0].size();
    
      for (size_t i = 0; i < this->size(); i++) {
        if (max < (*this)[i][0])
          max = (*this)[i][0];
      }

      for (size_t j = 0; j < max_col; j++) {
        if (max < (*this)[0][j])
          max = (*this)[0][j];
      }
    
      for (size_t i = 1; i < this->size(); i++) {
        for (size_t j = 1; j < max_col; j++) {
          (*this)[i][j] = triangle((*this)[i - 1][j], (*this)[i][j - 1], (*this)[i-1][j - 1], (*this)[i][j]);
        
          if (max < (*this)[i][j])
            max = (*this)[i][j];
        }
      }

      locally_solved = true;
      max_square_size = max * max;
    }
    
    return max_square_size;
  }

  /*
  Uint maximal_square(size_t row_start = 0, size_t row_end = size()) {
    Uint max = 0;
    size_t max_col = (*this)[row_start].size();
    
    for (size_t i = row_start; i < row_end; i++) {
      if (max < mat[i][0])
        max = mat[i][0];
    }

    for (size_t j = 0; j < max_col; j++) {
      if (max < mat[row_start][j])
        max = mat[row_start][j];
    }
    
    for (size_t i = row_start + 1; i < row_end; i++) {
      for (size_t j = 1; j < max_col; j++) {
        mat[i][j] = triangle(mat[i - 1][j], mat[i][j - 1], mat[i-1][j - 1], mat[i][j]);
        
        if (max < mat[i][j])
          max = mat[i][j];
      }
    }
    
    return max * max;
    } */
};

using io_util::operator<<;

/*int seq_main(int argc, char* argv[]) {
  auto m = io_util::read_matrix<unsigned>();
  Solution<unsigned> s;
  unsigned sol = s.maximal_square(m, m.size(), m[0].size());
  cout << "Solution: " << sol << endl;

  return 0;
  }*/


int main(int argc, char* argv[]) {
  std::ifstream in;

  if (argc <= 1) {
    return 1;
  }
  
  /*** Initialize MPI ***/
  //const int MASTER_RANK = 0;
  int size, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cout << "Soy rank " << rank << endl;

  std::stringstream ss;
  // partition format: matrixname.mat-10p2 (10 processes, current rank 2)
  ss << argv[1] << "-" << size << "p" << rank;
  std::string my_filename;
  ss >> my_filename;

  in = std::ifstream(my_filename);
  auto m = LocalSolution<unsigned>(io_util::read_matrix<unsigned>(in));
  
  cout << "Soy rank " << rank << " y mi matriz mide " << m.size() << endl;
  
  unsigned mi_sol = m.maximal_square();
  cout << "Soy rank " << rank << ", max_sq = " << mi_sol << endl;

  MPI_Finalize();

  return 0;
}
