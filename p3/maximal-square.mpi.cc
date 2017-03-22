#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

using std::vector;
using std::string;
using std::cin;
using std::cout;
using std::endl;

template <class T>
using matrix = vector<vector<T>>;

namespace type_util {
  matrix<int> chtoi(matrix<char> m) {
    matrix<int> r;
    std::transform(m.begin(), m.end(), std::back_inserter(r), [](vector<char> row) {
        vector<int> rr;
        std::transform(row.begin(), row.end(), std::back_inserter(rr), [](char c) { return c - '0'; });
        return rr;
      });

    return r;
  }
}

namespace io_util {
  template <class Uint>
  matrix<Uint> read_matrix(std::ifstream& in) {
    matrix<Uint> mat;
    vector<Uint> row;
    
    while (cin.good()) {
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

  template<class T>
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
template<class T>
class Proxy {
  matrix<T>& p_m;
  int row;

public:
  Proxy(matrix<T>& p_m, int row) :p_m(p_m), row(row) {}

  T& operator[](int col) {
    if (0 <= row && row < p_m.size() && 0 <= col && col < p_m[row].size())
      return p_m[row][col];
    else
      return T();
  }
  
  const T& operator[](int col) const {
    return this->operator[](col);
  }
};

template<class T>
class Wrap {
  matrix<T>& _m;

  
public:
  Wrap(matrix<T>& m) :_m(m) {}

  Proxy<T> operator[](int r) {
    return Proxy<T>(_m, r);
  }
};
*/

/*

Sobre el tipo Uint:
- debe ser un tipo de número entero sin signo (unsigned short, unsigned int, unsigned long).
- determina el tamaño máximo del cuadrado maximal: por ejemplo, usar unsigned short será válido si el lado del cuadrado maximal mide 65535 o menos.

*/
template<class Uint>
class Solution {
public:
  inline Uint triangle(Uint left, Uint right, Uint up, Uint down) {
    return (std::min(std::min(left, right), up) + 1) * down;
  }
  
  Uint maximal_square(matrix<Uint>& mat, size_t row_start, size_t row_end) {
    Uint max = 0;
    size_t max_col = mat[row_start].size();

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
  } 
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

  if (argc > 1) {
    in = std::ifstream(argv[1]);
  } else {
    return 1;
  }
  
  auto m = io_util::read_matrix<unsigned>(in);

  /*** Initialize MPI ***/
  //const int MASTER_RANK = 0;
  int size, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  cout << "Soy rank " << rank << " y mi matriz mide " << m.size() << endl;
  /*
  Solution<unsigned> s;
  size_t
    row_start = rank * m.size()/size,
    row_end = (rank < size - 1) ? (rank + 1) * m.size()/size : m[0].size();
  unsigned mi_sol = s.maximal_square(m, row_start, row_end);
  cout << "Soy rank " << rank << ", empezando desde " << row_start << ", max_sq = " << mi_sol << endl;*/

  MPI_Finalize();

  return 0;
}
