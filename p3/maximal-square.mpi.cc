#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>

using std::vector;
using std::string;
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
  std::ostream& operator<<(std::ostream& out, const vector<T>& row) {
    for (const T& e : row) {
      out << e << " ";
    }

    return out;
  }

  template<typename T>
  std::ostream& operator<<(std::ostream& out, const matrix<T>& mat) {
    out << mat.size() << "x" << mat.at(0).size() << " matrix:" << endl;
    
    for (const auto& row : mat) {
      out << row << endl;
    }

    return out;
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
  
  inline Uint triangle(const Uint& left, const Uint& right, const Uint& up, const Uint& down) const {
    return (std::min(std::min(left, right), up) + 1) * (down != 0);
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

  Uint use_top_row(const vector<Uint>& top) {
    if (!locally_solved) {
      maximal_square();
    }

    Uint max = 0;
    size_t max_col = (*this)[0].size();

    // vector<bool> is pretty slow. May consider changing it for vector<unsigned char> if storage constraints aren't very strong
    vector<bool> remaining(top.size(), true);
    //size_t count = std::count_if(remaining.begin(), remaining.end(), [](const auto& e) { return e != 0; });
    size_t count = top.size();
    
    for (size_t j = 0; j < max_col; ++j) {      
      bool dif = remaining[j] && ((*this)[0][j] == 0);
      count -= dif;
      remaining[j] = remaining[j] && ((*this)[0][j] != 0);
      
      if (remaining[j]) {
        (*this)[0][j] = triangle(top[j], (*this)[0][j - 1], top[j - 1], (*this)[0][j]);
        
        if (max < (*this)[0][j])
          max = (*this)[0][j];
      }
    }

    for (size_t i = 1; i < this->size() && count > 0; ++i) {
      for (size_t j = 0; j < max_col && count > 0; ++j) {
        bool dif = remaining[j] & ((*this)[i][j] == 0);
        count -= dif;
        remaining[j] = remaining[j] & ((*this)[i][j] != 0);
        
        if (remaining[j]) {
          (*this)[i][j] = triangle((*this)[i - 1][j], (*this)[i][j - 1], (*this)[i-1][j - 1], (*this)[i][j]);
        
          if (max < (*this)[i][j])
            max = (*this)[i][j];
        }
      }
    }
    
    max_square_size = std::max(max * max, max_square_size);

    return max_square_size;
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

  if (argc <= 1) {
    return 1;
  }
  
  /*** Initialize MPI ***/
  int size, rank;
  const int ROOT_RANK = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::stringstream ss;
  // partition format: matrixname.mat-10p2 (10 processes, current rank 2)
  ss << argv[1] << "-" << size << "p" << rank;
  std::string my_filename;
  ss >> my_filename;

  in = std::ifstream(my_filename);
  auto m = LocalSolution<unsigned>(io_util::read_matrix<unsigned>(in));

  auto start = std::chrono::high_resolution_clock::now();
  int local_sol = m.maximal_square();

  //MPI_Request handle;
  MPI_Status status;
  
  vector<unsigned> top(m[0].size());
  
  if (rank > 0) {
    MPI_Recv(&top[0], top.size(), MPI_UNSIGNED, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    local_sol = m.use_top_row(top);
  }

  if (rank < size - 1) {
    MPI_Send(&m[0][0], m[0].size(), MPI_UNSIGNED, rank + 1, 0, MPI_COMM_WORLD);
  }

  int global_sol;
  MPI_Reduce(&local_sol, &global_sol, 1, MPI_INT, MPI_MAX, ROOT_RANK, MPI_COMM_WORLD);

  auto end = std::chrono::high_resolution_clock::now();

  std::cerr << "Soy rank " << rank << ", max_sq = " << local_sol << endl;
  
  if (rank == ROOT_RANK) {
    std::cerr << "Solución global: " << global_sol << endl;
    std::cout << size << "," << size * m.size() << ","
              << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << endl;
  }

  MPI_Finalize();

  return 0;
}
