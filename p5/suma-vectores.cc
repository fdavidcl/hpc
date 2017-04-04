#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>

template <typename f>
std::vector<f> read_vector(std::ifstream& in) {
  size_t size;

  in >> size;

  std::vector<f> v(size);
  for (size_t i = 0; i < size; i++) {
    in >> v[i];
  }

  return v;
}

template<typename f>
void operator<<(std::ostream& out, const std::vector<f>& v) {
  out << v.size() << " size vector:" << std::endl;
    
  for (const f& e : v) {
    out << e << " ";
  }

  out << std::endl;
}

template <typename f>
std::vector<f> sum(const std::vector<f>& v1, const std::vector<f>& v2) {
  std::vector<f> r(v1.size());

  std::transform(v1.cbegin(), v1.cend(), v2.cbegin(), r.begin(), std::plus<f>());

  return r;
}

template <typename f>
f compare(const std::vector<f>& v1, const std::vector<f>& v2) {
  f er = 0;
  for (size_t i = 0; i < v1.size(); ++i)
    er = std::max(er, static_cast<f>(fabs(v1[i] - v2[i])));

  return er;
}

template <typename f>
bool equals(const std::vector<f>& v1, const std::vector<f>& v2, f tolerance = 1e-4) {
  bool eq = true;
  for (size_t i = 0; i < v1.size() && eq; ++i)
    eq = abs(v1[i] - v2[i]) < tolerance;

  return eq;
}

int main(int argc, char * argv[]) {
  if (argc < 4)
    return -1;

  auto
    filea = std::ifstream(argv[1]),
    fileb = std::ifstream(argv[2]),
    filec = std::ifstream(argv[3]);
  
  auto
    a = read_vector<float>(filea),
    b = read_vector<float>(fileb),
    c = read_vector<float>(filec),
    d = sum(a, b);
  
  std::cout << compare(c, d) << std::endl;
  
  return 0;
}
