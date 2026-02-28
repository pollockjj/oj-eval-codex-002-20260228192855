#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// Do not use "using namespace std;"

namespace sjtu {
class int2048 {
private:
  static const int BASE = 100;
  static const int BASE_DIGS = 2;

  std::vector<int> digits_;
  bool negative_;

  bool isZero() const;
  void normalize();

  static int absCompare(const int2048 &, const int2048 &);
  static std::vector<int> multiplySimple(const std::vector<int> &, const std::vector<int> &);
  static void fft(std::vector<std::complex<double>> &, bool);
  static std::vector<int> multiplyFFT(const std::vector<int> &, const std::vector<int> &);
  static std::vector<int> multiplyDigits(const std::vector<int> &, const std::vector<int> &);

  static void trimVector(std::vector<int> &);
  static std::vector<int> mulVectorInt(const std::vector<int> &, int);
  static std::vector<int> divVectorInt(const std::vector<int> &, int);
  static int cmpVector(const std::vector<int> &, const std::vector<int> &);
  static void divmodAbsVectors(const std::vector<int> &, const std::vector<int> &, std::vector<int> &, std::vector<int> &);
  static void divmodFloor(const int2048 &, const int2048 &, int2048 &, int2048 &);

  void addAbs(const int2048 &);
  void subAbs(const int2048 &);

public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif
