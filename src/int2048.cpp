#include "int2048.h"

namespace sjtu {

int2048::int2048() : digits_(1, 0), negative_(false) {}

int2048::int2048(long long value) : digits_(1, 0), negative_(false) {
  if (value < 0) {
    negative_ = true;
  }
  unsigned long long x = 0;
  if (value < 0) {
    x = static_cast<unsigned long long>(-(value + 1)) + 1ULL;
  } else {
    x = static_cast<unsigned long long>(value);
  }
  digits_.clear();
  if (x == 0) {
    digits_.push_back(0);
  } else {
    while (x > 0) {
      digits_.push_back(static_cast<int>(x % BASE));
      x /= BASE;
    }
  }
  normalize();
}

int2048::int2048(const std::string &s) : digits_(1, 0), negative_(false) { read(s); }

int2048::int2048(const int2048 &other) : digits_(other.digits_), negative_(other.negative_) {}

bool int2048::isZero() const { return digits_.size() == 1 && digits_[0] == 0; }

void int2048::normalize() {
  while (digits_.size() > 1 && digits_.back() == 0) {
    digits_.pop_back();
  }
  if (digits_.empty()) {
    digits_.push_back(0);
  }
  if (isZero()) {
    negative_ = false;
  }
}

void int2048::trimVector(std::vector<int> &v) {
  while (v.size() > 1 && v.back() == 0) {
    v.pop_back();
  }
  if (v.empty()) {
    v.push_back(0);
  }
}

int int2048::cmpVector(const std::vector<int> &a, const std::vector<int> &b) {
  if (a.size() != b.size()) {
    return a.size() < b.size() ? -1 : 1;
  }
  for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
    if (a[i] != b[i]) {
      return a[i] < b[i] ? -1 : 1;
    }
  }
  return 0;
}

std::vector<int> int2048::mulVectorInt(const std::vector<int> &a, int m) {
  if ((a.size() == 1 && a[0] == 0) || m == 0) {
    return std::vector<int>(1, 0);
  }
  std::vector<int> res;
  res.resize(a.size() + 1, 0);
  long long carry = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    long long cur = carry + 1LL * a[i] * m;
    res[i] = static_cast<int>(cur % BASE);
    carry = cur / BASE;
  }
  if (carry > 0) {
    res[a.size()] = static_cast<int>(carry);
  }
  trimVector(res);
  return res;
}

std::vector<int> int2048::divVectorInt(const std::vector<int> &a, int v) {
  std::vector<int> res(a.size(), 0);
  long long rem = 0;
  for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
    long long cur = rem * BASE + a[i];
    res[i] = static_cast<int>(cur / v);
    rem = cur % v;
  }
  trimVector(res);
  return res;
}

int int2048::absCompare(const int2048 &lhs, const int2048 &rhs) { return cmpVector(lhs.digits_, rhs.digits_); }

void int2048::addAbs(const int2048 &rhs) {
  if (rhs.isZero()) {
    return;
  }
  const size_t n = digits_.size() > rhs.digits_.size() ? digits_.size() : rhs.digits_.size();
  if (digits_.size() < n) {
    digits_.resize(n, 0);
  }
  int carry = 0;
  for (size_t i = 0; i < n; ++i) {
    int cur = digits_[i] + carry;
    if (i < rhs.digits_.size()) {
      cur += rhs.digits_[i];
    }
    if (cur >= BASE) {
      cur -= BASE;
      carry = 1;
    } else {
      carry = 0;
    }
    digits_[i] = cur;
  }
  if (carry != 0) {
    digits_.push_back(carry);
  }
}

void int2048::subAbs(const int2048 &rhs) {
  int borrow = 0;
  for (size_t i = 0; i < digits_.size(); ++i) {
    int cur = digits_[i] - borrow;
    if (i < rhs.digits_.size()) {
      cur -= rhs.digits_[i];
    }
    if (cur < 0) {
      cur += BASE;
      borrow = 1;
    } else {
      borrow = 0;
    }
    digits_[i] = cur;
  }
  normalize();
}

std::vector<int> int2048::multiplySimple(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> res(a.size() + b.size(), 0);
  for (size_t i = 0; i < a.size(); ++i) {
    long long carry = 0;
    for (size_t j = 0; j < b.size() || carry > 0; ++j) {
      long long cur = res[i + j] + carry;
      if (j < b.size()) {
        cur += 1LL * a[i] * b[j];
      }
      res[i + j] = static_cast<int>(cur % BASE);
      carry = cur / BASE;
    }
  }
  trimVector(res);
  return res;
}

void int2048::fft(std::vector<std::complex<double>> &a, bool invert) {
  const int n = static_cast<int>(a.size());
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::complex<double> tmp = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
  }

  const double PI = 3.141592653589793238462643383279502884;
  for (int len = 2; len <= n; len <<= 1) {
    const double ang = 2.0 * PI / len * (invert ? -1.0 : 1.0);
    const std::complex<double> wlen = std::polar(1.0, ang);
    for (int i = 0; i < n; i += len) {
      std::complex<double> w(1.0, 0.0);
      const int half = len >> 1;
      for (int j = 0; j < half; ++j) {
        std::complex<double> u = a[i + j];
        std::complex<double> v = a[i + j + half] * w;
        a[i + j] = u + v;
        a[i + j + half] = u - v;
        w *= wlen;
      }
    }
  }

  if (invert) {
    for (int i = 0; i < n; ++i) {
      a[i] /= static_cast<double>(n);
    }
  }
}

std::vector<int> int2048::multiplyFFT(const std::vector<int> &a, const std::vector<int> &b) {
  int n = 1;
  while (n < static_cast<int>(a.size() + b.size())) {
    n <<= 1;
  }

  std::vector<std::complex<double>> fa(n);
  std::vector<std::complex<double>> fb(n);
  for (size_t i = 0; i < a.size(); ++i) {
    fa[i] = std::complex<double>(a[i], 0.0);
  }
  for (size_t i = 0; i < b.size(); ++i) {
    fb[i] = std::complex<double>(b[i], 0.0);
  }

  fft(fa, false);
  fft(fb, false);
  for (int i = 0; i < n; ++i) {
    fa[i] *= fb[i];
  }
  fft(fa, true);

  std::vector<int> res(a.size() + b.size(), 0);
  long long carry = 0;
  for (size_t i = 0; i < res.size(); ++i) {
    double real_part = fa[i].real();
    long long cur = real_part >= 0.0 ? static_cast<long long>(real_part + 0.5) : static_cast<long long>(real_part - 0.5);
    cur += carry;
    res[i] = static_cast<int>(cur % BASE);
    carry = cur / BASE;
  }
  while (carry > 0) {
    res.push_back(static_cast<int>(carry % BASE));
    carry /= BASE;
  }

  trimVector(res);
  return res;
}

std::vector<int> int2048::multiplyDigits(const std::vector<int> &a, const std::vector<int> &b) {
  if ((a.size() == 1 && a[0] == 0) || (b.size() == 1 && b[0] == 0)) {
    return std::vector<int>(1, 0);
  }
  if (a.size() < 64 || b.size() < 64) {
    return multiplySimple(a, b);
  }
  return multiplyFFT(a, b);
}

void int2048::divmodAbsVectors(const std::vector<int> &a, const std::vector<int> &b, std::vector<int> &q, std::vector<int> &r) {
  if (cmpVector(a, b) < 0) {
    q.assign(1, 0);
    r = a;
    return;
  }

  if (b.size() == 1) {
    q = divVectorInt(a, b[0]);
    long long rem = 0;
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
      long long cur = rem * BASE + a[i];
      rem = cur % b[0];
    }
    r.assign(1, static_cast<int>(rem));
    trimVector(q);
    trimVector(r);
    return;
  }

  const int norm = BASE / (b.back() + 1);
  std::vector<int> u = mulVectorInt(a, norm);
  std::vector<int> v = mulVectorInt(b, norm);
  u.push_back(0);

  const int n = static_cast<int>(u.size());
  const int m = static_cast<int>(v.size());
  q.assign(n - m, 0);

  for (int j = n - m - 1; j >= 0; --j) {
    long long numerator = 1LL * u[j + m] * BASE + u[j + m - 1];
    long long qhat = numerator / v[m - 1];
    long long rhat = numerator % v[m - 1];

    if (qhat >= BASE) {
      qhat = BASE - 1;
    }

    while (qhat * v[m - 2] > 1LL * BASE * rhat + u[j + m - 2]) {
      --qhat;
      rhat += v[m - 1];
      if (rhat >= BASE) {
        break;
      }
    }

    long long borrow = 0;
    long long carry = 0;
    for (int i = 0; i < m; ++i) {
      long long prod = qhat * v[i] + carry;
      carry = prod / BASE;
      long long sub = static_cast<long long>(u[j + i]) - (prod % BASE) - borrow;
      if (sub < 0) {
        sub += BASE;
        borrow = 1;
      } else {
        borrow = 0;
      }
      u[j + i] = static_cast<int>(sub);
    }

    long long sub = static_cast<long long>(u[j + m]) - carry - borrow;
    if (sub < 0) {
      --qhat;
      long long carry2 = 0;
      for (int i = 0; i < m; ++i) {
        long long sum = static_cast<long long>(u[j + i]) + v[i] + carry2;
        if (sum >= BASE) {
          sum -= BASE;
          carry2 = 1;
        } else {
          carry2 = 0;
        }
        u[j + i] = static_cast<int>(sum);
      }
      u[j + m] += static_cast<int>(carry2);
    } else {
      u[j + m] = static_cast<int>(sub);
    }

    q[j] = static_cast<int>(qhat);
  }

  trimVector(q);
  r.assign(u.begin(), u.begin() + m);
  r = divVectorInt(r, norm);
  trimVector(r);
}

void int2048::divmodFloor(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r) {
  int2048 aa(a);
  int2048 bb(b);
  aa.negative_ = false;
  bb.negative_ = false;

  std::vector<int> qabs;
  std::vector<int> rabs;
  divmodAbsVectors(aa.digits_, bb.digits_, qabs, rabs);

  int2048 q0;
  q0.digits_ = qabs;
  q0.negative_ = false;
  q0.normalize();

  int2048 r0;
  r0.digits_ = rabs;
  r0.negative_ = false;
  r0.normalize();

  const bool sign_diff = (a.negative_ != b.negative_);
  if (!sign_diff) {
    q = q0;
    q.negative_ = false;
    q.normalize();

    r = r0;
    if (!r.isZero() && b.negative_) {
      r.negative_ = true;
    }
    r.normalize();
    return;
  }

  if (r0.isZero()) {
    q = q0;
    if (!q.isZero()) {
      q.negative_ = true;
    }
    q.normalize();

    r = r0;
    r.normalize();
    return;
  }

  q = q0;
  q += int2048(1);
  q.negative_ = true;
  q.normalize();

  r = bb;
  r -= r0;
  if (!r.isZero() && b.negative_) {
    r.negative_ = true;
  }
  r.normalize();
}

void int2048::read(const std::string &s) {
  digits_.assign(1, 0);
  negative_ = false;

  if (s.empty()) {
    return;
  }

  int pos = 0;
  if (s[0] == '-') {
    negative_ = true;
    pos = 1;
  } else if (s[0] == '+') {
    pos = 1;
  }

  while (pos < static_cast<int>(s.size()) && s[pos] == '0') {
    ++pos;
  }

  if (pos >= static_cast<int>(s.size())) {
    digits_[0] = 0;
    negative_ = false;
    return;
  }

  digits_.clear();
  for (int i = static_cast<int>(s.size()); i > pos; i -= BASE_DIGS) {
    int left = i - BASE_DIGS;
    if (left < pos) {
      left = pos;
    }
    int x = 0;
    for (int j = left; j < i; ++j) {
      x = x * 10 + (s[j] - '0');
    }
    digits_.push_back(x);
  }

  normalize();
}

void int2048::print() { std::cout << *this; }

int2048 &int2048::add(const int2048 &rhs) { return (*this) += rhs; }

int2048 add(int2048 lhs, const int2048 &rhs) { return lhs.add(rhs); }

int2048 &int2048::minus(const int2048 &rhs) { return (*this) -= rhs; }

int2048 minus(int2048 lhs, const int2048 &rhs) { return lhs.minus(rhs); }

int2048 int2048::operator+() const { return *this; }

int2048 int2048::operator-() const {
  int2048 res(*this);
  if (!res.isZero()) {
    res.negative_ = !res.negative_;
  }
  return res;
}

int2048 &int2048::operator=(const int2048 &rhs) {
  if (this == &rhs) {
    return *this;
  }
  digits_ = rhs.digits_;
  negative_ = rhs.negative_;
  return *this;
}

int2048 &int2048::operator+=(const int2048 &rhs) {
  if (rhs.isZero()) {
    return *this;
  }
  if (isZero()) {
    *this = rhs;
    return *this;
  }

  if (negative_ == rhs.negative_) {
    addAbs(rhs);
    normalize();
    return *this;
  }

  int cmp = absCompare(*this, rhs);
  if (cmp == 0) {
    digits_.assign(1, 0);
    negative_ = false;
    return *this;
  }

  if (cmp > 0) {
    subAbs(rhs);
  } else {
    int2048 tmp(rhs);
    tmp.subAbs(*this);
    *this = tmp;
  }
  normalize();
  return *this;
}

int2048 operator+(int2048 lhs, const int2048 &rhs) {
  lhs += rhs;
  return lhs;
}

int2048 &int2048::operator-=(const int2048 &rhs) {
  *this += (-rhs);
  return *this;
}

int2048 operator-(int2048 lhs, const int2048 &rhs) {
  lhs -= rhs;
  return lhs;
}

int2048 &int2048::operator*=(const int2048 &rhs) {
  if (isZero() || rhs.isZero()) {
    digits_.assign(1, 0);
    negative_ = false;
    return *this;
  }

  negative_ = (negative_ != rhs.negative_);
  digits_ = multiplyDigits(digits_, rhs.digits_);
  normalize();
  return *this;
}

int2048 operator*(int2048 lhs, const int2048 &rhs) {
  lhs *= rhs;
  return lhs;
}

int2048 &int2048::operator/=(const int2048 &rhs) {
  int2048 q;
  int2048 r;
  divmodFloor(*this, rhs, q, r);
  *this = q;
  return *this;
}

int2048 operator/(int2048 lhs, const int2048 &rhs) {
  lhs /= rhs;
  return lhs;
}

int2048 &int2048::operator%=(const int2048 &rhs) {
  int2048 q;
  int2048 r;
  divmodFloor(*this, rhs, q, r);
  *this = r;
  return *this;
}

int2048 operator%(int2048 lhs, const int2048 &rhs) {
  lhs %= rhs;
  return lhs;
}

std::istream &operator>>(std::istream &is, int2048 &value) {
  std::string s;
  is >> s;
  value.read(s);
  return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &value) {
  if (value.negative_ && !value.isZero()) {
    os << '-';
  }
  os << value.digits_.back();
  for (int i = static_cast<int>(value.digits_.size()) - 2; i >= 0; --i) {
    int x = value.digits_[i];
    if (x < 1000) {
      os << '0';
    }
    if (x < 100) {
      os << '0';
    }
    if (x < 10) {
      os << '0';
    }
    os << x;
  }
  return os;
}

bool operator==(const int2048 &lhs, const int2048 &rhs) { return lhs.negative_ == rhs.negative_ && lhs.digits_ == rhs.digits_; }

bool operator!=(const int2048 &lhs, const int2048 &rhs) { return !(lhs == rhs); }

bool operator<(const int2048 &lhs, const int2048 &rhs) {
  if (lhs.negative_ != rhs.negative_) {
    return lhs.negative_;
  }
  int cmp = int2048::absCompare(lhs, rhs);
  if (!lhs.negative_) {
    return cmp < 0;
  }
  return cmp > 0;
}

bool operator>(const int2048 &lhs, const int2048 &rhs) { return rhs < lhs; }

bool operator<=(const int2048 &lhs, const int2048 &rhs) { return !(rhs < lhs); }

bool operator>=(const int2048 &lhs, const int2048 &rhs) { return !(lhs < rhs); }

} // namespace sjtu
