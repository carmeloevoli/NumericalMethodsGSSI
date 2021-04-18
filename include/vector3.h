// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#ifndef INCLUDE_VECTOR3_H_
#define INCLUDE_VECTOR3_H_

#include <cmath>
#include <fstream>

namespace NM {

template <typename T>
class Vector3 {
 public:
  T x, y, z;

  Vector3() : x(0), y(0), z(0) {}

  // Provides implicit conversion
  template <typename U>
  Vector3(const Vector3<U> &v) : x(v.x), y(v.y), z(v.z) {}

  explicit Vector3(const T &X, const T &Y, const T &Z) : x(X), y(Y), z(Z) {}

  explicit Vector3(T t) : x(t), y(t), z(t) {}

  void setX(const T X) { x = X; }

  void setY(const T Y) { y = Y; }

  void setZ(const T Z) { z = Z; }

  void setXYZ(const T X, const T Y, const T Z) {
    x = X;
    y = Y;
    z = Z;
  }

  T getX() const { return (x); }

  T getY() const { return (y); }

  T getZ() const { return (z); }

  // magnitude (2-norm) of the vector
  T getModule() const { return (std::sqrt(x * x + y * y + z * z)); }

  // return the unit-vector e_r
  Vector3<T> getUnitVector() const { return (*this / getModule()); }

  // return vector with absolute values
  Vector3<T> abs() const { return (Vector3<T>(std::abs(x), std::abs(y), std::abs(z))); }

  // return vector with floored values
  Vector3<T> floor() const { return (Vector3<T>(std::floor(x), std::floor(y), std::floor(z))); }

  // return vector with ceiled values
  Vector3<T> ceil() const { return (Vector3<T>(std::ceil(x), std::ceil(y), std::ceil(z))); }

  // minimum element
  T min() const { return (std::min(x, std::min(y, z))); }

  // maximum element
  T max() const { return (std::max(x, std::max(y, z))); }

  // dot product
  T dot(const Vector3<T> &v) const { return (x * v.x + y * v.y + z * v.z); }

  // cross product
  Vector3<T> cross(const Vector3<T> &v) const {
    return (Vector3<T>(y * v.z - v.y * z, z * v.x - v.z * x, x * v.y - v.x * y));
  }

  // returns true if all elements of the first vector are smaller than those in
  // the second vector
  bool operator<(const Vector3<T> &v) const {
    if (x > v.x)
      return (false);
    else if (x < v.x)
      return (true);
    if (y > v.y)
      return (false);
    else if (y < v.y)
      return (true);
    if (z >= v.z)
      return (false);
    else
      return (true);
  }

  // returns true if all elements of the two vectors are equal
  bool operator==(const Vector3<T> &v) const {
    if (x != v.x) return (false);
    if (y != v.y) return (false);
    if (z != v.z) return (false);
    return (true);
  }

  Vector3<T> operator+(const Vector3<T> &v) const { return (Vector3(x + v.x, y + v.y, z + v.z)); }

  Vector3<T> operator+(const T &f) const { return (Vector3(x + f, y + f, z + f)); }

  Vector3<T> operator-(const Vector3<T> &v) const { return (Vector3(x - v.x, y - v.y, z - v.z)); }

  Vector3<T> operator-(const T &f) const { return (Vector3(x - f, y - f, z - f)); }

  // element-wise multiplication
  Vector3<T> operator*(const Vector3<T> &v) const { return (Vector3(x * v.x, y * v.y, z * v.z)); }

  Vector3<T> operator*(const T &v) const { return (Vector3(x * v, y * v, z * v)); }

  // element-wise division
  Vector3<T> operator/(const Vector3<T> &v) const { return (Vector3(x / v.x, y / v.y, z / v.z)); }

  Vector3<T> operator/(const T &f) const { return (Vector3(x / f, y / f, z / f)); }

  // element-wise modulo operation
  Vector3<T> operator%(const Vector3<T> &v) const { return (Vector3(fmod(x, v.x), fmod(y, v.y), fmod(z, v.z))); }

  Vector3<T> operator%(const T &f) const { return (Vector3(fmod(x, f), fmod(y, f), fmod(z, f))); }

  Vector3<T> &operator-=(const Vector3<T> &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return (*this);
  }

  Vector3<T> &operator-=(const T &f) {
    x -= f;
    y -= f;
    z -= f;
    return (*this);
  }

  Vector3<T> &operator+=(const Vector3<T> &v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return (*this);
  }

  Vector3<T> &operator+=(const T &f) {
    x += f;
    y += f;
    z += f;
    return (*this);
  }

  // element-wise multiplication
  Vector3<T> &operator*=(const Vector3<T> &v) {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return (*this);
  }

  Vector3<T> &operator*=(const T &f) {
    x *= f;
    y *= f;
    z *= f;
    return (*this);
  }

  // element-wise division
  Vector3<T> &operator/=(const Vector3<T> &v) {
    x /= v.x;
    y /= v.y;
    z /= v.z;
    return (*this);
  }

  Vector3<T> &operator/=(const T &f) {
    x /= f;
    y /= f;
    z /= f;
    return (*this);
  }

  // element-wise modulo operation
  Vector3<T> &operator%=(const Vector3<T> &v) {
    x = fmod(x, v.x);
    y = fmod(y, v.y);
    z = fmod(z, v.z);
    return (*this);
  }

  Vector3<T> &operator%=(const T &f) {
    x = fmod(x, f);
    y = fmod(y, f);
    z = fmod(z, f);
    return (*this);
  }

  Vector3<T> &operator=(const Vector3<T> &v) {
    x = v.x;
    y = v.y;
    z = v.z;
    return (*this);
  }

  Vector3<T> &operator=(const T &f) {
    x = f;
    y = f;
    z = f;
    return (*this);
  }
};

template <typename T>
inline std::ostream &operator<<(std::ostream &out, const Vector3<T> &v) {
  out << v.x << " " << v.y << " " << v.z;
  return (out);
}

template <typename T>
inline std::istream &operator>>(std::istream &in, Vector3<T> &v) {
  in >> v.x >> v.y >> v.z;
  return (in);
}

template <typename T>
inline Vector3<T> operator*(T f, const Vector3<T> &v) {
  return (Vector3<T>(v.x * f, v.y * f, v.z * f));
}

using Vector3d = Vector3<double>;
using Vector3f = Vector3<float>;

}  // namespace NM

#endif  // INCLUDE_VECTOR3_H_
