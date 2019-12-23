#pragma once
#include <math.h>

#if 0
#if __cplusplus > 199711L
#define constexpr  constexpr
#else
#define constexpr  const
#endif
#endif

//! @brief Set of mathematical static functions. 
//! @tparam _float_type The real type to use.
template <typename _float_type> struct Math
{
  //! @brief Compute the squared root.
  //! @param x_ The value from which we compute the squared root.
  //! @return The squared root.
  static constexpr _float_type Sqrt(const _float_type&x_) noexcept;
  //! @brief Compute the absolute value.
  //! @param x_ The value from which we compute the absolue value.
  //! @return The absolute value.
  static constexpr _float_type Abs(const _float_type&x_) noexcept;
  //! @brief Compute the power x^y.
  //! @param x_ The base of the power.
  //! @param y_ The exponent of the power.
  //! @return The power x^y.
  static constexpr _float_type Pow(const _float_type&x_,const _float_type&y_) noexcept;
  //! @brief Compute the maximum between two real.
  //! @param a_ The first operand.
  //! @param b_ The second operand.
  //! @return MAX(a,b)
  static constexpr _float_type Max(const _float_type&a_,const _float_type&b_) noexcept;
  //! @brief Compute the maximum between two real.
  //! @param a_ The first operand.
  //! @param b_ The second operand.
  //! @return MAX(a,b)
  static constexpr _float_type Min(const _float_type&a_,const _float_type&b_) noexcept;
#if 0
  //! @brief Get the maximum value of the _float_type.
  //! @return The maximum value of the _float_type.
  static constexpr _float_type MaxValue() noexcept;
  //! @brief Get the minimum value of the _float_type.
  //! @return The minimum value of the _float_type.
  static constexpr _float_type MinValue() noexcept;
#endif
  //! @brief Set the value to zero if the absolute value is less than a tolerance.
  //! @param The value.
  //! @param The tolerance.
  //! @return The filtered value.
  static constexpr _float_type Filter(const _float_type&value_,const _float_type&tolerance_) noexcept;

};

//! @copydoc Math
template <> struct Math<double>
{
  static constexpr unsigned int NumDigits = 15;
  static constexpr double Zero = 0.0;
  static constexpr double MachineEpsilon = 2.22044604925031308e-16;
  static constexpr double machineEpsilon = 2.22044604925031308e-16;
  //! @copydoc Math::Pow(const _float_type&x_,const _float_type&y_)
  static constexpr double Pow(const double&x_,const double&y_) noexcept
  {
    return pow(x_,y_);
  };
  //! @copydoc Math::Sqrt(const _float_type&x_)
  static constexpr double Sqrt(const double&x_) noexcept
  {
    return  sqrt(x_);
  };
  //! @copydoc Math::Abs(const _float_type&x_)
  static constexpr double Abs(const double&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };
#if 1
  //! @copydoc Math::MaxValue()
  static constexpr double MaxValue() noexcept { return DBL_MAX; };
  //! @copydoc Math::MinValue()
  static constexpr double MinValue() noexcept { return DBL_MIN; };
#endif
  //! @copydoc Math::Max(const _float_type&,const _float_type&)
  static constexpr double Max(const double&a_,const double&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc Math::Min(const _float_type&,const _float_type&)
  static constexpr double Min(const double&a_,const double&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc Math::Filter(const _float_type&,const _float_type&)
  static constexpr double Filter(const double&value_,const double&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0  );
  };

};
//constexpr double Math<double>::MachineEpsilon;
//constexpr double Math<double>::machineEpsilon;


//! @copydoc Math
template <> struct Math<float>
{
  static constexpr unsigned int NumDigits = 7;
  static constexpr float Zero = 0.0f;
  static constexpr float MachineEpsilon = 1.19209290e-07f;
  //! @copydoc Math::Pow(const _float_type&x_,const _float_type&y_)
  static constexpr float Pow(const float&x_,const float&y_) noexcept
  {
    return powf(x_,y_);
  };
  //! @copydoc Math::Sqrt(const _float_type&x_)
  static constexpr float Sqrt(const float&x_) noexcept
  {
    return sqrtf(x_);
  };
  //! @copydoc Math::Abs(const _float_type&x_)
  static constexpr float Abs(const float&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };
  //! @copydoc Math::MaxValue()
#if 0
  static constexpr float MaxValue() noexcept { return FLT_MAX; };
  //! @copydoc Math::MinValue()
  static constexpr float MinValue() noexcept { return FLT_MIN; };
#endif
  //! @copydoc Math::Max(const _float_type&,const _float_type&)
  static constexpr float Max(const float&a_,const float&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc Math::Min(const _float_type&,const _float_type&)
  static constexpr float Min(const float&a_,const float&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc Math::Filter(const _float_type&,const _float_type&)
  static constexpr float Filter(const float &value_,const float&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0f  );
  };

};
//constexpr float Math<float>::MachineEpsilon;

//! @copydoc Math
template <> struct Math<long double>
{
  static constexpr unsigned int NumDigits = 17;
  static constexpr long double Zero = 0.0L;
  static constexpr long double MachineEpsilon = 1.08420217248550443401e-19L;
  //! @copydoc Math::Pow(const _float_type&x_,const _float_type&y_)
  static constexpr long double Pow(const long double&x_,const long double&y_) noexcept
  {
    return powl(x_,y_);
  };
  //! @copydoc Math::Sqrt(const _float_type&x_)
  static constexpr long double Sqrt(const long double&x_) noexcept
  {
    return sqrtl(x_);
  };
  //! @copydoc Math::Abs(const _float_type&x_)
  static constexpr long double Abs(const long double&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };
#if 0
  //! @copydoc Math::MaxValue()
  static constexpr long double MaxValue() noexcept { return LDBL_MAX; };
  //! @copydoc Math::MinValue()
  static constexpr long double MinValue() noexcept { return LDBL_MIN; };
#endif
  //! @copydoc Math::Max(const _float_type&,const _float_type&)
  static constexpr long double Max(const long double&a_,const long double&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc Math::Min(const _float_type&,const _float_type&)
  static constexpr long double Min(const long double&a_,const long double&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc Math::Filter(const _float_type&,const _float_type&)
  static constexpr long double Filter(const long double &value_,const long double&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0L  );
  };

};
//constexpr long double Math<long double>::MachineEpsilon;
