#include "Timespan.h"

#include "Globals.h"

Timespan::Timespan()
: time_span_(0.0) {
}

Timespan::Timespan(const unsigned int days, const unsigned int hours,
        const unsigned int minutes, const double seconds) {

    SetValue(days, hours, minutes, seconds);
}

Timespan::Timespan(const double b) {

    time_span_ = b;
}

Timespan::Timespan(const Timespan& b) {
    time_span_ = b.time_span_;
}

Timespan::~Timespan(void) {
}

void Timespan::SetValue(const unsigned int days, const unsigned int hours,
        const unsigned int minutes, const double seconds) {

    time_span_ = static_cast<double> (days);
    AddHours(hours);
    AddMinutes(minutes);
    AddSeconds(seconds);
}

void Timespan::AddDays(const unsigned int days) {
    time_span_ += static_cast<double> (days);
}

void Timespan::AddHours(const unsigned int hours) {
    time_span_ += (static_cast<double> (hours) / kHOURS_PER_DAY);
}

void Timespan::AddMinutes(const unsigned int minutes) {
    time_span_ += (static_cast<double> (minutes) / kMINUTES_PER_DAY);
}

void Timespan::AddSeconds(const double seconds) {
    time_span_ += (seconds / kSECONDS_PER_DAY);
}

double Timespan::GetTotalDays() const {
    return time_span_;
}

double Timespan::GetTotalHours() const {
    return time_span_ * kHOURS_PER_DAY;
}

double Timespan::GetTotalMinutes() const {
    return time_span_ * kMINUTES_PER_DAY;
}

double Timespan::GetTotalSeconds() const {
    return time_span_ * kSECONDS_PER_DAY;
}

Timespan& Timespan::operator =(const Timespan& b) {

    if (this != &b) {
        time_span_ = b.time_span_;
    }
    return (*this);
}

Timespan Timespan::operator +(const Timespan& b) const {

    return Timespan(*this) += b;
}

Timespan Timespan::operator -(const Timespan& b) const {

    return Timespan(*this) -= b;
}

Timespan Timespan::operator/(const double b) const {

    return Timespan(*this) /= b;
}

Timespan Timespan::operator*(const double b) const {

    return Timespan(*this) *= b;
}

Timespan & Timespan::operator+=(const Timespan& b) {

    time_span_ += b.time_span_;
    return (*this);
}

Timespan & Timespan::operator-=(const Timespan& b) {

    time_span_ -= b.time_span_;
    return (*this);
}

Timespan & Timespan::operator/=(const double b) {

    time_span_ /= b;
    return (*this);
}

Timespan & Timespan::operator*=(const double b) {

    time_span_ *= b;
    return (*this);
}

bool Timespan::operator ==(const Timespan& b) const {

    if (time_span_ == b.time_span_)
        return true;
    else
        return false;
}

bool Timespan::operator !=(const Timespan& b) const {

    return !(*this == b);
}

bool Timespan::operator>(const Timespan& b) const {

    if (time_span_ > b.time_span_)
        return true;
    else
        return false;
}

bool Timespan::operator<(const Timespan& b) const {

    if (time_span_ < b.time_span_)
        return true;
    else
        return false;
}

bool Timespan::operator >=(const Timespan& b) const {

    if (time_span_ >= b.time_span_)
        return true;
    else
        return false;
}

bool Timespan::operator <=(const Timespan & b) const {

    if (time_span_ <= b.time_span_)
        return true;
    else
        return false;
}

double& operator +=(double& a, const Timespan& b) {

    a += b.time_span_;
    return a;
}

double& operator -=(double& a, const Timespan& b) {

    a -= b.time_span_;
    return a;
}