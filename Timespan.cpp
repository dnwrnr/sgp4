#include "Timespan.h"

#include "Globals.h"

Timespan::Timespan(void) {
}

Timespan::Timespan(const double time_span) {
}

Timespan::Timespan(const Timespan& b) {

    time_span_ = b.time_span_;
}

Timespan::~Timespan(void) {
}

Timespan& Timespan::operator =(const Timespan& b) {

    if (this != &b) {
        time_span_ = b.time_span_;

    }
    return (*this);
}

Timespan Timespan::operator +(const Timespan& b) const {

    Timespan result(*this);
    result.time_span_ += b.GetTotalDays();
    return result;
}

Timespan Timespan::operator -(const Timespan& b) const {

    Timespan result(*this);
    result.time_span_ -= b.GetTotalDays();
    return result;
}

const Timespan & Timespan::operator+=(const Timespan& b) {

    time_span_ += b.time_span_;
    return (*this);
}

const Timespan & Timespan::operator-=(const Timespan& b) {

    time_span_ -= b.time_span_;
    return (*this);
}

Timespan Timespan::operator +() const {

    Timespan result(*this);
    result.time_span_ = +result.time_span_;
    return result;
}

Timespan Timespan::operator -() const {

    Timespan result(*this);
    result.time_span_ = -result.time_span_;
    return result;
}

bool Timespan::operator ==(const Timespan& b) const {

    if (time_span_ == b.time_span_)
        return true;
    else
        return false;
}

bool Timespan::operator !=(const Timespan& b) const {

    if (time_span_ == b.time_span_)
        return false;
    else
        return true;
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

double Timespan::GetTotalDays() const{
    return time_span_ / Globals::MIN_PER_DAY();
}