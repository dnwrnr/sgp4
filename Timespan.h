#ifndef TIMESPAN_H_
#define TIMESPAN_H_

class Timespan {
public:
    Timespan();
    Timespan(const unsigned int days, const unsigned int hours,
            const unsigned int minutes, const double seconds);
    Timespan(const double b);
    Timespan(const Timespan& b);
    virtual ~Timespan(void);

    void SetValue(const unsigned int days, const unsigned int hours,
            const unsigned int minutes, const double seconds);

    void AddDays(const unsigned int days);
    void AddHours(const unsigned int hours);
    void AddMinutes(const unsigned int minutes);
    void AddSeconds(const double seconds);

    double GetTotalDays() const;
    double GetTotalHours() const;
    double GetTotalMinutes() const;
    double GetTotalSeconds() const;

    // assignment
    Timespan & operator=(const Timespan& b);
    // arithmetic
    Timespan operator+(const Timespan& b) const;
    Timespan operator-(const Timespan& b) const;
    Timespan operator/(const double b) const;
    Timespan operator*(const double b) const;
    // compound arithmetic
    Timespan & operator+=(const Timespan& b);
    Timespan & operator-=(const Timespan& b);
    Timespan & operator/=(const double b);
    Timespan & operator*=(const double b);
    // comparison
    bool operator==(const Timespan& b) const;
    bool operator!=(const Timespan& b) const;
    bool operator>(const Timespan& b) const;
    bool operator<(const Timespan& b) const;
    bool operator>=(const Timespan& b) const;
    bool operator<=(const Timespan& b) const;

    friend double& operator +=(double& a, const Timespan& b);
    friend double& operator -=(double& a, const Timespan& b);

private:
    /*
     * stores value in minutes
     */
    double time_span_;
};

#endif

