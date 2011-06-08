#ifndef TIMESPAN_H_
#define TIMESPAN_H_

class Timespan {
public:
    Timespan(void);
    Timespan(const double time_span);
    Timespan(const Timespan& b);
    virtual ~Timespan(void);

    double GetTotalDays() const;
    double GetTotalHours() const;
    double GetTotalMinutes() const;
    double GetTotalSeconds() const;

    /*
     * overloaded operators
     */
    Timespan & operator=(const Timespan& b);
    Timespan operator+(const Timespan& b) const;
    Timespan operator-(const Timespan& b) const;
    Timespan operator+() const;
    Timespan operator-() const;
    const Timespan & operator+=(const Timespan& b);
    const Timespan & operator-=(const Timespan& b);
    bool operator==(const Timespan& b) const;
    bool operator!=(const Timespan& b) const;
    bool operator>(const Timespan& b) const;
    bool operator<(const Timespan& b) const;
    bool operator>=(const Timespan& b) const;
    bool operator<=(const Timespan& b) const;

private:
    /*
     * stores value in minutes
     */
    double time_span_;
};

#endif

