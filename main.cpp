#include "Julian.h"

#include <stdio.h>

int main() {
    Julian julian_now;
    int year = 0;
    int month = 0;
    double day = 0.0;
    julian_now.GetComponent(year, month, day);
    printf("year: %i\nmonth: %i\nday: %lf\n", year, month, day);
    return 0;
}