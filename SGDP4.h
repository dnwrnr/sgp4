#ifndef SGDP4_H_
#define SGDP4_H_

#include "Tle.h"

class SGDP4
{
public:
  SGDP4(void);
  virtual ~SGDP4(void);

  void Initialize(const Tle& tle);

      struct TleData {
        double bstar;
        double eo;
        double omega;
        double xincl;
        double xmo;
        double xno;
        double xnodeo;
        Julian epoch;
    };
    
private:
  bool first_run_;

  struct TleData tle_data_;
};

#endif

