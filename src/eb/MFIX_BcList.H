#ifndef _MFIX_BCLIST_H_
#define _MFIX_BCLIST_H_

// Boundary conditions list auxiliary structure
//
class BcList {
  public:
    BcList(const int undefined = 0,
           const int ig = 9,
           const int pinf = 10,
           const int pout = 11,
           const int minf = 20,
           const int nsw = 100)
      : und_(undefined)
      , ig_(ig)
      , pinf_(pinf)
      , pout_(pout)
      , minf_(minf)
      , nsw_(nsw)
    {}

    const int get_undefined() const {return und_;};
    const int get_ig() const {return ig_;};
    const int get_pinf() const {return pinf_;};
    const int get_pout() const {return pout_;};
    const int get_minf() const {return minf_;};
    const int get_nsw() const {return nsw_;};

  private:
    const int und_;
    const int ig_;
    const int pinf_;
    const int pout_;
    const int minf_;
    const int nsw_;
};

#endif
