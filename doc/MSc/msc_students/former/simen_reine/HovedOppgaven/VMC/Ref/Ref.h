#ifndef Ref_IS_INCLUDED
#define Ref_IS_INCLUDED

#ifdef _DEBUG_
#define REF(TYPE) RefFund<TYPE>
#endif

#ifndef _DEBUG_
#define REF(TYPE) TYPE*
#endif

// ****************************************************************
// *                             REF                              *
// ****************************************************************
template <class type>
class Ref {

 protected:
  type* item;

 public:
  Ref() {}
  Ref(type& _item)            {item = &_item;}
  void attach(type& _item)    {item = &_item;}
  void operator=(type& _item) {item = &_item;}
  void operator=(type* _item) {item =  _item;}
  inline type& operator()()   {return *item;}
  inline type* getPtr()       {return item;}
};

// ****************************************************************
// *                           REFFUND                            *
// ****************************************************************
template <class type>
class RefFund : public Ref<type> {

 public:
  void operator=(type& _item) {item = &_item;}
  void operator=(type* _item) {item =  _item;}
  inline type& operator()()   {return *item;}
  inline operator type& ()    {return *item;}  

};

#endif
