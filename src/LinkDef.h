#ifdef __CLING__
#pragma link C++ class wood::GenericOut+;
#pragma link C++ class wood::GobbiOut+;
#pragma link C++ class wood::TexNeutOut+;
#pragma link C++ class wood::S800Out+;
#pragma link C++ class std::deque<wood::GobbiOut>+;
#pragma link C++ class std::deque<wood::TexNeutOut>+;
#pragma link C++ class std::deque<wood::S800Out>+;
#pragma link C++ class CKinematics+;
#if defined(rel) && rel == 1
#pragma link C++ class CEinstein+;
#else
#pragma link C++ class CNewton+;
#endif
#pragma link C++ class solution+;
#pragma link C++ class std::vector<solution>+;
#endif
